# Generated from the maintained revision.R; do not edit this distribution separately.
# The numeric model contains no patient-level expression or original fitted objects.
mira_default_model_path <- function() system.file("extdata", "MIRA_final_model.rds", package="MIRA", mustWork=TRUE)
mira_compact_model <- function(full, source_sha256=NA_character_) { if(!identical(full$format,"MIRA_final_numeric_v1")) stop("Use a compact frozen MIRA model."); full }
mira_validate_model <- function (model) 
{
    if (is.character(model) && length(model) == 1L) 
        model <- readRDS(model)
    model <- mira_compact_model(model)
    if (!identical(model$format, "MIRA_final_numeric_v1") || !identical(model$classes, c("C1a", "C1b", "C2")) || 
        !identical(unname(model$thresholds), c(0.5, 0.5)) || !identical(model$alpha, 0.5) || !identical(model$lambda, 
        0.003)) 
        stop("Unsupported final model specification.")
    valid_names <- function(v, n) is.character(v) && length(v) == n && !anyNA(v) && all(nzchar(v)) && !anyDuplicated(v)
    if (!valid_names(model$gene_universe, 18686L) || !valid_names(model$mhc_genes, 8L) || !valid_names(model$b_genes, 
        8L) || !all(c(model$mhc_genes, model$b_genes) %in% model$gene_universe) || !identical(model$first_features, 
        paste0("MHC_", model$mhc_genes)) || !identical(model$second_features, "rankMCP_B lineage")) 
        stop("Invalid frozen gene or feature definitions.")
    for (step in c("first", "second")) {
        f <- model[[paste0(step, "_features")]]
        b <- model[[paste0(step, "_coefficients")]]
        s <- model[[paste0(step, "_scaler")]]
        if (!identical(names(b), c("(Intercept)", f)) || any(!is.finite(b)) || !identical(names(s$mean), f) || 
            !identical(names(s$sd), f) || any(!is.finite(s$mean)) || any(!is.finite(s$sd)) || any(s$sd <= 0)) 
            stop("Invalid frozen coefficients or training transformations.")
    }
    model
}

mira_batch_ids <- function (batch_id, ids) 
{
    if (is.null(batch_id) || !length(batch_id)) 
        stop("Explicit study/assay batch_id is required.")
    if (length(batch_id) == 1L) 
        batch_id <- rep(as.character(batch_id), length(ids))
    else if (!is.null(names(batch_id))) {
        if (anyDuplicated(names(batch_id)) || !setequal(names(batch_id), ids)) 
            stop("Named batch_id must match the sample IDs exactly.")
        batch_id <- batch_id[ids]
    }
    batch_id <- as.character(batch_id)
    if (length(batch_id) != length(ids) || anyNA(batch_id) || any(!nzchar(trimws(batch_id)))) 
        stop("Invalid batch_id.")
    if (any(table(batch_id) < 2L)) 
        stop("Each batch requires at least two samples for calibration; a validated minimum effective batch size has not been established.")
    batch_id
}

mira_two_step_labels <- function (p_C1, p_C1a_given_C1) 
{
    if (length(p_C1) != length(p_C1a_given_C1) || any(!is.finite(c(p_C1, p_C1a_given_C1))) || any(c(p_C1, p_C1a_given_C1) < 
        0 | c(p_C1, p_C1a_given_C1) > 1)) 
        stop("Invalid component probabilities.")
    ifelse(p_C1 < 0.5, "C2", ifelse(p_C1a_given_C1 >= 0.5, "C1a", "C1b"))
}

mira_final_panel <- function (expression, model) 
{
    x <- as.matrix(expression)
    if (!is.numeric(x) || length(dim(x)) != 2L || !nrow(x) || !ncol(x)) 
        stop("Expression must be a nonempty numeric genes-by-samples matrix.")
    for (ids in dimnames(x)) if (anyNA(ids) || any(!nzchar(trimws(ids))) || anyDuplicated(ids)) 
        stop("Unique nonempty gene symbols and sample IDs are required.")
    if (is.null(rownames(x)) || is.null(colnames(x))) 
        stop("Gene symbols and sample IDs are required.")
    missing_genes <- setdiff(model$gene_universe, rownames(x))
    if (length(missing_genes)) 
        stop(sprintf("Missing %d/%d fixed rank-background genes (e.g. %s); no silent mapping or imputation.", length(missing_genes), 
            length(model$gene_universe), paste(head(missing_genes, 5), collapse = ", ")))
    x <- x[model$gene_universe, , drop = FALSE]
    if (any(!is.finite(x))) 
        stop("Required expression values must be finite.")
    ranks <- vapply(seq_len(ncol(x)), function(j) (rank(x[, j], ties.method = "average") - 0.5)/nrow(x), numeric(nrow(x)))
    dimnames(ranks) <- dimnames(x)
    panel <- t(rbind(ranks[model$mhc_genes, , drop = FALSE], colMeans(ranks[model$b_genes, , drop = FALSE])))
    colnames(panel) <- c(model$first_features, model$second_features)
    panel
}

predict_MIRA_final_panel <- function (panel, batch_id = NULL, model_path = mira_default_model_path(), model = NULL, decision_rule = c("two_step", 
    "legacy_joint_argmax")) 
{
    decision_rule <- match.arg(decision_rule)
    m <- mira_validate_model(if (is.null(model)) 
        model_path
    else model)
    x <- as.matrix(panel)
    features <- c(m$first_features, m$second_features)
    if (!is.numeric(x) || !nrow(x) || is.null(colnames(x)) || anyDuplicated(colnames(x)) || !all(features %in% 
        colnames(x))) 
        stop("Missing or invalid final rank-panel features.")
    ids <- rownames(x)
    if (is.null(ids) || anyNA(ids) || any(!nzchar(trimws(ids))) || anyDuplicated(ids)) 
        stop("Unique nonempty sample IDs required.")
    x <- x[, features, drop = FALSE]
    if (any(!is.finite(x)) || any(x < 0 | x > 1)) 
        stop("Panel must contain finite ranks in [0,1].")
    batch <- mira_batch_ids(batch_id, ids)
    z <- x
    calibration <- list()
    for (b in unique(batch)) {
        i <- which(batch == b)
        mu <- colMeans(x[i, , drop = FALSE])
        s <- apply(x[i, , drop = FALSE], 2, stats::sd)
        constant <- !is.finite(s) | s < 1e-10
        s[constant] <- 1
        z[i, ] <- sweep(sweep(x[i, , drop = FALSE], 2, mu, "-"), 2, s, "/")
        calibration[[b]] <- list(mean = mu, sd = s, n = length(i), constant = constant)
    }
    probability <- function(step) {
        f <- m[[paste0(step, "_features")]]
        s <- m[[paste0(step, "_scaler")]]
        coef <- m[[paste0(step, "_coefficients")]]
        a <- sweep(sweep(z[, f, drop = FALSE], 2, s$mean, "-"), 2, s$sd, "/")
        as.numeric(stats::plogis(drop(a %*% coef[f]) + coef["(Intercept)"]))
    }
    g <- probability("first")
    b <- probability("second")
    joint <- cbind(C1a = g * b, C1b = g * (1 - b), C2 = 1 - g)
    call <- if (decision_rule == "two_step") 
        mira_two_step_labels(g, b)
    else m$classes[max.col(joint, ties.method = "first")]
    ans <- data.frame(sample_id = ids, predicted_label = call, p_C1 = g, p_C1a_given_C1 = b, p_C1a = joint[, 1], 
        p_C1b = joint[, 2], p_C2 = joint[, 3], decision_rule = decision_rule, target_batch_id = batch, target_batch_n = as.integer(table(batch)[batch]), 
        model_id = m$model_id, prediction_mode = "fixed_coefficients_with_unlabeled_target_batch_calibration", 
        stringsAsFactors = FALSE, row.names = NULL)
    attr(ans, "unlabeled_batch_calibration") <- calibration
    attr(ans, "probability_note") <- paste("Component model probabilities and derived joint scores are not clinically calibrated confidence.", 
        if (decision_rule == "two_step") 
            "Labels use sequential 0.5/0.5 thresholds, not joint argmax."
        else "Explicit historical comparison: labels use joint argmax.")
    ans
}

.mira_predict_core <- function (expression, batch_id = NULL, model_path = mira_default_model_path(), model = NULL, decision_rule = c("two_step", 
    "legacy_joint_argmax")) 
{
    m <- mira_validate_model(if (is.null(model)) 
        model_path
    else model)
    predict_MIRA_final_panel(mira_final_panel(expression, m), batch_id, model = m, decision_rule = match.arg(decision_rule))
}

mira_read_expression <- function (path) 
{
    if (grepl("\\.rds$", path, ignore.case = TRUE)) {
        x <- readRDS(path)
        if (!is.matrix(x) && !is.data.frame(x)) 
            stop("Expression RDS must contain a matrix or data.frame.")
        return(as.matrix(x))
    }
    d <- utils::read.delim(path, check.names = FALSE, stringsAsFactors = FALSE, colClasses = c("character", rep("numeric", 
        max(1L, length(strsplit(readLines(path, n = 1L, warn = FALSE), "\t", fixed = TRUE)[[1]]) - 1L))))
    if (ncol(d) < 2L) 
        stop("TSV needs gene symbols in the first column and expression in all remaining columns.")
    x <- as.matrix(d[, -1, drop = FALSE])
    rownames(x) <- d[[1]]
    x
}

predict_mira <- function(expression, batch_id=NULL, model=NULL, decision_rule=c("two_step","legacy_joint_argmax")) {
  .mira_predict_core(expression, batch_id=batch_id, model=model, decision_rule=match.arg(decision_rule))
}
predict_MIRA <- predict_mira
predict_MIRA_v2 <- predict_mira
predict_MIRA_final <- predict_mira
mira_model <- function() mira_validate_model(mira_default_model_path())
mira_example <- function() readRDS(system.file("extdata", "synthetic_example.rds", package="MIRA", mustWork=TRUE))
read_mira_expression <- function(path) mira_read_expression(path)
