# MIRA 使用说明

安装后直接调用冻结模型，不需要另外下载参数或安装 `glmnet`：

```r
install.packages("remotes") # 已安装则跳过
remotes::install_github("xlucpu/MIRA")
library(MIRA)
example <- mira_example()
prediction <- predict_mira(example$expression, batch_id = example$batch_id)
prediction[, c("sample_id", "predicted_label", "p_C1", "p_C1a_given_C1")]
```

示例为人工合成数据，预期输出在 `example$expected`。自己的数据可以用 `read_mira_expression("expression.tsv")` 读取；首列是基因符号，其余列是样本表达。

输入必须包含固定的全部 **18,686 个背景基因**，不能只给16个marker，也不能给MCP-counter丰度。第一步使用8个MHC-II基因各自的表达秩，第二步使用8个B-lineage基因秩的均值。`mira_model()`可查看基因、冻结系数和训练变换。

最终规则为：第一步 `P(C1)<0.5` 判为C2；否则第二步 `P(C1a|C1)>=0.5` 判为C1a，其余为C1b。不要再对joint scores取最大值来替换标签。`predict_MIRA()`、`predict_MIRA_v2()`和`predict_MIRA_final()`与`predict_mira()`使用相同默认规则。

必须输入按study/assay预先定义的完整适用批次。程序在不使用标签或结局的情况下计算批次均值和样本标准差，然后应用固定训练变换和系数。第二步也使用完整批次进行校准，包括之后被判为C2的样本，不应仅在预测C1样本中重新校准。

批次构成变化可能改变预测。每批至少2例只是计算下限，并非已验证的有效最小样本量；独立单患者检测尚未验证。多个批次可提供与样本列对应的batch向量，具名向量按样本ID匹配。

缺背景基因、重复ID或非有限表达值均会报错。包不负责counts转TPM、基因名映射和缺失值填补。模型概率不是临床置信度或疗效概率。

包内模型是最终合并拟合模型。S11的539例预测来自三个留一队列模型，已经另行核验；不能把最终模型的训练样本预测当作S11留出验证。包内不包含患者表达、训练标签或原始拟合对象。维护代码统一在主`revision.R`，此包由它自动生成。详细说明见英文README和`?predict_mira`。

