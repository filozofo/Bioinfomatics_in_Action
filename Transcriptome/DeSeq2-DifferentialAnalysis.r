#使用所提供的count_matrix.txt和samples.txt文件，完成下列任务：

#读入两个文件。count_matrix.txt记录了两种小鼠细胞的全转录组基因表达量，共六个样品。sample.txt记录了这六个样品的样品信息。KO为一组，共三个重复。WT为对照组，共三个重复。
#使用读入的文件创建DESeq2分析所使用的DESeqDataSet。
#选取表达量在所有六个样本总和大于25的基因。
#使用R语言DESeq2软件包对这些基因进行基因差异表达分析。
#运用plotDispEst和plotMA两个命令分别绘制Dispersion和MA图，并把图像输出保存为PDF格式。
#绘制差异表达分析后所有基因的P-value的histogram，并把图像输出保存为PDF格式。
#统计共有多少基因差异表达的adjusted P-value (padj)小于0.01 。
#根据差异表达分析所得到的P-value对差异表达分析结果进行排序，并将排序过的结果输出到一个CSV文档内。
#差异表达最显著的5个基因分别是什么？Gapdh基因是否有差异表达？为什么？
#请把上述五个差异表达基因加上Gapdh的相关分析结果从总结果中抽提出来，单独存为一个文件。
#请把用于完成上述所有分析流程的R语句写入一个单独文档。
#请把上述所有的图、结果和分析流程的R语句文档打包到一个zip文件并通过邮件提交，邮箱地址l***@jnu.edu.cn。请于2021年7月9日23：59分前提交，逾期不交算缺考。
#温馨提示：分析过程可参考http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html和Beginner's guide to using the DESeq2 package。
#另外，本考试为开卷，但是严格杜绝互相抄袭，如经查重发现有互相抄袭现象，抄袭者和被抄袭者均按零分处理。

rm(list = ls())
Sys.setenv("LANGUAGE"="EN")
library(pacman)
p_load(DESeq2,tidyverse,data.table)
tb <- fread("count_matrix.txt",delim="\t")  %>% column_to_rownames(var="gene")
### 选取表达量在所有六个样本总和大于25的基因。
tb %>% dplyr::filter(if_all(1:6,~.x > 25))  %>% fwrite("express25.csv")
### 使用R语言DESeq2软件包对这些基因进行基因差异表达分析
coldata <- read_delim("samples.txt",delim="\t")

#第1步，构建 DESeqDataSet 对象
dds <- DESeqDataSetFromMatrix(countData = tb, colData = coldata, design= ~condition)
#第2步，计算差异倍数并获得 p 值
dds1 <- DESeq(dds, fitType = 'mean', minReplicatesForReplace = 7, parallel = T)
res <- results(dds1, contrast = c('condition', 'KO', 'WT'))
### 运用plotDispEst和plotMA两个命令分别绘制Dispersion和MA图，并把图像输出保存为PDF格式。
pdf("Dispersion.pdf",width = 30,height = 30)
pde = plotDispEsts(dds1)
dev.off()
pdf("pMA.pdf",width = 30,height = 30)
pma <- plotMA(dds1)
dev.off()
# 绘制差异表达分析后所有基因的P-value的histogram，并把图像输出保存为PDF格式。
dat <- res %>% as_tibble() %>% mutate(gene=rownames(res))
dat %>% ggplot(aes(x=pvalue)+geom_histogram()

phis = ggplot(dat, aes(x=pvalue)) + geom_histogram() %>% 
ggsave("phist.pdf",phis)

#7. 统计共有多少基因差异表达的adjusted P-value (padj)小于0.01 。
dat  %>% filter(padj<0.01) %>% nrow()
sum(res$padj < 0.01, na.rm=TRUE)
#8. 根据差异表达分析所得到的P-value对差异表达分析结果进行排序，并将排序过的结果输出到一个CSV文档内。
dat  %>% arrange(pvalue) %>% write_csv("expressarranged.csv")
#9. 差异表达最显著的5个基因分别是什么？Gapdh基因是否有差异表达？为什么？
dat %>% slice_min(order_by=padj,n=5) %>% pull(gene)
dat %>% slice_max(order_by=padj,n=5) %>% pull(gene)
dat %>% filter(str_detect(gene,"Gapdh$"))
