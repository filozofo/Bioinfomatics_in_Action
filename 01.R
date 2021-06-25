rm(list = ls())
options(stringsAsFactors = F)
Sys.setenv("LANGUAGE"="EN")
setwd("C:\\Users\\aubot\\Documents\\code\\RWork\\LYM")
library(DESeq2)
library(tidyverse)
tb <- read_delim("count_matrix.txt",delim="\t")  %>% column_to_rownames(var="gene")
### 选取表达量在所有六个样本总和大于25的基因。
tb %>% filter(if_all(1:6,~.x > 25))  %>% write_csv("express25.csv")
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

