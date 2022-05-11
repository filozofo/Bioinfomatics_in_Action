# https://mubu.com/doc/3L0wkgGUVg
# 基于R的可视化习题30个
options(stringsAsFactors = F)
library(airway)
data(airway)
library(tidyverse)
# 这里需要自行学习bioconductor里面的RangedSummarizedExperiment对象
airway 
RNAseq_expr=assay(airway)
colnames(RNAseq_expr) 
RNAseq_expr[1:4,1:4] 
# RNAseq_expr 是一个数值型矩阵，属于连续性变量，可以探索众数、分位数和平均数 ，极差，方差和标准差等统计学指标
RNAseq_gl=colData(airway)[,3]
table(RNAseq_gl)

boxplot(RNAseq_expr[,2])
density(RNAseq_expr[,2])
