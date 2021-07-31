library(tidyverse)

# 理解 **定性变量**（qualitative variable） 和 **定量变量**(quantitative variable）

# 定量数据的**集中趋势**指标主要是：众数、分位数和平均数

# 定量数据的**离散趋势**指标主要是：极差，方差和标准差，标准分数，相对离散系数（变异系数），偏态系数与峰态系数

## Q1: 载入R中自带的数据集 iris，指出其每列是定性还是定量数据

data(iris)
head(iris)

# 第1:4裂为定量变量，第5列为定性变量

## Q2: 对数据集 iris的所有定量数据列计算集中趋势指标：众数、分位数和平均数

apply(iris[,1:4],2,median)
apply(iris[,1:4],2,quantile)
colMeans(iris[,1:4])

## Q3:对数据集 iris的所有定性数据列计算水平及频次

table(iris[,5])

## Q4:对数据集 iris的所有定量数据列计算离散趋势指标：方差和标准差等

apply(iris[,1:4],2,var)
apply(iris[,1:4],2,sd)

## Q5:计算数据集 iris的前两列变量的相关性，提示cor函数可以选择3种methods

cor(iris[,1],iris[,2])
cor(iris[,1],iris[,2],method = "pearson" )
cor(iris[,1],iris[,2],method = "kendall" )
cor(iris[,1],iris[,2],method = "spearman" )


## Q6:对数据集 iris的所有定量数据列内部zcore标准化，并计算标准化后每列的平均值和标准差

iris_sc <- scale(iris[,1:4],center=T,scale=T)
colMeans(iris_sc[,1:4])
apply(iris_sc[,1:4],2,sd)


## Q7:计算列内部zcore标准化后 iris的前两列变量的相关性

cor(iris_sc[,1],iris_sc[,2])

## Q8: 根据数据集 iris的第五列拆分数据集后重复上面的Q2到Q7问题

ir <- iris
ir[,5] <- factor(ir[,5])
apply(ir[,1:4],2,median)
apply(ir[,1:4],2,quantile)
colMeans(iris[,1:4])

## Q9:载入R中自带的数据集 mtcars，重复上面的Q1到Q7个问题
dx = mtcars[,c(2,8:11)]
dl = mtcars[,c(1,3,4,5,6,7)]
data(mtcars)
apply(dl,2,median)
apply(dl,2,quantile)
colMeans(mtcars)
#定性数据列计算水平及频次
attach(mtcars)
table(cyl)
table(vs)
table(am)
table(gear)
table(carb)
detach(mtcars)
#定量统计
apply(dl,2,var)
apply(dl,2,sd)
cor(dl[,1],dl[,2])

mtcars_sc <- scale(dl,center=T,scale=T)
colMeans(mtcars_sc)
apply(mtcars_sc,2,sd)


## Q10: 载入r包airway并且通过assay函数拿到其表达矩阵后计算每列之间的相关性

library(airway)
data(airway)
RNAseq_expr=assay(airway)
colnames(RNAseq_expr) 
RNAseq_expr[1:4,1:4] 
# RNAseq_expr 是一个数值型矩阵，属于连续性变量，可以探索众数、分位数和平均数 ，极差，方差和标准差等统计学指标
RNAseq_gl=colData(airway)[,3]
table(RNAseq_gl)

t.test(t(RNAseq_expr)~RNAseq_gl)
