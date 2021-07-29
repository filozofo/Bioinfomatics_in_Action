#题目http://www.bio-info-trainee.com/3415.html
rm(list = ls())
options(stringsAsFactors = F)
Sys.setenv("LANGUAGE"="EN")
setwd("C:/Users/aubot/Documents/code/learn/jm-rexercise")
library(tidyverse)

#1.安装一些R包：
pkg <- c("ALL","CLL", "pasilla", "airway","limma","DESeq2","clusterProfiler","reshape2","ggplot2")

install.packages(pkg,chara)
if (!requireNamespace(pkg, quietly = TRUE)){
      BiocManager::install(pkg)}

#2.了解ExpressionSet对象，比如CLL包里面就有data(sCLLex) ，找到它包含的元素，提取其表达矩阵(使用exprs函数)，查看其大小
#参考：http://www.bio-info-trainee.com/bioconductor_China/software/limma.html
#参考：ExpressionSet对象 https://github.com/bioconductor-china/basic/blob/master/ExpressionSet.md
library(CLL)
data(sCLLex)
sCLLex
dim(sCLLex)
exprSet <- exprs(sCLLex)  %>% 
  as_tibble() %>% 
  mutate(probe_id = row.names(exprs(sCLLex))) 

#3.了解 str,head,help函数，作用于 第二步提取到的表达矩阵
str(exprSet)
head(exprSet)
#4.安装并了解 hgu95av2.db 包,看看 ls("package:hgu95av2.db") 后 显示的那些变量
if (!requireNamespace("hgu95av2.db", quietly = TRUE)){
      BiocManager::install("hgu95av2.db")}
library(hgu95av2.db)
ls("package:hgu95av2.db")

#5.理解 head(toTable(hgu95av2SYMBOL)) 的用法，找到 TP53 基因对应的探针ID
head(toTable(hgu95av2SYMBOL))
probe2symbol  <- toTable(hgu95av2SYMBOL) %>% as_tibble()
probe2symbol %>% filter(symbol=="TP53")
#6.理解探针与基因的对应关系，总共多少个基因，基因最多对应多少个探针，
#是哪些基因，是不是因为这些基因很长，所以在其上面设计多个探针呢？
#参考 https://www.jianshu.com/p/fa5f6ddb1543
#探针和基因的对应关系:
#在用基因芯片中,是通过一系列已知的核酸探针,
#然后和带有荧光标记的核酸序列互补配对,通过测定荧光强度最强的位置,
#获取一组互补的探针序列.后期的分析,通过探针组的结果,去分析对应的基因.一一对应.
#总共多少个基因，基因最多对应多少个探针，是哪些基因
length(unique(probe2symbol$symbol))  #查看长度,也即是个数
# [1] 8776   #一共有8584个基因symbol
probe2symbol %>% count(symbol,sort = TRUE) %>% head()
#统计频数，查看前6

#7第二步提取到的表达矩阵是12625个探针在22个样本的表达量矩阵，找到那些不在 hgu95av2.db 包收录的对应着SYMBOL的探针。

anti_join(exprSet,probe2symbol,by = "probe_id") %>% dim()

full_join(exprSet,probe2symbol,by = "probe_id") %>% 
    filter(is.na(symbol)) %>% 
    count(symbol)

#8.过滤表达矩阵，删除那1165个没有对应基因名字的探针。
#9.整合表达矩阵，多个探针对应一个基因的情况下，只保留在所有样本里面平均表达量最大的那个探针。
#提示，理解 tapply,by,aggregate,split 函数 , 首先对每个基因找到最大表达量的探针。
#然后根据得到探针去过滤原始表达矩阵
#10把过滤后的表达矩阵更改行名为基因的symbol，因为这个时候探针和基因是一对一关系了。

transid <- function(probe2symbol,exprSet,method="mean"){
  library(dplyr)
  ex <- exprSet %>% 
  dplyr::inner_join(probe2symbol,by="probe_id") %>% #合并探针的信息
  dplyr::select(-probe_id) %>% #去掉多余信息
  dplyr::select(symbol, everything()) %>% #重新排列，
  dplyr::mutate(ref=apply(across(where(is.numeric)),1,method)) %>% #求出平均数(这边的.真的是画龙点睛)
  dplyr::arrange(desc(ref)) %>% #把表达量的平均值按从大到小排序
  dplyr::select(-ref) %>%  #反向选择去除rowMean这一列
  dplyr::distinct(symbol,.keep_all = T)
  return(ex)
}

ex <- transid(probe2symbol,exprSet)

#11.对第10步得到的表达矩阵进行探索，先画第一个样本的所有基因的表达量的boxplot,hist,density ， 然后画所有样本的 这些图
#参考：http://bio-info-trainee.com/tmp/basic_visualization_for_expression_matrix.html
#理解ggplot2的绘图语法，数据和图形元素的映射关系

pd <- pData(sCLLex) %>% rownames_to_column("sample")

ex1 <- ex %>% 
    pivot_longer(cols = 2:23,names_to = "sample", values_to = "value") %>% 
    full_join(pd,by = "sample") %>% 
    group_by(Disease)
#boxplot
ggplot(ex1,mapping = aes(Disease,value,fill=Disease)) + geom_boxplot()
#histogram
ggplot(ex1,mapping = aes(value,fill=Disease)) + geom_histogram()
#density
ggplot(ex1,aes(value,col=Disease))+geom_density() 
ggplot(ex1,aes(value,col=Disease))+geom_density()+facet_wrap(~sample, nrow = 4)
#12.#理解统计学指标mean,median,max,min,sd,var,mad并计算出每个基因在所有样本的这些统计学指标，最后按照mad值排序，取top 50 mad值的基因，得到列表。
#注意：这个题目出的并不合规，请仔细看。

top50 <- ex %>% mutate(ref=apply(across(where(is.numeric)),1,mad)) %>% slice_max(ref,n=50)

#13.#根据第12步骤得到top 50 mad值的基因列表来取表达矩阵的子集，并且热图可视化子表达矩阵。试试看其它5种热图的包的不同效果。
top50_mat <- top50 %>% column_to_rownames("symbol") %>% as.matrix()
library(pheatmap)
pheatmap(top50_mat,scale = "row")
library(gplots)
install.packages("heatmap")
heatmap.2(top50_mat,scale = "row")
heatmap(top50_mat,scale = "row")
#14.#取不同统计学指标mean,median,max,mean,sd,var,mad的各top50基因列表，使用UpSetR包来看他们之间的overlap情况。

select_top <- function(expr,method="mean",nu=50){
  top50 <- expr %>% 
    mutate(ref=apply(across(where(is.numeric)),1,method)) %>% 
    slice_max(ref,n=nu) %>% 
    pull(symbol) %>% 
  return(top50)
}
top50_mean <- select_top(expr=ex)
top50_median <- select_top(expr=ex,method = "median")
top50_max <- select_top(expr=ex,method = "max")
top50_sd <- select_top(expr=ex,method = "sd")
top50_var <- select_top(expr=ex,method = "var")
top50_mad <- select_top(expr=ex,method = "mad")

top50_ups <- list(top50_mean=top50_mean,top50_median=top50_median,
  top50_max=top50_max,top50_sd=top50_sd,top50_var=top50_var,top50_mad=top50_mad)

library(UpSetR)

upset(fromList(top50_ups),
  nsets=20,
  nintersects = 1000, 
  sets = c("top50_mean","top50_median","top50_max","top50_sd","top50_var","top50_mad"),
  point.size = 3,
  line.size = 1,
  number.angles =0, 
  text.scale = c(1.5, 1.2, 1.2, 1, 1.5,1),
  order.by="freq",
  matrix.color="green",
  main.bar.color = 'yellow',
  mainbar.y.label = 'Intersection Size', 
  sets.bar.color="red")

#15.#在第二步的基础上面提取CLL包里面的data(sCLLex) 数据对象的样本的表型数据。
#参考11题
#16.#对所有样本的表达矩阵进行聚类并且绘图，然后添加样本的临床表型数据信息(更改样本名)
exps <- transid(probe2symbol,exprSet,method = "mean") %>% 
    column_to_rownames("symbol") %>%
    as.matrix()
out.dist=dist(t(exps),method='euclidean')
out.hclust=hclust(out.dist,method='complete')
plot(out.hclust)
#17.#对所有样本的表达矩阵进行PCA分析并且绘图，同样要添加表型信息。
pc <- prcomp(t(exps),scale=TRUE)
pcx=data.frame(pc$x)
group_list=sCLLex$Disease
pcr=cbind(samples=rownames(pcx),group_list, pcx) 
p=ggplot(pcr, aes(PC1, PC2))+geom_point(size=5, aes(color=group_list)) +
  geom_text(aes(label=samples),hjust=-0.1, vjust=-0.3)
print(p)
#18.#根据表达矩阵及样本分组信息进行批量T检验，得到检验结果表格
#按照Disease分为两组检验
dat = exprSet
group_list <- pd$Disease
group_list=as.factor(group_list)
group1 = which(group_list == levels(group_list)[1])
group2 = which(group_list == levels(group_list)[2])
dat1 = dat[, group1]
dat2 = dat[, group2]
dat = cbind(dat1, dat2)
tb <- exprSet %>% dplyr::select(-probe_id)
pvals = apply(tb, 1, function(x){
  t.test(as.numeric(x)~group_list)$p.value
})
p.adj = p.adjust(pvals, method = "BH")
avg_1 = rowMeans(dat1)
avg_2 = rowMeans(dat2)
log2FC = avg_2-avg_1
DEG_t.test = cbind(avg_1, avg_2, log2FC, pvals, p.adj)
DEG_t.test=DEG_t.test[order(DEG_t.test[,4]),]
DEG_t.test=as.data.frame(DEG_t.test)
head(DEG_t.test)

#19.#使用limma包对表达矩阵及样本分组信息进行差异分析，得到差异分析表格，重点看logFC和P值，画个火山图(就是logFC和-log10(P值)的散点图。)。
# DEG by limma 
suppressMessages(library(limma)) 
design <- model.matrix(~0+group_list)
colnames(design)=levels(group_list)
rownames(design)=pd$sample
design

## 下面的 contrast.matrix 矩阵非常重要，制定了谁比谁这个规则
contrast.matrix<-makeContrasts(paste0(unique(group_list),collapse = "-"),levels = design)
contrast.matrix 
##这个矩阵声明，我们要把progres.组跟stable进行差异分析比较
##step1
fit <- lmFit(tb,design)
##step2
fit2 <- contrasts.fit(fit, contrast.matrix) ##这一步很重要，大家可以自行看看效果
fit2 <- eBayes(fit2)  ## default no trend !!!
##eBayes() with trend=TRUE
##step3
DEG  <-  topTable(fit2, coef=1, n=Inf) %>% drop_na()

## volcano plot
logFC_cutoff <- with(DEG,mean(abs(logFC)) + 2*sd(abs(logFC)) )
DEG$change = as.factor(ifelse(DEG$P.Value < 0.05 & abs(DEG$logFC) > logFC_cutoff,
                              ifelse(DEG$logFC > logFC_cutoff ,'UP','DOWN'),'NS')
)
this_tile <- paste0('Cutoff for logFC is ',
  round(logFC_cutoff,3),
  '\nThe number of up gene is ',
  nrow(DEG[DEG$change =='UP',]),
  '\nThe number of down gene is ',
  nrow(DEG[DEG$change =='DOWN',]))
this_tile
head(DEG)
gp  <-  ggplot(data=DEG, aes(x=logFC, y=-log10(P.Value), color=change))

gp + geom_point(alpha=0.4, size=1.75) +
  theme_set(theme_set(theme_bw(base_size=20)))+
  xlab("log2 fold change") + 
  ylab("-log10 p-value") +
  ggtitle(this_tile) + 
  theme(plot.title = element_text(size=15,hjust = 0.5))+
  scale_colour_manual(values = c('blue','black','red'))  ## corresponding to the levels(res$change)

#20.#对T检验结果的P值和limma包差异分析的P值画散点图。看看哪些基因相差很大。
### different P values 
head(DEG)
head(DEG_t.test)
DEG_t.test=DEG_t.test[rownames(DEG),]

plot(DEG_t.test[,3],DEG[,1]) ## 可以看到logFC是相反的
plot(DEG_t.test[,4],DEG[,4]) # 可以看到使用limma包和t.test本身的p值差异尚可接受
plot(-log10(DEG_t.test[,4]),-log10(nrDEG[,4]))

library(ggplot2)
library(ggpubr)
my_comparisons <- pd %>% pull(Disease) %>% unique() %>% unfactor() %>% list()
ex1 %>% filter(symbol=="DLEU1") %>% 
  ggboxplot(x = "Disease", y = "value",
    color = "Disease",
    add = "jitter") +
    stat_compare_means(comparisons = my_comparisons, method = "t.test")

## heatmap 
library(pheatmap)
pheatmap(exps)

