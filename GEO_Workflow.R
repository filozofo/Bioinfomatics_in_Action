###
rm(list=ls())
options(stringsAsFactors = F)
Sys.setenv("LANGUAGE"="EN")
setwd("C:/Users/aubot/Documents/code/project/POKEMON")
library(pacman)
p_load(GEOquery, tidyverse, limma, AnnoProbe, ggsignif,statmod,parallel)
###

gset=AnnoProbe::geoChina("GSE70680")
# check the ExpressionSet
eSet=gset[[1]]

# extract the expression matrix and phenotype data

probes_expr <- exprs(eSet);dim(probes_expr)
head(probes_expr[,1:4])
boxplot(probes_expr,las=2)
#判断是否log2处理
log2if <- function(expr){
  qx <- as.numeric(quantile(expr, c(0, 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
  LogC <- (qx[5] > 100) ||
    (qx[6]-qx[1] > 50 && qx[2] > 0) ||
    (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
  if (LogC) {
    expr[which(expr <= 0)] <- NaN
    expr <- log2(expr)
    }
  return(expr)
}
probes_expr <- log2if(probes_expr)

boxplot(probes_expr,outline=FALSE, notch=T, las=2)
###标准化
probes_expr=normalizeBetweenArrays(probes_expr)
boxplot(probes_expr,outline=T, notch=T, las=2)

## pheno info
phenoDat <- pData(eSet)
#######-------------然后对表达芯片的探针进行基因注释--------------#############
## check GPL and annotate the probes to genes.
gpl=eSet@annotation
checkGPL(gpl)
printGPLInfo(gpl)
probe2gene=idmap(gpl)
head(probe2gene)
#如果没有GPL注释信息，用一下方法
probes_exprs <- probes_expr %>% 
  as_tibble() %>% 
  mutate(probe_id = row.names(probes_expr)) 
#获取探针信息
idmaps <- function(ann_file,skip){
  probe2symbol = read_delim(ann_file,delim="\t",skip=skip) %>% 
    select(ID,"Gene symbol") %>% 
    rename(probe_id=ID,symbol="Gene symbol")
  return(probe2symbol)
}
probe2symbol = idmaps(ann_file="GPL8759.annot",skip=27)

#ID转换
transid <- function(probe2symbol,exprSet,method){
  library(tidyverse)
  ex <- exprSet %>% 
  dplyr::inner_join(probe2symbol,by="probe_id") %>% #合并探针的信息
  dplyr::select(-probe_id) %>% #去掉多余信息
  dplyr::select(symbol, everything()) %>% #重新排列，
  dplyr::mutate(ref=apply(across(where(is.numeric)),1,method)) %>%
  dplyr::arrange(desc(ref)) %>% #把表达量的平均值按从大到小排序
  dplyr::select(-ref) %>%  #反向选择去除rowMean这一列
  dplyr::distinct(symbol,.keep_all = T)
  return(ex)
}

genes_expr = transid(probe2symbol,probes_exprs,"mean") %>% 
  tidyr::drop_na(any_of(symbol)) %>% 
  column_to_rownames("symbol")

write.table(genes_expr,"GSE70680_express.txt",sep="\t")
head(genes_expr)
#走limma的经???2组差异分???
# do DEG
## define the group
group_list=factor(c(rep('CT',2),rep('TR',2)))
table(group_list)
design=model.matrix(~factor(group_list))
design
fit=lmFit(genes_expr,design)
fit=eBayes(fit)
DEG=topTable(fit,coef=2,n=Inf)
head(DEG)
####对差异分析结果进行一些检???
## visualization
need_deg=data.frame(symbols=rownames(DEG), logFC=DEG$logFC, p=DEG$P.Value)
deg_volcano(need_deg,2)

deg_heatmap(DEG,genes_expr,group_list,30)

check_diff_genes('Rfc4',genes_expr,group_list)
#如果你做了GO/KEGG注释后也可以挑选基因集进行可视???
# 假设我这里对hsa03410感兴???
library(KEGGREST)
cg <- KEGGREST::keggGet("Rfc4")[[1]]$GENE
cg=as.character(sapply(cg[seq(2,length(cg),by=2)], function(x) strsplit(x,';')[[1]][1]))
check_diff_genes( cg ,genes_expr,group_list)


