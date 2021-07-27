#题目 http://www.bio-info-trainee.com/3750.html
rm(list = ls())
options(stringsAsFactors = F)
Sys.setenv("LANGUAGE"="EN")
setwd("C:/Users/aubot/Documents/code/learn/jm-rexercise")
library(tidyverse)



# 作业 1
#根据R包org.Hs.eg.db找到下面ensembl 基因ID 对应的基因名(symbol)

## 基因ID转换使用select方法
library(org.Hs.eg.db)

x <- org.Hs.eg.db

ids <- c("ENSG00000000003.13", "ENSG00000000005.5", "ENSG00000000419.11", 
"ENSG00000000457.12", "ENSG00000000460.15", "ENSG00000000938.11")


ids <- map_chr(ids,~strsplit(.x,"[.]")[[1]][1])

select(x,keys = ids, columns = c('SYMBOL','GENENAME'), keytype = 'ENSEMBL')

# 作业 2
#根据R包hgu133a.db找到下面探针对应的基因名(symbol)

library(hgu133a.db)

ids=toTable(hgu133aSYMBOL)

prob <- c("1053_at","117_at","121_at","1255_g_at","1316_at","1320_at",
"1405_i_at","1431_at","1438_at","1487_at","1494_f_at","1598_g_at",
"160020_at","1729_at","177_at") 

prob <- tibble(prob)

left_join(prob,ids,by=c("prob"="probe_id"))

#作业 3
#找到R包CLL内置的数据集的表达矩阵里面的TP53基因的表达量，并且绘制在 progres.-stable分组的boxplot图

suppressPackageStartupMessages(library(CLL))

#数据集
data(sCLLex)

#获得表达矩阵

exprSet <- exprs(sCLLex)  %>% as_tibble() %>% mutate(probe_id = row.names(exprs(sCLLex))) 


#探针转换为基因名
library(hgu95av2.db)

probe2entrezID=toTable(hgu95av2ENTREZID)
probe2symbol <- toTable(hgu95av2SYMBOL)  %>% as_tibble()
probe2genename=toTable(hgu95av2GENENAME)

transid <- function(probe2symbol,exprSet,method){
  library(tidyverse)
  ex <- exprSet %>% 
  dplyr::inner_join(probe2symbol,by="probe_id") %>% #合并探针的信息
  dplyr::select(-probe_id) %>% #去掉多余信息
  dplyr::select(symbol, everything()) %>% #重新排列，
  dplyr::mutate(ref=apply(across(is.numeric),1,method)) %>% #求出平均数(这边的.真的是画龙点睛)
  dplyr::arrange(desc(ref)) %>% #把表达量的平均值按从大到小排序
  dplyr::select(-ref) %>%  #反向选择去除rowMean这一列
  dplyr::distinct(symbol,.keep_all = T)
  return(ex)
}

ex <- transid(probe2symbol,exprSet,"mean")

ex2 <- ex %>% filter(symbol=="TP53")

data(disease)

ds <- disease  %>% group_by(Disease)


ex3  <- ex2 %>% rename_with(~substr(.x,1,5))

ex4 <- ex3 %>% pivot_longer(cols = 2:23,names_to = "sample", values_to = "value")

ex5 <- full_join(ex4,ds,by = c("sample"="SampleID")) %>% group_by(Disease)

boxplot(ex5$value~ex5$Disease)


#作业 4
#找到BRCA1基因在TCGA数据库的乳腺癌数据集(Breast Invasive Carcinoma (TCGA, PanCancer Atlas))的表达情况

#提示：使用http://www.cbioportal.org/index.do 定位数据集：http://www.cbioportal.org/datasets

dt <- read_delim("plot.txt",delim = "\t") 
dt %>% colnames()

dt <- dt %>% rename("id","subtype","expression","mutant","CNA")
names(dt) <- c("id","subtype","expression","mutant","CNA")
library(ggpubr)
p <- ggboxplot(data =dt, x = 'subtype',  y = 'expression', 
          color = "subtype",
          add = "jitter", shape = "subtype",
          palette = "lancet" )

p + theme(axis.text.x = element_text(angle = 45, hjust = 1))

#作业 5
#找到TP53基因在TCGA数据库的乳腺癌数据集的表达量分组看其是否影响生存
#提示使用：http://www.oncolnc.org/

#作业6
#下载数据集GSE17215的表达矩阵并且提取下面的基因画热图

geodown <- function(geo){
    library(GEOquery)
    gd <- getGEO(geo,GSEMatrix=TRUE)
    eSet=gd[[1]]
    return(eSet)}

eSet <- geodown('GSE17215')

probes_expr <- exprs(eSet)

head(probes_expr[,1:4])

boxplot(probes_expr,las=2)

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

boxplot(probes_expr,las=2)

genes <- read_lines("06.txt") %>% strsplit(" ") 
  
genes <- matrix(data=unlist(genes), ncol=1) %>% 
  as_tibble_col(column_name="symbol")

#获得探针和symbol两列对应信息
getprobe2symb <- function(eSet){
  index <- eSet@annotation
  nm <- read_delim("platformMap.txt",delim = ", ") %>% 
    filter(gpl==index) %>% pull(4)
  pkg <- paste0(nm,".db")
  symb <- paste0(nm,"SYMBOL")
  if (!requireNamespace(pkg, quietly = TRUE)){
      BiocManager::install(pkg)
      }
  library(pkg, character.only = TRUE)
  probe2symbol <- toTable(get(symb))
  return(probe2symbol)
}

probe2symbol <- getprobe2symb(eSet)

exprSet <- probes_expr %>% as_tibble() %>% mutate(probe_id = row.names(probes_expr)) 

#探针转换成symbol,获得表达矩阵用于绘制热图
exps <- transid(probe2symbol,exprSet,"mean") %>% 
    inner_join(genes,by="symbol") %>% 
    drop_na()  %>% 
    column_to_rownames("symbol") %>% 
    as.matrix()

library(pheatmap)

pheatmap(exps)

#作业7
#下载数据集GSE24673的表达矩阵计算样本的相关性并且绘制热图，需要标记上样本分组信息
#题目 http://www.bio-info-trainee.com/3750.html

eSet <- geodown('GSE24673')

exps <- transid(probe2symbol,exprSet,method = "mean") %>% 
    column_to_rownames("symbol") %>%
    as.matrix()
grp <- pData(eSet) %>% dplyr::select(source_name_ch1)

library(pheatmap)

pheatmap(exps,scale = 'row', show_rownames=F,annotation_col = grp)

#作业8
#找到 GPL6244 platform of Affymetrix Human Gene 1.0 ST Array 对应的R的bioconductor注释包，并且安装它！
index <- "GPL6244"
nm <- read_delim("platformMap.txt",delim = ", ") %>% filter(gpl==index) %>% pull(4)
pkg <- paste0(nm,".db")
if (!requireNamespace(pkg, quietly = TRUE)){
    BiocManager::install(pkg)
    }

#作业9
#下载数据集GSE42872的表达矩阵，并且分别挑选出 所有样本的(平均表达量/sd/mad/)最大的探针，并且找到它们对应的基因。

eSet <- geodown("GSE42872")

probes_expr <- exprs(eSet)

probe2symbol <- getprobe2symb(eSet)

exprs_mean <- transid(probe2symbol,exprSet,"mean")

exprs_sd <- transid(probe2symbol,exprSet,"sd")

exprs_mad <- transid(probe2symbol,exprSet,"mad")

#作业10
#下载数据集GSE42872的表达矩阵，并且根据分组使用limma做差异分析，得到差异结果矩阵

counts <- exprs_mean %>% column_to_rownames("symbol")
library(limma)
library(edgeR)
dge <- DGEList(counts = counts)
dge <- calcNormFactors(dge)
logCPM <- cpm(dge, log=F, prior.count=3)
#设置分组信息
group_list <- factor(c(rep("control",3), rep("treatment",3)))
design <- model.matrix(~group_list)
colnames(design) <- levels(group_list)
rownames(design) <- colnames(counts)
#常规的差异分析
fit <- lmFit(logCPM, design)
fit <- eBayes(fit, trend=TRUE)
output <- topTable(fit, coef=2,n=Inf)