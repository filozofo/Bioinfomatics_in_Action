#题目 http://www.bio-info-trainee.com/3793.html

rm(list = ls())
options(stringsAsFactors = F)
Sys.setenv("LANGUAGE"="EN")
library(tidyverse)
setwd("C:/Users/aubot/Documents/code/learn/jm-rexercise")
a <- read_delim('SraRunTable.txt',delim ="\t")

b=read_csv('sample.csv')

colnames(a)
colnames(b)

d <- full_join(a,b,by=c('Sample_Name'='Accession'))

e <- d %>% select(MBases,Title)

save(e,file = 'input.Rdata')

load(file = 'input.Rdata')

plate <- e %>% pull(2) %>% map_chr(~strsplit(.x,'_')[[1]][3])

f <- e  %>% mutate(grp=plate)

boxplot(pull(e,1)~plate)

t.test(pull(e,1)~plate)

ggplot(f,mapping = aes(grp,MBases))+geom_boxplot()

library(ggpubr)

ggboxplot(f,x="grp",y="MBases") + stat_compare_means(method = 't.test')
