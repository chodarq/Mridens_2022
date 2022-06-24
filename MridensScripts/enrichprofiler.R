library(tidyr)
library(dplyr)
library(GO.db)
library(clusterProfiler)
library(GOSemSim)
library(enrichplot)
library(GOplot)
library(AnnotationDbi)
library(DOSE)
library(AnnotationHub)
library(stringr)
library(ggplot2)

eggann<-read.csv("Mridens.emapper.annotations",sep="\t",header=T,skip =4)
eggann<-head(eggann,-10)
names(eggann)[1]<-"GID"
lista<-read.csv("lista.isoformas.csv",sep="\t",header=T)
names(lista)[1]<-"GID"

anno<-inner_join(lista,eggann,by="GID")
write.table(anno,"anotacion.isoforma.mas.larga.csv",sep = "\t")

annoGO<-anno[c(1,2,11)]
columnGOs<-separate_rows(annoGO,GOs,sep=",")
listago<-columnGOs %>% 
  filter(!grepl("-", GOs))

term2gene<-listago[c(3,2)]
names(term2gene)<-c("GOID","GenID")
term2name<-select(GO.db,gosunicos,c("TERM"))
gosunicos<-unique(term2gene$GOs)
GO <- as.list(GOTERM)

# 1. Read DE tables ####
l35l39 <-read.table("resultadosDEs/L35_vs_L39_DE_annotated.txt",sep="\t",header=T)
l35mlib <-read.table("resultadosDEs/L35_vs_MLIB_DE_annotated.txt",sep="\t",header=T)
l35mmix <-read.table("resultadosDEs/L35_vs_MMIX_DE_annotated.txt",sep="\t",header=T)
mlibmmix <-read.table("resultadosDEs/MLIB_vs_MMIX_DE_annotated.txt",sep="\t",header=T)
l35l39 <- l35l39 %>% select(1:7,20)
l35mlib <- l35mlib %>% select(1:7,20)
l35mmix <- l35mmix %>% select(1:7,20)
mlibmmix <-mlibmmix %>% select(1:7,20)



