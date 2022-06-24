install.packages("org.Mridens.eg.db",repos = NULL, type="source")
library(org.Mridens.eg.db)
library(clusterProfiler)
library(GOSemSim)
library(enrichplot)
 library(GOplot)
library(AnnotationDbi)
library(DOSE)
library(dplyr)
library(AnnotationHub)
library(stringr)
library(ggplot2)
options(ggrepel.max.overlaps = Inf)
set.seed(997766)
# 1. Read DE tables ####
l35l39 <-read.table("resultadosDEs/L35_vs_L39_DE_annotated.txt",sep="\t",header=T)
l35mlib <-read.table("resultadosDEs/L35_vs_MLIB_DE_annotated.txt",sep="\t",header=T)
l35mmix <-read.table("resultadosDEs/L35_vs_MMIX_DE_annotated.txt",sep="\t",header=T)
mlibmmix <-read.table("resultadosDEs/MLIB_vs_MMIX_DE_annotated.txt",sep="\t",header=T)
l35l39 <- l35l39 %>% select(1:7,20)
l35mlib <- l35mlib %>% select(1:7,20)
l35mmix <- l35mmix %>% select(1:7,20)
mlibmmix <-mlibmmix %>% select(1:7,20)

# 2. Functions ####
enrich<-function(x){
  enrichKEGG(gene          = x$transcript_id,
           organism         = "ko",
           pAdjustMethod = "bonferroni",
           qvalueCutoff  = 0.1,
           keyType = "kegg")  
}

# 3. Up/Down separated list ####

## L35 vs L39 ####
upl35 <- l35l39 %>%  filter(logFC > 2)
upl39 <- l35l39 %>%  filter(logFC < -2)

## L35 vs MLIB ####
upl35_2 <- l35mlib %>%  filter(logFC > 2)
upmlib <- l35mlib %>%  filter(logFC < -2)

## L35 vs MMIX ####
upl35_3 <- l35mmix %>%  filter(logFC > 2)
upmmix <- l35mmix %>%  filter(logFC < -2)

## MLIB vs MMIX ####
upmlib_2 <- mlibmmix %>%  filter(logFC > 2)
upmmix_2 <- mlibmmix %>%  filter(logFC < -2)


# 4. Up/Down separated enrichGO BP ####
upl35bp   <-enrich(upl35)
upl39bp   <-enrichBP(upl39)
upl35_2bp <-enrichBP(upl35_2)
upmlibbp  <-enrichBP(upmlib)
upl35_3bp <-enrichBP(upl35_3)
upmmixbp  <-enrichBP(upmmix)
upmlib_2bp  <-enrichBP(upmlib_2)
upmmix_2bp  <-enrichBP(upmmix_2)

# 5. Up/Down separated enrichGO MF ####
upl35mf   <-enrichMF(upl35)
upl39mf   <-enrichMF(upl39)
upl35_2mf <-enrichMF(upl35_2)
upmlibmf  <-enrichMF(upmlib)
upl35_3mf <-enrichMF(upl35_3)
upmmixmf  <-enrichMF(upmmix)
upmlib_2mf  <-enrichMF(upmlib_2)
upmmix_2mf  <-enrichMF(upmmix_2)

# 6. Export data ####
a<-data.frame(upl35bp)
b<-data.frame(upl39bp)
c<-data.frame(upl35_2bp)
d<-data.frame(upmlibbp)
e<-data.frame(upl35_3bp)
f<-data.frame(upmmixbp)
g<-data.frame(upmlib_2bp)
h<-data.frame(upmmix_2bp)
i<-data.frame(upl35mf)
j<-data.frame(upl39mf)
k<-data.frame(upl35_2mf)
l<-data.frame(upmlibmf)
m<-data.frame(upl35_3mf)
n<-data.frame(upmmixmf)
o<-data.frame(upmlib_2mf)
p<-data.frame(upmmix_2mf)

write.table(a,file="enrichGO/upl35bp.csv",sep="\t",row.names=F)
write.table(b,file="enrichGO/upl39bp.csv",sep="\t",row.names=F)
write.table(c,file="enrichGO/upl35_2bp.csv",sep="\t",row.names=F)
write.table(d,file="enrichGO/upmlibbp.csv",sep="\t",row.names=F)
write.table(e,file="enrichGO/upl35_3bp.csv",sep="\t",row.names=F)
write.table(f,file="enrichGO/upmmixbp.csv",sep="\t",row.names=F)
write.table(g,file="enrichGO/upmlib_2bp.csv",sep="\t",row.names=F)
write.table(h,file="enrichGO/upmmix_2bp",sep="\t",row.names=F)
write.table(i,file="enrichGO/upl35mf",sep="\t",row.names=F)
write.table(j,file="enrichGO/upl39mf",sep="\t",row.names=F)
write.table(k,file="enrichGO/upl35_2mf.csv",sep="\t",row.names=F)
write.table(l,file="enrichGO/upmlibmf.csv",sep="\t",row.names=F)
write.table(m,file="enrichGO/upl35_3mf.csv",sep="\t",row.names=F)
write.table(n,file="enrichGO/upmmixmf",sep="\t",row.names=F)
write.table(o,file="enrichGO/upmlib_2mf",sep="\t",row.names=F)
write.table(p,file="enrichGO/upmmix_2mf",sep="\t",row.names=F)

#ojo: cambiar al directorio de trabajo donde estan los csv
bp<-Sys.glob("*bp.csv")
mf<-Sys.glob("*mf.csv")

library(openxlsx)
wb<-createWorkbook()
for (i in 1:length(bp))
{
  test<-read.csv(bp[i],header=T,sep="\t")
  # create a sheet in the workbook
  sheetname<-bp[i]
  addWorksheet(wb, sheetName=sheetname)
  # add the data to the new sheet
  writeData(wb, sheetname, test)
}
saveWorkbook(wb,
             file= "enrichGO_BiologicalProcess.xlsx",
             overwrite = TRUE)


wb<-createWorkbook()
for (i in 1:length(mf))
{
  test<-read.csv(mf[i],header=T,sep="\t")
  # create a sheet in the workbook
  sheetname<-mf[i]
  addWorksheet(wb, sheetName=sheetname)
  # add the data to the new sheet
  writeData(wb, sheetname, test)
}
saveWorkbook(wb,
             file= "enrichGO_MolecularFunction.xlsx",
             overwrite = TRUE)

# 7. Graficos ####

bpg<-function(x){
  x <- pairwise_termsim(x)
  p<-emapplot(x,
          showCategory = 50,
           cex_category=1,
           cex_line = 0.3,
           label_style = "ggforce",
           node_label= "group",
           repel=F,
           min_edge=0.2,
           group_category=T,
           cex_label_group=1.2,
           group_legend=T,
            ) 
  return(p)
  }

listas<-list("upl35bp"=upl35bp,
             "upl39bp"=upl39bp,
             "upl35_2bp"=upl35_2bp,
             "upmlibbp"=upmlibbp,
             "upl35_3bp"=upl35_3bp,
             "upmmixbp"=upmmixbp,
             "upmlib_2bp"=upmlib_2bp,
             "upmmix_2bp"=upmmix_2bp)


for(i in 1:length(listas)){
  pdf(paste0("enrichGO/",names(listas[i]),".pdf"),height=8.5,width=14)
  print(bpg(listas[[i]]))
  dev.off()
}


listasmf<-list("upl35mf"=upl35mf,
             "upl39mf"=upl39mf,
             "upl35_2mf"=upl35_2mf,
             "upmlibmf"=upmlibmf,
             "upl35_3mf"=upl35_3mf,
             "upmmixmf"=upmmixmf,
             "upmlib_2mf"=upmlib_2mf,
             "upmmix_2mf"=upmmix_2mf)

for(i in 1:length(listasmf)){
  pdf(paste0("enrichGO/",names(listasmf[i]),".pdf"),height=8.5,width=14)
  print(bpg(listasmf[[i]]))
  dev.off()
}




