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
# 1. Read DE tables ####
l35l39 <-read.table("resultadosDEs/L35_vs_L39_DE_annotated.txt",sep="\t",header=T)
l35mlib <-read.table("resultadosDEs/L35_vs_MLIB_DE_annotated.txt",sep="\t",header=T)
l35mmix <-read.table("resultadosDEs/L35_vs_MMIX_DE_annotated.txt",sep="\t",header=T)
mlibmmix <-read.table("resultadosDEs/MLIB_vs_MMIX_DE_annotated.txt",sep="\t",header=T)
l35l39 <- l35l39 %>% select(1:7,20)
l35mlib <- l35mlib %>% select(1:7,20)
l35mmix <- l35mmix %>% select(1:7,20)
mlibmmix <-mlibmmix %>% select(1:7,20)

# 2. Clusterprofiler GO enrichment ####
## L35 vs L39 ####
upl35 <- l35l39 %>%
  filter(logFC>2)
upl39 <- l35l39 %>%
  filter(logFC < -2)

l35u <- enrichGO(gene          = upl35$transcript_id,
                #universe      = gene2go$GID,
                OrgDb         = org.Mridens.eg.db,
                ont           = "all",
                pAdjustMethod = "bonferroni",
                pvalueCutoff  = 0.05,
                qvalueCutoff  = 0.05,
                readable      = F,
                keyType = "GID")

l35us <- pairwise_termsim(l35u)
p1<-emapplot(l35us,
             cex_category=0.8,
             cex_line = 0.5,
             label_style = "ggforce") 

l39u <- enrichGO(gene          = upl39$transcript_id,
                 #universe      = gene2go$GID,
                 OrgDb         = org.Mridens.eg.db,
                 ont           = "BP",
                 pAdjustMethod = "bonferroni",
                 pvalueCutoff  = 0.05,
                 qvalueCutoff  = 0.05,
                 readable      = F,
                 keyType = "GID")
l39us <- pairwise_termsim(l39u)
p2<-emapplot(l39us,
             cex_category=0.8,
             cex_line = 0.5,
             label_style = "ggforce") 

l35u2 <- enrichGO(gene          = upl35$transcript_id,
                 #universe      = gene2go$GID,
                 OrgDb         = org.Mridens.eg.db,
                 ont           = "MF",
                 pAdjustMethod = "bonferroni",
                 pvalueCutoff  = 0.05,
                 qvalueCutoff  = 0.05,
                 readable      = F,
                 keyType = "GID")
l35us2 <- pairwise_termsim(l35u2)
p3<-emapplot(l35us2,
             cex_category=0.8,
             cex_line = 0.5,
             label_style = "ggforce") 

l39u2 <- enrichGO(gene          = upl39$transcript_id,
                 #universe      = gene2go$GID,
                 OrgDb         = org.Mridens.eg.db,
                 ont           = "BP",
                 pAdjustMethod = "bonferroni",
                 pvalueCutoff  = 0.05,
                 qvalueCutoff  = 0.05,
                 readable      = F,
                 keyType = "GID")
l39us <- pairwise_termsim(l39u)
p2<-emapplot(l39us,
             cex_category=0.8,
             cex_line = 0.5,
             label_style = "ggforce") 

## L35 vs MLIB ####
l35mlib <- enrichGO(gene          = l35mlib$transcript_id,
                  #universe      = gene2go$GID,
                  OrgDb         = org.Mridens.eg.db,
                  ont           = "all",
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 0.05,
                  qvalueCutoff  = 0.05,
                  readable      = F,
                  keyType = "GID")

l35mlibs <- pairwise_termsim(l35mlib)
p1<-emapplot(l35mlibs,
             cex_category=0.8,
             cex_line = 0.5,
             label_style = "ggforce") 

## L35 vs MMIX ####
l35mmix <- enrichGO(gene          = l35mmix$transcript_id,
                    #universe      = gene2go$GID,
                    OrgDb         = org.Mridens.eg.db,
                    ont           = "all",
                    pAdjustMethod = "BH",
                    pvalueCutoff  = 0.05,
                    qvalueCutoff  = 0.05,
                    readable      = F,
                    keyType = "GID")

l35mmixs <- pairwise_termsim(l35mmix)
p1<-emapplot(l35mmixs,
             cex_category=0.8,
             cex_line = 0.5,
             label_style = "ggforce") 

## MLIB vs MMIX ####
mm <- enrichGO(gene          = mlibmmix$transcript_id,
                    #universe      = gene2go$GID,
                    OrgDb         = org.Mridens.eg.db,
                    ont           = "all",
                    pAdjustMethod = "BH",
                    pvalueCutoff  = 0.05,
                    qvalueCutoff  = 0.05,
                    readable      = F,
                    keyType = "GID")

mms <- pairwise_termsim(mm)
p1<-emapplot(l35mmixs,
             cex_category=0.8,
             cex_line = 0.5,
             label_style = "ggforce") 

# 3. Net GO-enrich ####
## L35 vs L39 ####
l3539 <- enrichGO(gene          = l35l39$transcript_id,
                  #universe      = gene2go$GID,
                  OrgDb         = org.Mridens.eg.db,
                  ont           = "all",
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 0.05,
                  qvalueCutoff  = 0.05,
                  readable      = F,
                  keyType = "GID")
fc<-l35l39$logFC
names(fc) <- l35l39$transcript_id
cnetplot(l3539, 
         categorySize="pvalue", 
         showCategory = 20, 
         foldChange=fc, 
         vertex.label.font=6,layout="circle")