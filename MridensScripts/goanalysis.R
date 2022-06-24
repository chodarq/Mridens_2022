install.packages("org.Mridens.eg.db",repos = NULL, type="source")
library(org.Mridens.eg.db)
library(clusterProfiler)
library(GOSemSim)
library(enrichplot)
library(GOplot)
library(AnnotationDbi)
library(DOSE)

# 1. Read DE tables ####
l35l39 <-read.table("resultadosDEs/L35_vs_L39_DE_annotated.txt",sep="\t",header=T)
l35mlib <-read.table("resultadosDEs/L35_vs_MLIB_DE_annotated.txt",sep="\t",header=T)
l35mmix <-read.table("resultadosDEs/L35_vs_MMIX_DE_annotated.txt",sep="\t",header=T)
mlibmmix <-read.table("resultadosDEs/MLIB_vs_MMIX_DE_annotated.txt",sep="\t",header=T)
l35l39 <- l35l39 %>% select(1:7,20)
l35mlib <- l35mlib %>% select(1:7,20)
l35mmix <- l35mmix %>% select(1:7,20)
mlibmmix <-mlibmmix %>% select(1:7,20)

# 2. Analisis semántico de términos GO ####
mrGO<-godata("org.Mridens.eg.db",ont="BP",keytype = "GID")

gs1<-l35l39$transcript_id
gs2<-l35mlib$transcript_id
gs3<-l35mmix$transcript_id
gs4<-mlibmmix$transcript_id
lista1<-list(L35vsL39=gs1,L35vsMLIB=gs2,
            L35vsMMIX=gs3,MLIBvsMMIX=gs4)

sc<-mclusterSim(lista1, semData=mrGO, measure="Wang", combine="BMA")


# 3. clusterprofiler gsea analysis ####
l35mliblist<-l35mlib$logFC
names(l35mliblist)<-l35mlib$transcript_id
l35mliblist = sort(l35mliblist, decreasing = TRUE)

gse3539 <- gseGO(geneList=l35mliblist, 
             ont ="ALL", 
             keyType = "GID", 
             minGSSize = 3, 
             maxGSSize = 800, 
             pvalueCutoff = 0.05, 
             verbose = TRUE, 
             OrgDb = org.Mridens.eg.db, 
             pAdjustMethod = "none")

l3539list<-l35l39$logFC
names(l3539list)<-l35l39$transcript_id
l3539list = sort(l3539list, decreasing = TRUE)

gsel35mlib <- gseGO(geneList=l3539list, 
                 ont ="ALL", 
                 keyType = "GID", 
                 minGSSize = 3, 
                 maxGSSize = 800, 
                 pvalueCutoff = 0.05, 
                 verbose = TRUE, 
                 OrgDb = org.Mridens.eg.db, 
                 pAdjustMethod = "none")

l35mmixlist<-l35mmix$logFC
names(l35mmixlist)<-l35mmix$transcript_id
l35mmixlist = sort(l35mmixlist, decreasing = TRUE)

gsel35mmix <- gseGO(geneList=l35mmixlist, 
                 ont ="ALL", 
                 keyType = "GID", 
                 minGSSize = 3, 
                 maxGSSize = 800, 
                 pvalueCutoff = 0.05, 
                 verbose = TRUE, 
                 OrgDb = org.Mridens.eg.db, 
                 pAdjustMethod = "none")

mlibmmixlist<-mlibmmix$logFC
names(mlibmmixlist)<-mlibmmix$transcript_id
mlibmmixlist = sort(mlibmmixlist, decreasing = TRUE)

gsemlibmmix <- gseGO(geneList=mlibmmixlist, 
                    ont ="ALL", 
                    keyType = "GID", 
                    minGSSize = 3, 
                    maxGSSize = 800, 
                    pvalueCutoff = 0.05, 
                    verbose = TRUE, 
                    OrgDb = org.Mridens.eg.db, 
                    pAdjustMethod = "none")

dotplot(gse3539, showCategory=10, split=".sign") + facet_grid(.~.sign)
dotplot(gsel35mlib, showCategory=10, split=".sign") + facet_grid(.~.sign)
dotplot(gsel35mmix , showCategory=10, split=".sign") + facet_grid(.~.sign)
dotplot(gsemlibmmix, showCategory=10, split=".sign") + facet_grid(.~.sign)

# 4. clusterprofiler GO enrichment ####

ego <- enrichGO(gene          = l35l39$transcript_id,
                #universe      = gene2go$GID,
                OrgDb         = org.Mridens.eg.db,
                ont           = "all",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.05,
                qvalueCutoff  = 0.05,
                readable      = F,
                keyType = "GID")

head(ego)
dotplot(ego, split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale="free")

