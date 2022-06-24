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
## L35 vs L39 ####
l35l39list<-l35l39$logFC
names(l35l39list)<-l35l39$transcript_id
l35l39list = sort(l35l39list, decreasing = TRUE)

gse3539 <- gseGO(geneList=l35l39list, 
             ont ="ALL", 
             keyType = "GID", 
             minGSSize = 3, 
             maxGSSize = 800, 
             pvalueCutoff = 0.05, 
             verbose = TRUE, 
             OrgDb = org.Mridens.eg.db, 
             pAdjustMethod = "none")

simMatrix <- calculateSimMatrix(gse3539$ID,
                                orgdb="org.Mridens.eg.db",
                                ont="BP",
                                method="Rel",
                                keytype = "GID")
scores <- setNames(-log10(gse3539$qvalues), gse3539$ID)
reducedTerms <- reduceSimMatrix(simMatrix,
                                scores,
                                threshold=0.7,
                                orgdb="org.Mridens.eg.db",
                                keytype = "GID")
heatmapPlot(simMatrix,
            reducedTerms,
            annotateParent=TRUE,
            annotationLabel="parentTerm",
            fontsize=10
)
treemapPlot(reducedTerms)

## Up in L39 (vs L35) ####
upl35l39<-l35l39[l35l39$logFC <= -2,]
l39up<-upl35l39$logFC
names(l39up)<-upl35l39$transcript_id
l39up = sort(l39up, decreasing = TRUE)

gse39 <- gseGO(geneList=l39up, 
                 ont ="ALL", 
                 keyType = "GID", 
                 minGSSize = 3, 
                 maxGSSize = 800, 
                 pvalueCutoff = 0.05, 
                 verbose = TRUE, 
                 OrgDb = org.Mridens.eg.db, 
                 pAdjustMethod = "none")

simMatrix <- calculateSimMatrix(gse39$ID,
                                orgdb="org.Mridens.eg.db",
                                ont="BP",
                                method="Rel",
                                keytype = "GID")
scores <- setNames(-log10(gse39$qvalues), gse39$ID)
reducedTerms <- reduceSimMatrix(simMatrix,
                                scores,
                                threshold=0.7,
                                orgdb="org.Mridens.eg.db",
                                keytype = "GID")
heatmapPlot(simMatrix,
            reducedTerms,
            annotateParent=TRUE,
            annotationLabel="parentTerm",
            fontsize=10
)
treemapPlot(reducedTerms)

## Up in L35 (vs L39) ####
upl35l39<-l35l39[l35l39$logFC >= 2,]
l35up<-upl35l39$logFC
names(l35up)<-upl35l39$transcript_id
l35up = sort(l35up, decreasing = TRUE)

gse35 <- gseGO(geneList=l35up, 
               ont ="ALL", 
               keyType = "GID", 
               minGSSize = 3, 
               maxGSSize = 800, 
               pvalueCutoff = 0.05, 
               verbose = TRUE, 
               OrgDb = org.Mridens.eg.db, 
               pAdjustMethod = "none")

simMatrix <- calculateSimMatrix(gse35$ID,
                                orgdb="org.Mridens.eg.db",
                                ont="BP",
                                method="Rel",
                                keytype = "GID")
scores <- setNames(-log10(gse35$qvalues), gse35$ID)
reducedTerms <- reduceSimMatrix(simMatrix,
                                scores,
                                threshold=0.7,
                                orgdb="org.Mridens.eg.db",
                                keytype = "GID")
heatmapPlot(simMatrix,
            reducedTerms,
            annotateParent=TRUE,
            annotationLabel="parentTerm",
            fontsize=10
)
treemapPlot(reducedTerms)

## L35 vs MLIB ####
l35mliblist<-l35mlib$logFC
names(l35mliblist)<-l35mlib$transcript_id
l35mliblist = sort(l35mliblist, decreasing = TRUE)

gse35mlib <- gseGO(geneList=l35mliblist, 
                 ont ="ALL", 
                 keyType = "GID", 
                 minGSSize = 3, 
                 maxGSSize = 800, 
                 pvalueCutoff = 0.05, 
                 verbose = TRUE, 
                 OrgDb = org.Mridens.eg.db, 
                 pAdjustMethod = "none")

simMatrix <- calculateSimMatrix(gse35mlib$ID,
                                orgdb="org.Mridens.eg.db",
                                ont="BP",
                                method="Rel",
                                keytype = "GID")
scores <- setNames(-log10(gse35mlib$qvalues), gse35mlib$ID)
reducedTerms <- reduceSimMatrix(simMatrix,
                                scores,
                                threshold=0.7,
                                orgdb="org.Mridens.eg.db",
                                keytype = "GID")
heatmapPlot(simMatrix,
            reducedTerms,
            annotateParent=TRUE,
            annotationLabel="parentTerm",
            fontsize=6
)
treemapPlot(reducedTerms)



## Up in L35 (vs MLIB) ####
upl35ml<-l35mlib[l35mlib$logFC >= 2,]
l35up<-upl35ml$logFC
names(l35up)<-upl35ml$transcript_id
l35up = sort(l35up, decreasing = TRUE)

gse35 <- gseGO(geneList=l35up, 
               ont ="ALL", 
               keyType = "GID", 
               minGSSize = 3, 
               maxGSSize = 800, 
               pvalueCutoff = 0.05, 
               verbose = TRUE, 
               OrgDb = org.Mridens.eg.db, 
               pAdjustMethod = "none")

simMatrix <- calculateSimMatrix(gse35$ID,
                                orgdb="org.Mridens.eg.db",
                                ont="BP",
                                method="Rel",
                                keytype = "GID")
scores <- setNames(-log10(gse35$qvalues), gse35$ID)
reducedTerms <- reduceSimMatrix(simMatrix,
                                scores,
                                threshold=0.7,
                                orgdb="org.Mridens.eg.db",
                                keytype = "GID")
heatmapPlot(simMatrix,
            reducedTerms,
            annotateParent=TRUE,
            annotationLabel="parentTerm",
            fontsize=7
)
treemapPlot(reducedTerms)

## Up in MLIB (vs L35) ####
upl35ml<-l35mlib[l35mlib$logFC <= -2,]
mlup<-upl35ml$logFC
names(mlup)<-upl35ml$transcript_id
mlup = sort(mlup, decreasing = TRUE)

gseml <- gseGO(geneList=mlup, 
               ont ="ALL", 
               keyType = "GID", 
               minGSSize = 3, 
               maxGSSize = 800, 
               pvalueCutoff = 0.05, 
               verbose = TRUE, 
               OrgDb = org.Mridens.eg.db, 
               pAdjustMethod = "none")

simMatrix <- calculateSimMatrix(gseml$ID,
                                orgdb="org.Mridens.eg.db",
                                ont="BP",
                                method="Rel",
                                keytype = "GID")
scores <- setNames(-log10(gseml$qvalues), gseml$ID)
reducedTerms <- reduceSimMatrix(simMatrix,
                                scores,
                                threshold=0.7,
                                orgdb="org.Mridens.eg.db",
                                keytype = "GID")
heatmapPlot(simMatrix,
            reducedTerms,
            annotateParent=TRUE,
            annotationLabel="parentTerm",
            fontsize=7
)
treemapPlot(reducedTerms)


## L35 vs MMIX ####
l35mmlist<-l35mmix$logFC
names(l35mmlist)<-l35mmix$transcript_id
l35mmlist = sort(l35mmlist, decreasing = TRUE)

gse35mm <- gseGO(geneList=l35mmlist, 
                   ont ="ALL", 
                   keyType = "GID", 
                   minGSSize = 3, 
                   maxGSSize = 800, 
                   pvalueCutoff = 0.05, 
                   verbose = TRUE, 
                   OrgDb = org.Mridens.eg.db, 
                   pAdjustMethod = "none")

simMatrix <- calculateSimMatrix(gse35mm$ID,
                                orgdb="org.Mridens.eg.db",
                                ont="BP",
                                method="Rel",
                                keytype = "GID")
scores <- setNames(-log10(gse35mm$qvalues), gse35mm$ID)
reducedTerms <- reduceSimMatrix(simMatrix,
                                scores,
                                threshold=0.7,
                                orgdb="org.Mridens.eg.db",
                                keytype = "GID")
heatmapPlot(simMatrix,
            reducedTerms,
            annotateParent=TRUE,
            annotationLabel="parentTerm",
            fontsize=6
)
treemapPlot(reducedTerms)



## Up in L35 (vs MLIB) ####
upl35mm<-l35mmix[l35mmix$logFC >= 2,]
l35up<-upl35mm$logFC
names(l35up)<-upl35mm$transcript_id
l35up = sort(l35up, decreasing = TRUE)

gse35 <- gseGO(geneList=l35up, 
               ont ="ALL", 
               keyType = "GID", 
               minGSSize = 3, 
               maxGSSize = 800, 
               pvalueCutoff = 0.05, 
               verbose = TRUE, 
               OrgDb = org.Mridens.eg.db, 
               pAdjustMethod = "none")

simMatrix <- calculateSimMatrix(gse35$ID,
                                orgdb="org.Mridens.eg.db",
                                ont="BP",
                                method="Rel",
                                keytype = "GID")
scores <- setNames(-log10(gse35$qvalues), gse35$ID)
reducedTerms <- reduceSimMatrix(simMatrix,
                                scores,
                                threshold=0.7,
                                orgdb="org.Mridens.eg.db",
                                keytype = "GID")
heatmapPlot(simMatrix,
            reducedTerms,
            annotateParent=TRUE,
            annotationLabel="parentTerm",
            fontsize=7
)
treemapPlot(reducedTerms)


## Up in MMIX (vs L35) ####
upl35mm<-l35mmix[l35mmix$logFC <= -2,]
upmm<-upl35mm$logFC
names(upmm)<-upl35mm$transcript_id
upmm = sort(upmm, decreasing = TRUE)

gsemm <- gseGO(geneList=upmm, 
               ont ="ALL", 
               keyType = "GID", 
               minGSSize = 3, 
               maxGSSize = 800, 
               pvalueCutoff = 0.05, 
               verbose = TRUE, 
               OrgDb = org.Mridens.eg.db, 
               pAdjustMethod = "none")

simMatrix <- calculateSimMatrix(gsemm$ID,
                                orgdb="org.Mridens.eg.db",
                                ont="BP",
                                method="Rel",
                                keytype = "GID")
scores <- setNames(-log10(gsemm$qvalues), gsemm$ID)
reducedTerms <- reduceSimMatrix(simMatrix,
                                scores,
                                threshold=0.7,
                                orgdb="org.Mridens.eg.db",
                                keytype = "GID")
heatmapPlot(simMatrix,
            reducedTerms,
            annotateParent=TRUE,
            annotationLabel="parentTerm",
            fontsize=7
)
treemapPlot(reducedTerms)

## MLIB vs MMIX ####
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


simMatrix <- calculateSimMatrix(gsemlibmmix$ID,
                                orgdb="org.Mridens.eg.db",
                                ont="BP",
                                method="Rel",
                                keytype = "GID")
scores <- setNames(-log10(gsemlibmmix$qvalues), gsemlibmmix$ID)
reducedTerms <- reduceSimMatrix(simMatrix,
                                scores,
                                threshold=0.7,
                                orgdb="org.Mridens.eg.db",
                                keytype = "GID")
heatmapPlot(simMatrix,
            reducedTerms,
            annotateParent=TRUE,
            annotationLabel="parentTerm",
            fontsize=6
)
treemapPlot(reducedTerms)




