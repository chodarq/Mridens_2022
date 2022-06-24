library(rrvgo)
#L35 vs L39
go_analysis <- read.delim(system.file("extdata/example.txt", package="rrvgo"))
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
            fontsize=6
            )
scatterPlot(simMatrix, reducedTerms,size = "size")
treemapPlot(reducedTerms)

#L35 vs LMLIB
go_analysis <- read.delim(system.file("extdata/example.txt", package="rrvgo"))
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
            fontsize=6)
scatterPlot(simMatrix, reducedTerms)
treemapPlot(reducedTerms)
