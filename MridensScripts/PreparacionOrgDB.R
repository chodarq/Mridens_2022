library(dplyr)
library(stringr)
library(AnnotationForge)
library(jsonlite)
library(purrr)
library(RCurl)
library(tibble)

# 1. Comparación resultados trinotate vs eggnog ####
# resultados eggnog
idsegg<-read.csv("eggnog.csv",header=F)
colnames(idsegg)<-"txid"
# resultados tri
idtri<-read.csv("trinotate.csv",header=F)
colnames(idtri)<-"txid"

# comunes ambas listas
comunes<-intersect(idsegg,idtri)
exEgg<-setdiff(idsegg,idtri)
exTri<-setdiff(idtri,idsegg)
# Decidimos seguir con eggnog

# 2. Carga anotacion eggnog ####
eggann<-read.csv("Mridens.emapper.annotations",sep="\t",header=T,skip =4)
eggann<-head(eggann,-10)
names(eggann)[1]<-"GID"

## 2.1 Basar anotación en isoforma más larga 
lista<-read.csv("lista.isoformas.csv",sep="\t",header=T)
names(lista)[1]<-"GID"
anno<-inner_join(lista,eggann,by="GID")

# 3. Preparación de frames para archivo anotación ####
## 3.1 Informacion de los genes ####
mridgen<-anno[c(2,9)]
mridgen <- mridgen[mridgen[,2]!="-",]
colnames(mridgen)<-c("GID","DESCRIPTION")

## 3.2 Selección columnas KEGG_ko ####
update_kegg <- function(json = "ko00001.json",file=NULL) {
  pathway2name <- tibble(Pathway = character(), Name = character())
  ko2pathway <- tibble(Ko = character(), Pathway = character())
  kegg <- fromJSON(json)
  
  for (a in seq_along(kegg[["children"]][["children"]])) {
    A <- kegg[["children"]][["name"]][[a]]
    
    for (b in seq_along(kegg[["children"]][["children"]][[a]][["children"]])) {
      B <- kegg[["children"]][["children"]][[a]][["name"]][[b]] 
      
      for (c in seq_along(kegg[["children"]][["children"]][[a]][["children"]][[b]][["children"]])) {
        pathway_info <- kegg[["children"]][["children"]][[a]][["children"]][[b]][["name"]][[c]]
        
        pathway_id <- str_match(pathway_info, "ko[0-9]{5}")[1]
        pathway_name <- str_replace(pathway_info, " \\[PATH:ko[0-9]{5}\\]", "") %>% str_replace("[0-9]{5} ", "")
        pathway2name <- rbind(pathway2name, tibble(Pathway = pathway_id, Name = pathway_name))
        
        kos_info <- kegg[["children"]][["children"]][[a]][["children"]][[b]][["children"]][[c]][["name"]]
        
        kos <- str_match(kos_info, "K[0-9]*")[,1]
        
        ko2pathway <- rbind(ko2pathway, tibble(Ko = kos, Pathway = rep(pathway_id, length(kos))))
      }
    }
  }
  
  save(pathway2name, ko2pathway, file = file)
}

update_kegg(json = "ko00001.json",file="kegg_info.RData")
load("/media/chodar/DATOS/Github.Projects/PIrles_2022/kegg_info.RData")

gene2ko <- anno %>%
  dplyr::select(GID = Gen, KO = KEGG_ko) %>%
  na.omit()

colnames(ko2pathway)=c("KO","Pathway")
gene2ko$KO=str_replace(gene2ko$KO,"ko:","")

gene2pathway <- gene2ko %>% left_join(ko2pathway, by = "KO") %>% 
  dplyr::select(GID, Pathway) %>%
  na.omit()

## 3.3 Selección columna GOs ####
mridGO<-anno[c(2,11)]
names(mridGO)<-c("GID","GO","evidence")
mridGO <- mridGO[mridGO[,2]!="-",]
mridGO$evidence<-rep("IEA",dim(mridGO)[1])
all_go_list=str_split(mridGO$GO,",")
gene2go <- data.frame(GID = rep(mridGO$GID,
                                times = sapply(all_go_list, length)),
                      GO = unlist(all_go_list),
                      EVIDENCE = "IEA")

## 3.4 Creación de paquete Org.db ####
makeOrgPackage(gene_info=mridgen,
               go=gene2go,
               ko=gene2ko,
               maintainer='chodar <chodar@inta.uchile.cl>',
               author='chodar <chodar@inta.uchile.cl>',
               pathway=gene2pathway,
               version="0.0.1",
               outputDir = ".",
               tax_id=1567044,
               genus="Mastrus",
               species="ridens",
               goTable="go")