library(gatom)
library(org.Mm.eg.db)
org.Mm.eg.gatom.anno <- makeOrgGatomAnnotation(org.Mm.eg.db)
saveRDS(org.Mm.eg.gatom.anno, file="org.Mm.eg.gatom.anno.rds")

library(org.Hs.eg.db)
org.Hs.eg.gatom.anno <- makeOrgGatomAnnotation(org.Hs.eg.db)
saveRDS(org.Hs.eg.gatom.anno, file="org.Hs.eg.gatom.anno.rds")

library(org.At.tair.db)
org.At.tair.gatom.anno <- makeOrgGatomAnnotation(org.At.tair.db, idColumns = c("Tair"="TAIR",
                                                                               "Entrez"="ENTREZID",
                                                                               "Symbol"="SYMBOL",
                                                                               "Entrez"="ENTREZID"))
saveRDS(org.At.tair.gatom.anno, file="org.At.tair.gatom.anno.rds")

library(org.Sc.sgd.db)
org.Sc.sgd.gatom.anno <- makeOrgGatomAnnotation(org.Sc.sgd.db, idColumns = c("Orf"="ORF",
                                                                             "Entrez"="ENTREZID",
                                                                             "Symbol"="GENENAME",
                                                                             "Ensembl"="ENSEMBL",
                                                                             "Entrez"="ENTREZID"),
                                                nameColumn="GENENAME")
saveRDS(org.Sc.sgd.gatom.anno, file="org.Sc.sgd.gatom.anno.rds")
