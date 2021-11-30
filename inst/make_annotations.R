library(gatom)
library(org.Mm.eg.db)
org.Mm.eg.gatom.anno <- makeOrgGatomAnnotation(org.Mm.eg.db)
org.Mm.eg.gatom.anno$pathways <- getPathways2annotation(org.gatom.anno=org.Mm.eg.gatom.anno,
                                                        organism="mmu")
save(org.Mm.eg.gatom.anno, file="org.Mm.eg.gatom.anno.rda")

library(org.Hs.eg.db)
org.Hs.eg.gatom.anno <- makeOrgGatomAnnotation(org.Hs.eg.db)
org.Hs.eg.gatom.anno$pathways <- getPathways2annotation(org.gatom.anno=org.Hs.eg.gatom.anno,
                                                        organism="hsa")
save(org.Hs.eg.gatom.anno, file="org.Hs.eg.gatom.anno.rda")
