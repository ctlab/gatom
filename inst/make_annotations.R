library(gatom)
library(org.Mm.eg.db)
org.Mm.eg.gatom.anno <- makeOrgGatomAnnotation(org.Mm.eg.db)
saveRDS(org.Mm.eg.gatom.anno, file="org.Mm.eg.gatom.anno.rds")

library(org.Hs.eg.db)
org.Hs.eg.gatom.anno <- makeOrgGatomAnnotation(org.Hs.eg.db)
saveRDS(org.Hs.eg.gatom.anno, file="org.Hs.eg.gatom.anno.rds")
