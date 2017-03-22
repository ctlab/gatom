library(gatom)
library(org.Mm.eg.db)
org.Mm.eg.gatom.anno <- makeOrgGatomAnnotation(org.Mm.eg.db)
save(org.Mm.eg.gatom.anno, file="org.Mm.eg.gatom.anno.rda")

library(org.Hs.eg.db)
org.Hs.eg.gatom.anno <- makeOrgGatomAnnotation(org.Hs.eg.db)
save(org.Hs.eg.gatom.anno, file="org.Hs.eg.gatom.anno.rda")
