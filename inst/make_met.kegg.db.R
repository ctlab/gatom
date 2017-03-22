library(KEGGREST)
library(data.table)

metabolites <- keggList("compound")
metabolites <- data.table(
        metabolite=gsub("cpd:", "", names(metabolites), fixed = T),
        metabolite_name=gsub(";.*$", "", metabolites))
metabolites[, metabolite_url := sprintf("http://www.kegg.jp/entry/%s", metabolite)]
setkey(metabolites, metabolite)


HMDB2metabolite <- fread("hmdb2kegg.tsv")
HMDB2metabolite <- HMDB2metabolite[, list(HMDB=HMDB, metabolite=KEGG)]
setkey(HMDB2metabolite, HMDB)

anomers <- fread("mets2collapse.tsv")[, list(metabolite=from, base_metabolite=to)]
anomers <- unique(anomers)

met.kegg.db <- list()
met.kegg.db$metabolites <- metabolites
met.kegg.db$mapFrom <- list("HMDB"=HMDB2metabolite)
met.kegg.db$baseId <- "KEGG"
met.kegg.db$anomers <- list(
    metabolite2base_metabolite=copy(anomers),
    base_metabolite2metabolite=copy(anomers))
setkey(met.kegg.db$anomers$metabolite2base_metabolite, metabolite)
setkey(met.kegg.db$anomers$base_metabolite2metabolite, base_metabolite)

save(met.kegg.db, file="met.kegg.db.rda")
