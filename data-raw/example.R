library(readr)
library(KEGGREST)
library(devtools)

load(url("http://artyomovlab.wustl.edu/publications/supp_materials/GATOM/network.rda"))
load(url("http://artyomovlab.wustl.edu/publications/supp_materials/GATOM/org.Mm.eg.gatom.anno.rda"))
load(url("http://artyomovlab.wustl.edu/publications/supp_materials/GATOM/met.kegg.db.rda"))

met.de.raw <- read_tsv("http://artyomovlab.wustl.edu/publications/supp_materials/GAM/Ctrl.vs.MandLPSandIFNg.met.de.tsv.gz")
gene.de.raw <- read_tsv("http://artyomovlab.wustl.edu/publications/supp_materials/GAM/Ctrl.vs.MandLPSandIFNg.gene.de.tsv.gz")

reactionsEx <- c(gsub("^rn:", "", unname(keggLink("reaction", "rn01200"))),
                 "R02243")

networkEx <- list()
networkEx$reactions <- network$reactions[reaction %in% reactionsEx]
networkEx$enzyme2reaction <- network$enzyme2reaction[reaction %in% reactionsEx]
networkEx$reaction2rpair <- network$reaction2rpair[reaction %in% reactionsEx]
networkEx$rpairs <- network$rpairs[rpair %in% networkEx$reaction2rpair$rpair]
networkEx$rpair2align <- network$rpair2align[rpair %in% networkEx$reaction2rpair$rpair]
networkEx$atoms <- network$atoms[atom %in% networkEx$rpair2align[, c(atom.x, atom.y)]]
networkEx$metabolite2atom <- network$metabolite2atom[atom %in% networkEx$atoms$atom]

org.Mm.eg.gatom.annoEx <- list()
org.Mm.eg.gatom.annoEx$baseId <- org.Mm.eg.gatom.anno$baseId
org.Mm.eg.gatom.annoEx$gene2enzyme <- org.Mm.eg.gatom.anno$gene2enzyme[enzyme %in% networkEx$enzyme2reaction$enzyme]
org.Mm.eg.gatom.annoEx$genes <- org.Mm.eg.gatom.anno$genes[gene %in% org.Mm.eg.gatom.annoEx$gene2enzyme$gene]
org.Mm.eg.gatom.annoEx$mapFrom <- lapply(org.Mm.eg.gatom.anno$mapFrom,
                                          function(t) t[gene %in% org.Mm.eg.gatom.annoEx$genes$gene])


met.kegg.dbEx <- list()
met.kegg.dbEx$baseId <- met.kegg.db$baseId
met.kegg.dbEx$metabolites <- met.kegg.db$metabolites[metabolite %in% networkEx$metabolite2atom$metabolite]
met.kegg.dbEx$anomers <- met.kegg.db$anomers
met.kegg.dbEx$mapFrom <- lapply(met.kegg.db$mapFrom,
                                function(t) t[metabolite %in% met.kegg.dbEx$metabolites$metabolite])


use_data(networkEx, org.Mm.eg.gatom.annoEx, met.kegg.dbEx)

genesEx <- union(intersect(org.Mm.eg.gatom.annoEx$mapFrom$RefSeq$RefSeq,
                           gene.de.raw$ID),
                 sample(gene.de.raw$ID, 500))
gene.de.rawEx <- gene.de.raw[gene.de.raw$ID %in% genesEx, ]
met.de.rawEx <- met.de.raw[met.de.raw$ID %in% met.kegg.dbEx$mapFrom$HMDB$HMDB, ]
use_data(met.de.rawEx, gene.de.rawEx)
