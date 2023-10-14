# The script is used to generate example input and intermediate objects of the GATOM pipeline

library(data.table)
library(KEGGREST)
library(devtools)
library(mwcsr)
load_all()



# Loading full org.Mm.eg.gatom.anno organism annotation object for mouse, as created with following script:
# https://github.com/ctlab/gatom/blob/master/inst/scripts/make_annotation.R
org.Mm.eg.gatom.anno <- readRDS(url("http://artyomovlab.wustl.edu/publications/supp_materials/GATOM/org.Mm.eg.gatom.anno.rds"))

# Loading full KEGG-based network object and corresponding KEGG-based metabolite database object,
# as created with the pipeline https://github.com/ctlab/KEGG-network-pipeline
# (Emelllianova et al, 2022, https://doi.org/10.1093/nar/gkac427)
network <- readRDS(url("http://artyomovlab.wustl.edu/publications/supp_materials/GATOM/network.kegg.rds"))
met.kegg.db <- readRDS(url("http://artyomovlab.wustl.edu/publications/supp_materials/GATOM/met.kegg.db.rds"))

# Example differential gene expression and metabolite abundance data from
# M0 vs M1 macrophage activation comparison
# (Jha et al, 2015, http://doi.org/10.1016/j.immuni.2015.02.005)
met.de.raw <- fread("http://artyomovlab.wustl.edu/publications/supp_materials/GAM/Ctrl.vs.MandLPSandIFNg.met.de.tsv.gz")
gene.de.raw <- fread("http://artyomovlab.wustl.edu/publications/supp_materials/GAM/Ctrl.vs.MandLPSandIFNg.gene.de.tsv.gz")


# Selecting reactions from carbon metabolism pathway, plus Itaconate production reaction
reactionsEx <- c(gsub("^rn:", "", unname(keggLink("reaction", "rn01200"))), "R02243")


# Subsetting the full network to contain only the selected reactions
networkEx <- list()
networkEx$reactions <- network$reactions[reaction %in% reactionsEx]
networkEx$enzyme2reaction <- network$enzyme2reaction[reaction %in% reactionsEx]
networkEx$reaction2align <- network$reaction2align[reaction %in% reactionsEx]
networkEx$atoms <- network$atoms[atom %in% networkEx$reaction2align[, c(atom.x, atom.y)]]
networkEx$metabolite2atom <- network$metabolite2atom[atom %in% networkEx$atoms$atom]

# Subsetting the annotation object to contain only the relevant enzymes
org.Mm.eg.gatom.annoEx <- list()
org.Mm.eg.gatom.annoEx$baseId <- org.Mm.eg.gatom.anno$baseId
org.Mm.eg.gatom.annoEx$gene2enzyme <- org.Mm.eg.gatom.anno$gene2enzyme[enzyme %in% networkEx$enzyme2reaction$enzyme]
org.Mm.eg.gatom.annoEx$genes <- org.Mm.eg.gatom.anno$genes[gene %in% org.Mm.eg.gatom.annoEx$gene2enzyme$gene]
org.Mm.eg.gatom.annoEx$mapFrom <- lapply(org.Mm.eg.gatom.anno$mapFrom,
                                          function(t) t[gene %in% org.Mm.eg.gatom.annoEx$genes$gene])

# Subsetting the metabolite database object to contain only the relevant metabolites
met.kegg.dbEx <- list()
met.kegg.dbEx$baseId <- met.kegg.db$baseId
met.kegg.dbEx$metabolites <- met.kegg.db$metabolites[metabolite %in% networkEx$metabolite2atom$metabolite]
met.kegg.dbEx$anomers <- met.kegg.db$anomers
met.kegg.dbEx$mapFrom <- lapply(met.kegg.db$mapFrom,
                                function(t) t[metabolite %in% met.kegg.dbEx$metabolites$metabolite])

# Subsetting input differential gene expression and metabolite abundance data
genesEx <- union(intersect(org.Mm.eg.gatom.annoEx$mapFrom$RefSeq$RefSeq,
                           gene.de.raw$ID),
                 sample(gene.de.raw$ID, 500))
gene.de.rawEx <- gene.de.raw[gene.de.raw$ID %in% genesEx, ]
met.de.rawEx <- met.de.raw[met.de.raw$ID %in% met.kegg.dbEx$mapFrom$HMDB$HMDB, ]
use_data(met.de.rawEx, gene.de.rawEx)

# Running GATOM pipeline on the example data
gEx <- makeMetabolicGraph(network=networkEx,
                          topology = "atoms",
                          org.gatom.anno=org.Mm.eg.gatom.annoEx,
                          gene.de=gene.de.rawEx,
                          met.db=met.kegg.dbEx,
                          met.de=met.de.rawEx)

gsEx <- scoreGraph(gEx, k.gene=25, k.met=25)
solver <- rnc_solver()
set.seed(42)
res <- solve_mwcsp(solver, gsEx)
mEx <- res$graph

use_data(networkEx, org.Mm.eg.gatom.annoEx, met.kegg.dbEx)
use_data(gEx, gsEx, mEx)
