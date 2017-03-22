library(KEGGREST)
library(data.table)


reaction2enzyme <- keggLink("enzyme", "reaction")
reaction2enzyme <- data.table(
    reaction=gsub("rn:", "", names(reaction2enzyme), fixed = T),
    enzyme=gsub("ec:", "", reaction2enzyme, fixed = T))

reactions <- keggList("reaction")
reactions <- data.table(
        reaction=gsub("rn:", "", names(reactions), fixed = T),
        reaction_name=gsub(";.*$", "", reactions),
        reaction_equation=gsub("^.*; ", "", reactions))

reactions[, reaction_url := sprintf("http://www.kegg.jp/entry/%s", reaction)]

enzyme2reaction <- copy(reaction2enzyme)
setkey(enzyme2reaction, enzyme)


network <- list()
network$reactions <- reactions
network$enzyme2reaction <- enzyme2reaction

network$rpairs <- fread("./network/rpairs.tsv")
network$rpair2align <- fread("./network/rpair2align.tsv")
setkey(network$rpair2align, rpair)
network$reaction2rpair <- fread("./network/rpair2reaction.tsv")
setkey(network$reaction2rpair, reaction)

network$atoms <- data.table(atom=union(network$rpair2align$atom.x,
                                       network$rpair2align$atom.y))

network$atoms[, metabolite := gsub("_.*", "", atom)]
network$atoms[, element := "C"]
setkey(network$atoms, atom)

network$metabolite2atom <- network$atoms[, list(metabolite, atom)]
setkey(network$metabolite2atom, metabolite)


save(network, file="network.rda")

