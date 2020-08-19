
[![Travis-CI Build
Status](https://travis-ci.org/ctlab/gatom.svg?branch=master)](https://travis-ci.org/ctlab/gatom)
[![codecov](https://codecov.io/gh/ctlab/gatom/branch/master/graph/badge.svg)](https://codecov.io/gh/ctlab/gatom)

# gatom

An R-package for finding active metabolic modules in atom transition
network.

Full vignette can be found
[here](https://rawgit.com/ctlab/gatom/master/inst/Using_gatom_package.html).

### Installation

``` r
library(devtools)
install_github("ctlab/mwcsr")
install_github("ctlab/gatom")
```

### Quick start

``` r
library(gatom)
library(data.table)
library(igraph)
library(mwcsr)
```

First let’s load data with atom mappings (`network` object), enzyme
annotations for mouse (`org.Mm.eg.gatom`) and metabolite annotations
(`met.kegg.db.rda`):

``` r
data("networkEx")
data("org.Mm.eg.gatom.annoEx")
data("met.kegg.dbEx")
```

Loading input data:

``` r
data("met.de.rawEx")
data("gene.de.rawEx")
```

Getting atom graph:

``` r
g <- makeAtomGraph(network=networkEx,
                   org.gatom.anno=org.Mm.eg.gatom.annoEx,
                   gene.de=gene.de.rawEx,
                   met.db=met.kegg.dbEx,
                   met.de=met.de.rawEx)
```

    ## Registered S3 method overwritten by 'pryr':
    ##   method      from
    ##   print.bytes Rcpp

    ## Found DE table for genes with RefSeq IDs

    ## Found DE table for metabolites with HMDB IDs

``` r
print(g)
```

    ## IGRAPH 710d1d9 UN-- 194 209 -- 
    ## + attr: name (v/c), metabolite (v/c), element (v/c), label (v/c), url
    ## | (v/c), pval (v/n), origin (v/n), HMDB (v/c), log2FC (v/n), baseMean
    ## | (v/n), logPval (v/n), signal (v/c), signalRank (v/n), label (e/c),
    ## | pval (e/n), origin (e/n), RefSeq (e/c), gene (e/c), enzyme (e/c),
    ## | reaction_name (e/c), reaction_equation (e/c), url (e/c), reaction
    ## | (e/c), rpair (e/c), log2FC (e/n), baseMean (e/n), logPval (e/n),
    ## | signal (e/c), signalRank (e/n)
    ## + edges from 710d1d9 (vertex names):
    ## [1] C00022_2--C00024_0 C00022_0--C00024_1 C00025_0--C00026_0 C00025_1--C00026_1
    ## [5] C00025_2--C00026_2 C00025_4--C00026_4 C00025_7--C00026_7 C00024_1--C00033_0
    ## + ... omitted several edges

Scoring graph, obtaining an instance of SGMWCS (Signal Generalize
Maximum Weight Subgraph) problem instance:

``` r
gs <- scoreGraph(g, k.gene = 25, k.met=25)
```

Initialize an SMGWCS solver (an heuristic Virgo solver is used for
simplicity, check out `mwcsr` package documentation for more options):

``` r
vhsolver <- virgo_solver(cplex_dir=NULL)
```

Finding a module:

``` r
res <- solve_mwcsp(vhsolver, gs)
m <- res$graph
```

``` r
print(m)
```

    ## IGRAPH eb68808 UNW- 103 102 -- 
    ## + attr: signals (g/n), name (v/c), metabolite (v/c), element (v/c),
    ## | label (v/c), url (v/c), pval (v/n), origin (v/n), HMDB (v/c), log2FC
    ## | (v/n), baseMean (v/n), logPval (v/n), signal (v/c), signalRank (v/n),
    ## | weight (v/n), index (v/n), label (e/c), pval (e/n), origin (e/n),
    ## | RefSeq (e/c), gene (e/c), enzyme (e/c), reaction_name (e/c),
    ## | reaction_equation (e/c), url (e/c), reaction (e/c), rpair (e/c),
    ## | log2FC (e/n), baseMean (e/n), logPval (e/n), signal (e/c), signalRank
    ## | (e/n), weight (e/n), index (e/n)
    ## + edges from eb68808 (vertex names):
    ## [1] C00025_0--C00026_0 C00025_2--C00026_2 C00025_4--C00026_4 C00025_7--C00026_7
    ## + ... omitted several edges

``` r
head(E(m)$label)
```

    ## [1] "Psat1" "Psat1" "Psat1" "Psat1" "Gpt2"  "Got2"

``` r
head(V(m)$label)
```

    ## [1] "Pyruvate"       "L-Glutamate"    "L-Glutamate"    "L-Glutamate"   
    ## [5] "L-Glutamate"    "2-Oxoglutarate"

We can save the module to different formats (dot, xgmml, svg,
pdf):

``` r
saveModuleToPdf(m, file="M0.vs.M1.pdf", name="M0.vs.M1", n_iter=100, force=1e-5, seed=1)
```

![Module](https://rawgit.com/ctlab/gatom/master/inst/M0.vs.M1.pdf.png)
