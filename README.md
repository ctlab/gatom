[![Travis-CI Build Status](https://travis-ci.org/ctlab/gatom.svg?branch=master)](https://travis-ci.org/ctlab/gatom)
[![codecov](https://codecov.io/gh/ctlab/gatom/branch/master/graph/badge.svg)](https://codecov.io/gh/ctlab/gatom)


# gatom

An R-package for finding active metabolic modules in atom transition network.

---

Full vignette can be found here [https://github.com/ctlab/gatom/blob/master/inst/Using_gatom_package.html]: 

### Installation 

```{r}
library(devtools)
install_github("ctlab/gatom")
```

### Quick start

```{r message=FALSE}
library(gatom)
library(data.table)
library(igraph)
```

First let's load data with atom mappings (`network` object),
enzyme annotations for mouse (`org.Mm.eg.gatom`)
and metabolite annotations (`met.kegg.db.rda`):

```{r}
data("networkEx")
data("org.Mm.eg.gatom.annoEx")
data("met.kegg.dbEx")
```

Loading input data:

```{r message=F}
data("met.de.rawEx")
data("gene.de.rawEx")

```

Getting atom graph:

```{r}
g <- makeAtomGraph(network=networkEx, 
                   org.gatom.anno=org.Mm.eg.gatom.annoEx, 
                   gene.de=gene.de.rawEx,
                   met.db=met.kegg.dbEx, 
                   met.de=met.de.rawEx)
print(g)
```

Scoring graph:

```{r}
gs <- scoreGraph(g, k.gene = 25, k.met=25)
```

Finding a module:

```{r}
set.seed(42)
m <- solveSgmwcsRandHeur(gs, max.iterations = 2000)
```

```{r}
print(m)
E(m)$label
V(m)$label
```

We can save the module to different formats (dot, xgmml, svg, pdf):

```{r echo=-2, message=F, warning=F, out.height=550}
saveModuleToPdf(m, file="M0.vs.M1.pdf", name="M0.vs.M1", n_iter=100, force=1e-5, seed=1)
knitr::include_graphics("M0.vs.M1.pdf")
```
