
#' gatom: a package for finding an active metabolic module in atom transition network
#'
#' This package implements a metabolic network analysis pipeline to
#' identify an active metabolic module based on high throughput data.
#' The pipeline takes as input transcriptional and/or metabolic data
#' and finds a metabolic subnetwork (module) most regulated between the two
#' conditions of interest. The package further provides functions for module
#' post-processing, annotation and visualization.
#'
#' @section Functions:
#'
#' Data preprocessing:
#' \code{\link{prepareDE}}, \code{\link{getMetDEMeta}}, \code{\link{getGeneDEMeta}}
#'
#' Graph creation:
#' \code{\link{makeMetabolicGraph}}
#'
#' Graph scoring:
#' \code{\link{scoreGraph}}
#'
#' Module postprocessing:
#' \code{\link{collapseAtomsIntoMetabolites}}, \code{\link{connectAtomsInsideMetabolite}},
#' \code{\link{addHighlyExpressedEdges}}, \code{\link{abbreviateLabels}}
#'
#' Plotting module:
#' \code{\link{createShinyCyJSWidget}}
#'
#' Exporting module:
#' \code{\link{saveModuleToHtml}}, \code{\link{saveModuleToDot}},
#' \code{\link{saveModuleToPdf}}, \code{\link{saveModuleToXgmml}}
#'
#' For detailed pipeline analysis, see gatom vignette.
#'
#' @section Example Data:
#' Example data provided by gatom consists of:
#' metabolite differential abundance data (\code{\link{met.de.rawEx}}),
#' gene differential expression data (\code{\link{gene.de.rawEx}}),
#' KEGG-based network object (\code{\link{networkEx}}),
#' KEGG-based metabolite database object (\code{\link{met.kegg.dbEx}}),
#' Example organism annotation object (\code{\link{org.Mm.eg.gatom.annoEx}}),
#' metabolic graph with atom topology (\code{\link{gEx}}),
#' scored metabolic graph with atom topology (\code{\link{gsEx}}),
#' and metabolic module (\code{\link{mEx}}).
#'
#' @docType package
#' @name gatom
NULL

#' Example metabolic graph with atom topology.
#'
#' See file \url{https://github.com/ctlab/gatom/blob/master/inst/scripts/example.R}
#' for details.
#'
#' @docType data
#' @name gEx
#' @format igraph object
NULL

#' Example scored metabolic graph with atom topology.
#'
#' See file \url{https://github.com/ctlab/gatom/blob/master/inst/scripts/example.R}
#' for details.
#'
#' @docType data
#' @name gsEx
#' @format igraph object
NULL

#' Example metabolic module.
#'
#' See file \url{https://github.com/ctlab/gatom/blob/master/inst/scripts/example.R}
#' for details.
#'
#' @docType data
#' @name mEx
#' @format igraph object
NULL

#' Example gene differential expression data.
#'
#' See file \url{https://github.com/ctlab/gatom/blob/master/inst/scripts/example.R}
#' for details.
#'
#' @docType data
#' @name gene.de.rawEx
#' @format tibble/data.frame object
NULL

#' Example metabolite differential abundance data.
#'
#' See file \url{https://github.com/ctlab/gatom/blob/master/inst/scripts/example.R}
#' for details.
#'
#' @docType data
#' @name met.de.rawEx
#' @format tibble/data.frame object
NULL

#' Example KEGG-based network object
#'
#' See file \url{https://github.com/ctlab/gatom/blob/master/inst/scripts/example.R}
#' for details.
#'
#' @docType data
#' @name networkEx
#' @format list object
NULL

#' Example KEGG-based metabolite database object
#'
#' See file \url{https://github.com/ctlab/gatom/blob/master/inst/scripts/example.R}
#' for details.
#'
#' @docType data
#' @name met.kegg.dbEx
#' @format list object
NULL

#' Example organism annotation object
#'
#' See file \url{https://github.com/ctlab/gatom/blob/master/inst/scripts/example.R}
#' for details.
#'
#' @docType data
#' @name org.Mm.eg.gatom.annoEx
#' @format list object
NULL

