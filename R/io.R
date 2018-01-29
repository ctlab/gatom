#' Save module to a nice pdf file
#' @param module Module to save
#' @param layout Edge layout
#' @param n_iter Number of repel algorithm iterations
#' @param force Value of repel force
#' @param file File to save to
#' @param name Name of the module
#' @export
#' @import GGally
#' @import igraph
saveModuleToPdf <- function(module, layout="sna_kk", n_iter=100, force=1e-5, file, name=NULL) {

    pdflayout <- getModulePdfLayout(module, layout, n_iter, force)
    layout2 <- pdflayout$layout2
    node_attrs <- getPdfModuleAttrs(module)$produce_node_attrs
    edge_attrs <- getPdfModuleAttrs(module)$produce_edge_attrs

    pdf(file=file, width=pdflayout$gwidth, height=pdflayout$gheight)
    plot(ggnet2(module, mode=layout2$layouts[[length(layout2$layouts)]],
                layout.exp=0.3 * (60 / length(V(module))),
                size=node_attrs$width, max_size=25, node.color=node_attrs$color,
                node.label=V(module)$label, label.size=node_attrs$fontsize, label.color="grey13",
                edge.size=edge_attrs$penwidth, edge.color=edge_attrs$color,
                edge.label.fill=NA, edge.label=E(module)$label, edge.label.size=edge_attrs$fontsize,
                legend.size=0, legend.position="up") +
         ggplot2::ggtitle(name) +
         ggplot2::theme(plot.title=
                            ggplot2::element_text(size=max(c(node_attrs$fontsize, edge_attrs$fontsize)) * 5)))
    dev.off()
}

getPdfModuleAttrs <- function(module) {
    produce_node_attrs <- getPdfNodeStyleAttributes(as_data_frame(module, what="vertices"))
    fn <- sapply(produce_node_attrs, is.factor)
    produce_node_attrs[fn] <- lapply(produce_node_attrs[fn], as.character)
    produce_edge_attrs <- getPdfEdgeStyleAttributes(as_data_frame(module))
    fe <- sapply(produce_edge_attrs, is.factor)
    produce_edge_attrs[fe] <- lapply(produce_edge_attrs[fe], as.character)

    produce_node_attrs$fontsize <- produce_node_attrs$fontsize
    produce_edge_attrs$fontsize <- produce_edge_attrs$fontsize

    produce_node_attrs$width <- produce_node_attrs$width
    produce_edge_attrs$penwidth <- produce_edge_attrs$penwidth

    return(list(
        produce_node_attrs=produce_node_attrs,
        produce_edge_attrs=produce_edge_attrs
    ))
}


#' @import grid
getModulePdfLayout <- function(module, layout, n_iter, force) {

    layout1 <- make_layout(module, layout=layout)

    # produce labels' grobs on canvas of a particular size
    # by modifying geom_text_repel algorithm (see function force_alg)
    gwidth <- 2.5 * (max(layout1[,1]) - min(layout1[,1]))
    gheight <- 2.5 * (max(layout1[,2]) - min(layout1[,2]))
    layout1 <- range01(layout1)

    pdf(file=paste0(tempdir(), "/device.pdf"), width=gwidth, height=gheight)

    node_attrs <- getPdfModuleAttrs(module)$produce_node_attrs
    edge_attrs <- getPdfModuleAttrs(module)$produce_edge_attrs

    box.padding <- unit(0.07 * node_attrs$fontsize, "lines")
    box_padding_x <- convertWidth(box.padding, "npc", valueOnly=TRUE)
    box_padding_y <- convertHeight(box.padding, "npc", valueOnly=TRUE)

    edges <- cbind(store_all_info(module, layout1)$lines_to_check$IDs_num,
                   store_all_info(module, layout1)$lines_to_check$IDe_num)

    # produce grobs
    nlabel_semisizes <- lapply(seq(node_attrs$label), function(i) {
        t <- textGrob(node_attrs$label[i],
                      x=store_all_info(module, layout1)$nt$x[i],
                      y=store_all_info(module, layout1)$nt$y[i],
                      default.units="npc",

                      gp=gpar(fontsize=node_attrs$fontsize[i] * ggplot2::.pt
                      ),
                      name="text"
        )

        gw <- convertWidth(grobWidth(t), "npc", TRUE) / 2
        gh <- convertHeight(grobHeight(t), "npc", TRUE) / 2

        c(
            "semi_w" = gw + box_padding_x[i],
            "semi_h" = gh + box_padding_y[i]
        )
    })
    nlabel_semisizes <- data.frame(do.call(rbind, nlabel_semisizes))

    elabel_semisizes <- lapply(seq(edge_attrs$label), function(i) {
        t <- textGrob(edge_attrs$label[i],
                      x=store_all_info(module, layout1)$lines_to_check$middle_x[i],
                      y=store_all_info(module, layout1)$lines_to_check$middle_y[i],
                      default.units="npc",

                      gp=gpar(fontsize=edge_attrs$fontsize[i] * ggplot2::.pt
                      ),
                      name="text"
        )

        gw <- convertWidth(grobWidth(t), "npc", TRUE) / 2
        gh <- convertHeight(grobHeight(t), "npc", TRUE) / 2

        c(
            "semi_w" = gw + box_padding_x[i],
            "semi_h" = gh + box_padding_y[i]
        )
    })
    elabel_semisizes <- data.frame(do.call(rbind, elabel_semisizes))

    # rescale layout1 to put all boxes inside [0:1]
    # new_borders_x:
    new_borders_x_bot <- max(nlabel_semisizes$semi_w)
    new_borders_x_top <- 1 - max(nlabel_semisizes$semi_w)
    rescale_par_x <- new_borders_x_top - new_borders_x_bot

    # new_borders_h:
    new_borders_y_bot <- max(nlabel_semisizes$semi_h)
    new_borders_y_top <- 1 - max(nlabel_semisizes$semi_h)
    rescale_par_y <- new_borders_y_top - new_borders_y_bot

    layout1[,1] <- new_borders_x_bot + layout1[,1] * rescale_par_x
    layout1[,2] <- new_borders_y_bot + layout1[,2] * rescale_par_y

    dev.off()

    layout2 <- force_alg(layout1,
                         nlabel_semisizes,
                         elabel_semisizes,
                         edges,
                         xlim=new_borders_x_bot,
                         ylim=new_borders_y_bot,
                         n_iter, force)

    return(list(
        layout2=layout2,
        gwidth=gwidth,
        gheight=gheight)
    )
}

getPdfSize <- function(logPval) {
    logPval[is.na(logPval)] <- -5
    return(pmin(0.2 - logPval/100/0.3, 0.5))
}

getPdfNodeStyleAttributes <- function(attrs) {
    logPval <- if (!is.null(attrs$logPval)) attrs$logPval else -5
    with(attrs, data.frame(
        label=if (!is.null(attrs$label)) label else "",
        shape=if (!is.null(attrs$nodeType)) nodeShapeMap[nodeType] else "circle",
        fixedsize="true",
        style="filled",
        width=sapply(logPval, getPdfSize) * 45,
        fontsize=sapply(logPval, getPdfSize) * 25,
        color=if (!is.null(attrs$log2FC)) sapply(log2FC, getDotColor) else "slategrey",
        fillcolor=if (!is.null(attrs$log2FC)) sapply(log2FC, getDotColor) else "white"
    ))
}

getPdfEdgeStyleAttributes <- function(attrs) {
    logPval <- if (!is.null(attrs$logPval)) attrs$logPval else -5
    with(attrs, data.frame(
        label=if (!is.null(attrs$label)) label else "",
        style=if (!is.null(attrs$rptype)) edgeStyleMap[rptype] else "solid",
        penwidth=sapply(logPval, getPdfSize) * 25,
        fontsize=sapply(logPval, getPdfSize) * 25,
        color=if (!is.null(attrs$log2FC)) sapply(log2FC, getDotColor) else "grey"
    ))
}



#' @import XML
sanitizeForXml <- function (string) {
    xmlValue(xmlTextNode(string))
}

#' Save module to an XGMML file
#' @param module Module to save
#' @param file File to save to
#' @param name Name of the module
#' @export
saveModuleToXgmml <- function(module, file, name=NULL) {
    if (is.null(name)) {
        name <- deparse(substitute(module))
    }
    s <- getModuleXmlString(module, name)
    write(s, file)
}

getModuleXmlString <- function(module, name) {
    res <- c()
    res <- c(res, '<?xml version="1.0"?>\n')
    res <- c(res, paste0('<graph label="', name, '" xmlns:dc="http://purl.org/dc/elements/1.1/" xmlns:xlink="http://www.w3.org/1999/xlink" xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#" xmlns:cy="http://www.cytoscape.org" xmlns="http://www.cs.rpi.edu/XGMML">\n'))
    res <- c(res, '  <att name="documentVersion" value="1.1"/>\n')
    res <- c(res, '  <att name="networkMetadata">\n')
    res <- c(res, getRDFXmlString(module, name, indent="    "))
    res <- c(res, '  </att>\n')
    res <- c(res, getNodeXmlStrings(module, indent="  "))
    res <- c(res, getEdgeXmlStrings(module, indent="  "))
    res <- c(res, '</graph>\n')

    paste(res, collapse="")
}

xmlNodeString <- function(name, text) {
    paste0("<", name, ">", text, "</", name, ">")
}

getRDFXmlString <- function(module, name, indent="") {
    res <- c()
    res <- c(res, paste0("<rdf:RDF>\n"))
    res <- c(res, paste0("  ", "<rdf:Description rdf:about=\"http://www.cytoscape.org/\">\n"))
    res <- c(res, paste0("    ", xmlNodeString("dc:type", "N/A"), "\n"))
    res <- c(res, paste0("    ", xmlNodeString("dc:description", "N/A"), "\n"))
    res <- c(res, paste0("    ", xmlNodeString("dc:identifier", "N/A"), "\n"))
    res <- c(res, paste0("    ", xmlNodeString("dc:date", Sys.time()), "\n"))
    res <- c(res, paste0("    ", xmlNodeString("dc:title", name), "\n"))
    res <- c(res, paste0("    ", xmlNodeString("dc:format", "Cytoscape-XGMML"), "\n"))
    res <- c(res, paste0("  ", "</rdf:Description>\n"))
    res <- c(res, paste0("</rdf:RDF>\n"))

    res <- paste0(indent, res)
    paste(res, collapse="")
}


getAttrXmlStrings <- function(attr.values, indent="") {
    attr.xmlStrings <- list()
    attr.names <- names(attr.values)
    for(i in seq_along(attr.values)) {
        attr.rtype <- is(attr.values[[i]])[1]
        if (attr.rtype == "character") {
            type <- "string"
        } else if(attr.rtype == "integer") {
            type <- "integer"
        } else if(attr.rtype == "numeric") {
            type <- "real"
            attr.values[[i]] <- sub("Inf", "Infinity", as.character(attr.values[[i]]), fixed=TRUE)
        } else {
            type <- "string"
            attr.values[[i]] <- as.character(attr.values[[i]])
        }

        if(type=="string"){
            attr.values[[i]] <- sapply(as.vector(attr.values[[i]]), sanitizeForXml)
        }

        attr.xmlStrings[[i]] <- paste0(indent, "<att type=", "\"", type, "\"", " name=", "\"", attr.names[i], "\"", " value=", "\"", attr.values[[i]], "\"/>\n")
        attr.xmlStrings[[i]][is.na(attr.values[[i]])] <- NA
    }
    names(attr.xmlStrings) <- attr.names
    do.call("cbind", attr.xmlStrings)
}

getNodeXmlStrings <- function(module, indent="") {
    if (length(V(module)) == 0) {
        return(NULL)
    }
    attr.values <- as_data_frame(module, what="vertices")
    attr.xmlStrings <- getAttrXmlStrings(attr.values, paste0(indent, "  "))

    if(is.null(V(module)$name))
    {
        V(module)$name <- as.character(V(module))
    }
    node.label <- V(module)$name
    node.id <- as.vector(V(module))
    xmlHeaders <- paste0(indent,
                         "<node",
                         " label=", "\"", node.label, "\"",
                         " id=", "\"", node.id, "\"",
                         ">\n")
    xmlFooter <- paste0(indent, "</node>\n")
    xmlStrings <- paste0(xmlHeaders,
                         apply(attr.xmlStrings, 1, function(x) paste(na.omit(x), collapse="")),
                         xmlFooter)

    xmlStrings
}

getEdgeXmlStrings <- function(module, indent="") {
    if (length(E(module)) == 0) {
        return(NULL)
    }
    attr.values <- as_data_frame(module, what="edges")
    attr.xmlStrings <- getAttrXmlStrings(attr.values, paste0(indent, "  "))

    edgelist.names <- get.edgelist(module, names=TRUE)
    edgelist.names <- paste(edgelist.names[,1], edgelist.names[,2], sep=" (pp) ")
    edgelist.ids <- get.edgelist(module, names=FALSE)

    xmlHeaders <- paste0(indent,
                         "<edge",
                         " label=", "\"", edgelist.names, "\"",
                         " source=", "\"", edgelist.ids[,1], "\"",
                         " target=", "\"", edgelist.ids[,2], "\"",
                         ">\n")
    xmlFooter <- paste0(indent, "</edge>\n")
    xmlStrings <- paste0(xmlHeaders,
                         if (is.null(attr.xmlStrings)) {
                             NULL
                         } else {
                             apply(attr.xmlStrings, 1, function(x) paste(na.omit(x), collapse=""))
                         },
                         xmlFooter)

    xmlStrings
}

#' Save module to a graphviz dot file
#' @param module Module to save
#' @param file File to save to
#' @param name Name of the module
#' @export
saveModuleToDot <- function(module, file, name=NULL) {
    if (is.null(name)) {
        name <- deparse(substitute(module))
    }
    s <- getGraphDotString(module, name)
    write(s, file)
}

getGraphDotString <- function(module, name) {
    res <- c()
    res <- c(res, sprintf('graph "%s" {\n', name))
    res <- c(res, "outputorder=edgesfirst\n")
    res <- c(res, "bgcolor=transparent\n")
    res <- c(res, "pad=\"2,0.25\"\n")
    res <- c(res, getEdgeDotStrings(module, indent="  "))
    res <- c(res, getNodeDotStrings(module, indent="  "))
    res <- c(res, '}\n')

    paste(res, collapse="")
}


getAttrDotStrings <- function(attr.values) {
    attr.dotStrings <- list()
    attr.names <- names(attr.values)
    for(i in seq_along(attr.values)) {
        attr.dotStrings[[i]] <- sprintf("%s = \"%s\"", attr.names[i], attr.values[[i]])
        attr.dotStrings[[i]][is.na(attr.values[[i]])] <- NA
    }
    names(attr.dotStrings) <- attr.names
    do.call("cbind", attr.dotStrings)
}

upRamp <- colorRamp(c("#cccccc", "#ff3333"))
downRamp <- colorRamp(c("#cccccc", "green"))

getDotColor <- function(log2FC) {
    if (is.na(log2FC)) {
        return("#7777ff")
    }

    if (log2FC > 0) {
        col <- upRamp(min(1, abs(log2FC) / 2))
    } else {
        col <- downRamp(min(1, abs(log2FC) / 2))
    }
    rgb(col[,1], col[,2], col[,3], maxColorValue=255)
}

getDotSize <- function(logPval) {
    logPval[is.na(logPval)] <- 0
    return(pmin(0.2 - logPval/100/0.3, 0.5))
}


nodeShapeMap <- c(met="circle", rxn="square")
edgeStyleMap <- c(main="solid", trans="dashed")

getDotNodeStyleAttributes <- function(attrs) {
    logPval <- if (!is.null(attrs$logPval)) attrs$logPval else 1
    with(attrs, data.frame(
        label=if (!is.null(attrs$label)) label else "",
        shape=if (!is.null(attrs$nodeType)) nodeShapeMap[nodeType] else "circle",
        fixedsize="true",
        style="filled",
        width=sapply(logPval, getDotSize),
        fontsize=sapply(logPval, getDotSize) * 45,
        color=if (!is.null(attrs$log2FC)) sapply(log2FC, getDotColor) else "black",
        fillcolor=if (!is.null(attrs$log2FC)) sapply(log2FC, getDotColor) else "white"
    ))
}

getDotEdgeStyleAttributes <- function(attrs) {
    logPval <- if (!is.null(attrs$logPval)) attrs$logPval else 1
    with(attrs, data.frame(
        label=if (!is.null(attrs$label)) label else "",
        style=if (!is.null(attrs$rptype)) edgeStyleMap[rptype] else "solid",
        penwidth=sapply(logPval, getDotSize) * 20,
        fontsize=sapply(logPval, getDotSize) * 45,
        color=if (!is.null(attrs$log2FC)) sapply(log2FC, getDotColor) else "grey"
    ))
}

getDotTooltip <- function(attr.values) {
    attr.strings <- list()
    attr.names <- names(attr.values)
    for(i in seq_along(attr.values)) {
        attr.strings[[i]] <- sprintf("%s: %s", attr.names[i], attr.values[[i]])
    }
    names(attr.strings) <- attr.names
    apply(do.call("cbind", attr.strings), 1, paste0, collapse="&#10;")
}

getNodeDotStrings <- function(module, indent="") {
    if (length(V(module)) == 0) {
        return(NULL)
    }
    attr.values <- as_data_frame(module, what="vertices")
    style.attr.values <- getDotNodeStyleAttributes(attr.values)
    # ignoring technical nodeType and big pathway attributes
    tooltip <- getDotTooltip(attr.values[, !colnames(attr.values) %in% c("pathway", "nodeType")])
    url <- attr.values$url

    attr.dotStrings <- getAttrDotStrings(cbind(style.attr.values, tooltip=tooltip))
    if (length(url) != 0) {
        attr.dotStrings <- cbind(attr.dotStrings,
                                 getAttrDotStrings(data.frame(URL=url, target="_blank")))
    }


    if(is.null(V(module)$name))
    {
        V(module)$name <- as.character(V(module))
    }


    node.label <- V(module)$name
    node.id <- as.vector(V(module))
    node.attrs <- apply(attr.dotStrings, 1, function(x) paste(na.omit(x), collapse=", "))
    nodeStrings <- sprintf("%sn%s [ %s ];\n", indent, node.id, node.attrs)
    nodeStrings
}

getEdgeDotStrings <- function(module, indent="") {
    if (length(E(module)) == 0) {
        return(NULL)
    }
    attr.values <- as_data_frame(module, what="edges")
    style.attr.values <- getDotEdgeStyleAttributes(attr.values)

    # ignoring big pathway attribute
    tooltip.values <- attr.values[, !colnames(attr.values) %in% c("pathway")]
    # for readability replacing node IDs with labels
    tooltip.values$from <- V(module)[tooltip.values$from]$label
    tooltip.values$to <- V(module)[tooltip.values$to]$label
    tooltip <- getDotTooltip(tooltip.values)

    url <- attr.values$url


    attr.dotStrings <- getAttrDotStrings(cbind(style.attr.values, tooltip=tooltip, labeltooltip=tooltip))

    if (length(url) != 0) {
        attr.dotStrings <- cbind(attr.dotStrings,
                                 getAttrDotStrings(data.frame(URL=url, target="_blank")))
    }

    edgelist.names <- get.edgelist(module, names=TRUE)
    edgelist.names <- paste(edgelist.names[,1], edgelist.names[,2], sep=" (pp) ")
    edgelist.ids <- get.edgelist(module, names=FALSE)

    edge.attrs <- apply(attr.dotStrings, 1, function(x) paste(na.omit(x), collapse=", "))
    edgeStrings <- sprintf("%sn%s -- n%s [ %s ];\n",
                           indent,
                           edgelist.ids[,1],
                           edgelist.ids[,2],
                           edge.attrs)
    edgeStrings
}
