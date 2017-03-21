#' @export
sgmwcs.solver <- function (sgmwcs, nthreads = 1, timeLimit = -1, nodes.group.by=NULL, edges.group.by=NULL,
                           group.only.positive=F,
                           minimize.size=F,
                           other.args=NULL) {
    function(network) {
        network.orig <- network
        score.edges <- "score" %in% list.edge.attributes(network)
        score.nodes <- "score" %in% list.vertex.attributes(network)
        if (!score.nodes) {
            V(network)$score <- 0
        }
        if (!score.edges) {
            E(network)$score <- 0
        }

        graph.dir <- tempfile("graph")

        instance <- writeSgmwcsInstance(graph.dir = graph.dir,
                                        network=network,
                                        nodes.group.by = nodes.group.by,
                                        edges.group.by = edges.group.by,
                                        group.only.positive = group.only.positive
        )
        log.file <- file.path(graph.dir, "log")

        system2(sgmwcs, c("--nodes", instance$nodes.file,
                          "--edges", instance$edges.file,
                          "--synonyms", instance$synonyms.file,
                          "--threads", nthreads,
                          "--timelimit", timeLimit,
                          if (minimize.size) c("-p", 1e-3) else NULL,
                          other.args
        ), stdout = log.file, stderr = log.file)

        solution.file <- paste0(instance$nodes.file, ".out")
        if (!file.exists(solution.file)) {
            warning("Solution file not found")
            return(NULL)
        }
        res <- GAM:::readGraph(node.file = solution.file,
                               edge.file = paste0(instance$edges.file, ".out"),
                               network = network)
        solutionLog <- readLines(log.file)
        attr(res, "optimal") <- any(grepl("SOLVED TO OPTIMALITY", solutionLog))
        return(res)
    }
}

writeSgmwcsInstance <- function(graph.dir, network,
                                nodes.group.by=NULL,
                                edges.group.by=NULL,
                                group.only.positive=F) {



    dir.create(graph.dir, showWarnings = FALSE)
    edges.file <- file.path(graph.dir, "edges.txt")
    nodes.file <- file.path(graph.dir, "nodes.txt")
    synonyms.file <- file.path(graph.dir, "synonyms.txt")

    synonyms <- c()

    nt <- get.vertex.attributes(network)
    if (!is.null(nodes.group.by)) {
        f <- as.formula(sprintf("name ~ %s", nodes.group.by))
        if (all(all.vars(f) %in% colnames(nt))) {
            synonyms <- c(synonyms, aggregate(f, data=nt, paste0, collapse=" ")$name)
        } else {
            warningf("Can't collapse nodes, not all fields present: %s",
                     paste0(setdiff(all.vars(f), colnames(nt)), collapse=", "))
            synonyms <- c(synonyms, nt$name)
        }

    } else {
        synonyms <- c(synonyms, nt$name)
    }
    nt <- rename(nt[, c("name", "score")], c("name"="#name"))

    et <- get.edge.attributes(network, include.ends = T)
    if (!is.null(edges.group.by)) {
        etx <- if (group.only.positive) {
            synonyms <- c(synonyms, with(et[et$score <= 0,], sprintf("%s -- %s", from, to)))
            et[et$score > 0,]
        } else {
            et
        }
        if (nrow(etx) > 0) {
            synonyms <- c(synonyms,
                          aggregate(name ~ edges.group.by,
                                    data=list(
                                        name=sprintf("%s -- %s", etx$from, etx$to),
                                        edges.group.by=etx[[edges.group.by]]),
                                    paste0, collapse=" ")$name
            )
        }


    }
    et <- rename(et[, c("from", "to", "score")], c("from"="#from"))

    write.tsv(nt, file=nodes.file)
    write.tsv(et, file=edges.file)
    writeLines(sprintf("%s", synonyms), con=synonyms.file)

    if (length(synonyms) == 0) {
        synonyms.file <- NULL
    }

    list(nodes.file=nodes.file,
         edges.file=edges.file,
         synonyms.file=synonyms.file)
}

#' @export
solveSgmwcsRandHeur <- function(g,
                                nodes.group.by=NULL,
                                edges.group.by=NULL,
                                max.iterations = 10000) {
    n <- length(V(g))
    m <- length(E(g))
    vscores <- V(g)$score
    escores <- E(g)$score

    if (!is.null(nodes.group.by)) {
        vsignals <- ifelse(vscores > 0,
                           get.vertex.attribute(g, name=nodes.group.by),
                           seq_len(n))
    } else {
        vsignals <- seq_len(n)
    }

    if (!is.null(edges.group.by)) {
        esignals <- ifelse(escores > 0,
                           get.edge.attribute(g, name=edges.group.by),
                           seq_len(m))
    } else {
        esignals <- seq_len(m)
    }

    score <- function(nodes, edges) {
        nodes1 <- nodes[!duplicated(vsignals[nodes])]
        edges1 <- edges[!duplicated(esignals[edges])]
        sum(vscores[nodes1]) + sum(escores[edges1])
    }
    best.solutions <- list()
    best.solutions.edges <- list()
    edges <- get.edgelist(g, names = FALSE)
    solution <- list()
    solution$score <- -Inf
    makeSolution <- function(id, nodes, edges, score = NULL) {
        res <- newEmptyObject()
        res$id <- id
        res$nodes <- nodes
        res$edges <- edges
        res$score <- score
        if (is.null(res$score)) {
            res$score <- score(nodes, edges)
        }
        res
    }
    updateBest <- function(newsolution) {
        if (newsolution$score > solution$score) {
            solution <<- newsolution
        }
    }
    for (v in V(g)) {
        newsolution <- makeSolution(v, v, NULL)
        updateBest(newsolution)
        best.solutions[[v]] <- newsolution
    }
    best.scores <- vscores
    checked <- rep(FALSE, m)
    for (i in seq_len(max.iterations)) {
        edges.weights <- pmax(best.scores[edges[, 1]], best.scores[edges[,
                                                                         2]]) + escores
        edges.unchecked <- which(!checked)
        if (length(edges.unchecked) == 0) {
            break
        }
        candidates <- sapply(1:2, function(x) {
            max(ceiling(runif(1, max = length(edges.unchecked))),
                1)
        })
        edge <- candidates[which.max(edges.weights[candidates])]
        checked[edge] <- TRUE
        start <- edges[edge, 1]
        end <- edges[edge, 2]
        startsolution <- best.solutions[[start]]
        endsolution <- best.solutions[[end]]
        score.start <- startsolution$score
        score.end <- endsolution$score
        if (startsolution$id == endsolution$id && edge %in% startsolution$edges) {
            next
        }
        newnodes <- union(startsolution$nodes, endsolution$nodes)
        newedges <- union(edge, union(startsolution$edges, endsolution$edges))
        newsolution <- makeSolution(n + i, newnodes, newedges)
        to_update <- c()
        if (newsolution$score > score.start) {
            start.nodes <- startsolution$nodes
            to_update <- union(to_update, start.nodes[best.scores[start.nodes] <
                                                          newsolution$score])
        }
        if (newsolution$score > score.end) {
            end.nodes <- endsolution$nodes
            to_update <- union(to_update, end.nodes[best.scores[end.nodes] <
                                                        newsolution$score])
        }
        if (length(to_update) > 0) {
            updateBest(newsolution)
            best.solutions[to_update] <- list(newsolution)
            best.scores[to_update] <- newsolution$score
            to_check <- sapply(E(g)[adj(to_update)], identity)
            checked[to_check] <- FALSE
        }
    }
    solution
    .messagef("Score: %s", solution$score)
    subgraph.edges(g, E(g)[solution$edges])
}
