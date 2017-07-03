#' @export
sgmwcs.solver <- function (sgmwcs, nthreads = 1, timeLimit = -1,
                           nodes.group.by="signal",
                           edges.group.by="signal",
                           group.only.positive=T,
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
                          "--signals", instance$signals.file,
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
        res <- readGraph(node.file = solution.file,
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
                                group.only.positive=T) {



    dir.create(graph.dir, showWarnings = FALSE)
    edges.file <- file.path(graph.dir, "edges.txt")
    nodes.file <- file.path(graph.dir, "nodes.txt")
    signals.file <- file.path(graph.dir, "signals.txt")

    nt <- as_data_frame(network, what="vertices")
    if (!is.null(nodes.group.by)) {
        if (all(nodes.group.by %in% colnames(nt))) {
            nt$signal <- do.call(paste, c(nt[, nodes.group.by, drop=F], sep="\r"))
            if (group.only.positive) {
                nt$signal <- paste(nt$signal, ifelse(nt$score > 0,
                                                    "",
                                                    seq_len(nrow(nt))), sep="\r")
            }
        } else {
            .warningf("Can't collapse nodes, not all fields present: %s",
                     paste0(setdiff(all.vars(f), colnames(nt)), collapse=", "))
            nt$signal <- seq_len(nrow(nt))
        }

    } else {
        nt$signal <- seq_len(nrow(nt))
    }
    nt$signal <- factor(nt$signal)
    levels(nt$signal) <- paste0("ns", seq_len(length(levels(nt$signal))))


    et <- as_data_frame(network, what="edges")
    if (!is.null(edges.group.by)) {
        if (all(edges.group.by %in% colnames(et))) {
            et$signal <- do.call(paste, c(et[, edges.group.by, drop=F], sep="\r"))
            if (group.only.positive) {
                et$signal <- paste(et$signal, ifelse(et$score > 0,
                                                     "",
                                                     seq_len(nrow(et))), sep="\r")
            }
        } else {
            .warningf("Can't collapse edges, not all fields present: %s",
                      paste0(setdiff(all.vars(f), colnames(et)), collapse=", "))
            et$signal <- seq_len(nrow(et))
        }

    } else {
        et$signal <- seq_len(nrow(nt))

    }
    et$signal <- factor(et$signal)
    levels(et$signal) <- paste0("es", seq_len(length(levels(et$signal))))

    st <- rbind(nt[, c("signal", "score")],
                et[, c("signal", "score")])
    st <- st[!duplicated(st$signal),]
    rownames(st) <- NULL
    colnames(st) <- c("#signal", "score")

    nt <- nt[, c("name", "signal")]
    colnames(nt) <- c("#name", "signal")

    et <- et[, c("from", "to", "signal")]
    colnames(et) <- c("#from", "to", "signal")



    write.table(nt, file=nodes.file, sep="\t", row.names=F, quote=F, col.names=T)
    write.table(et, file=edges.file, sep="\t", row.names=F, quote=F, col.names=T)
    write.table(st, file=signals.file, sep="\t", row.names=F, quote=F, col.names=T)

    list(nodes.file=nodes.file,
         edges.file=edges.file,
         signals.file=signals.file)
}

#' @export
#' @importFrom methods newEmptyObject
solveSgmwcsRandHeur <- function(g,
                                nodes.group.by="signal",
                                edges.group.by="signal",
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

readGraph <- function(node.file, edge.file, network) {
    nodes <- readLines(node.file)

    if (length(nodes) <= 1) {
        return(induced.subgraph(network, vids = nodes))
    }

    edges <- as.matrix(read.table(file = edge.file,
                                  colClasses = "character"))

    eids <- get.edge.ids(network, t(edges))
    res <- subgraph.edges(network, eids = eids, delete.vertices = T)

    res
}
