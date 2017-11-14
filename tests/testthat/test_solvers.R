context("Solvers")
library(igraph)

testSgmwcsInstance <- function() {
    et <- data.frame(
        from=  c("a1", "a2", "b" , "c" ),
        to=    c("b" , "b" , "c" , "d" ),
        signal=c("p1", "p1", "p2", "p3"),
        score=c(1,    1,    -1.5,  2  ),
        stringsAsFactors=F)
    res <- graph.data.frame(et, directed = FALSE)
    V(res)$signal <- seq_len(vcount(res))
    V(res)$score <- 0
    res
}


test_that("solveSgmwcsRandHeur works", {
    g <- testSgmwcsInstance()

    set.seed(42)
    m <- solveSgmwcsRandHeur(g)
    expect_false("a1" %in% V(m)$name)
    expect_true("d" %in% V(m)$name)
})

test_that("sgmwcs.solver works", {
    skip_if_not(system("sgmwcs -h", ignore.stdout = T) == 0)

    g <- testSgmwcsInstance()

    solve <- sgmwcs.solver("sgmwcs")
    m <- solve(g)
    expect_false("a1" %in% V(m)$name)
    expect_true("d" %in% V(m)$name)

})

test_that("sgmwcs.solver works when grouping by absent attributes", {
    skip_if_not(system("sgmwcs -h", ignore.stdout = T) == 0)

    g <- testSgmwcsInstance()

    solve <- sgmwcs.solver("sgmwcs", nodes.group.by = "bla")
    expect_warning(m <- solve(g))

    solve <- sgmwcs.solver("sgmwcs", edges.group.by = "bla")
    expect_warning(m <- solve(g))
})

test_that("sgmwcs.solver works for single-vertex solution", {
    skip_if_not(system("sgmwcs -h", ignore.stdout = T) == 0)

    g <- testSgmwcsInstance()
    V(g)$score <- c(1, rep(-1, vcount(g) - 1))
    E(g)$score <- 0

    solve <- sgmwcs.solver("sgmwcs")
    m <- solve(g)
    expect_equivalent(V(m)$name, V(g)[1]$name)
})
