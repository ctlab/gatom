convertPvalDT <- function(de.pvals, map) {
    from <- key(map)
    to <- setdiff(colnames(map), from)
    res.pvals <- map[de.pvals, nomatch=0]
    # setnames(res.pvals, c("ID", to), c(from, "ID"))
    res.pvals <- res.pvals[order(pval)]
    res.pvals <- res.pvals[!duplicated(res.pvals[, to, with=F])]
    setkeyv(res.pvals, to)
    setcolorder(res.pvals, c(to,
                             setdiff(colnames(res.pvals),
                                     c(to, from)),
                             from))
    res.pvals
}
