.messagef <- function (...)  { message(sprintf(...)) }
.warningf <- function (...)  { warning(sprintf(...)) }

.replaceNA <- function(x, y) { ifelse(is.na(x), y, x) }

.pairwiseCompare <- function (FUN, list1, list2 = list1, ...)
{
    additionalArguments <- list(...)
    f1 <- function(...) {
        mapply(FUN = function(x, y) {
            do.call(FUN, c(list(list1[[x]], list2[[y]]), additionalArguments))
        }, ...)
    }
    z <- outer(seq_along(list1), seq_along(list2), FUN = f1)
    rownames(z) <- names(list1)
    colnames(z) <- names(list2)
    z
}
