.messagef <- function (...)  { message(sprintf(...)) }

.replaceNA <- function(x, y) { ifelse(is.na(x), y, x) }
