.messagef <- function (...)  { message(sprintf(...)) }
.warningf <- function (...)  { warning(sprintf(...)) }

.replaceNA <- function(x, y) { ifelse(is.na(x), y, x) }
