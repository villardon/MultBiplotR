# cbreaks <- function (range, breaks = extended_breaks(), labels = scientific_format()) 
# {
#   if (zero_range(range)) {
#     return(list(breaks = range[1], labels = format(range[1])))
#   }
#   if (is.function(breaks)) {
#     breaks <- breaks(range)
#     if (!is.function(labels)) {
#       stop("Labels can only be manually specified in conjunction with breaks", 
#            call. = FALSE)
#     }
#   }
#   if (is.function(labels)) {
#     labels <- labels(breaks)
#   }
#   else {
#     if (length(labels) != length(breaks)) {
#       stop("Labels and breaks must be same length")
#     }
#     if (is.expression(labels)) {
#       labels <- as.list(labels)
#     }
#     else {
#       labels <- as.character(labels)
#     }
#   }
#   list(breaks = breaks, labels = labels)
# }
# 
# 
# zero_range <- function (x, tol = 2.220446e-14) 
# {
#   if (length(x) == 1) 
#     return(TRUE)
#   if (length(x) != 2) 
#     stop("x must be length 1 or 2")
#   if (any(is.na(x))) 
#     return(NA)
#   if (x[1] == x[2]) 
#     return(TRUE)
#   if (all(is.infinite(x))) 
#     return(FALSE)
#   m <- min(abs(x))
#   if (m == 0) 
#     return(FALSE)
#   abs((x[1] - x[2])/m) < tol
# }
# 
# 
# pretty_breaks <-function (n = 5, ...) 
# {
#   function(x) {
#     breaks <- pretty(x, n, ...)
#     names(breaks) <- attr(breaks, "labels")
#     breaks
#   }
# }
# 
# extended_breaks <- function (n = 5, ...) 
# {
#   function(x) {
#     labeling::extended(min(x), max(x), n, only.loose = FALSE, 
#                        ...)
#   }
# }
# 
# scientific_format <- function (digits = 3, ...) 
# {
#   function(x) scientific(x, digits, ...)
# }
# 
# scientific <-function (x, digits = 3, ...) 
# {
#   x <- signif(x, digits)
#   format(x, trim = TRUE, scientific = TRUE, ...)
# }