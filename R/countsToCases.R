countsToCases <- function(x, countcol = "Freq") {
  # Get the row indices to pull from x
  idx <- rep.int(seq_len(nrow(x)), x[[countcol]])
  # Drop count column
  x[[countcol]] <- NULL
  # Get the rows from x
  x[idx, ]
}