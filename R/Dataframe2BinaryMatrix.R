# Autor: Jose Luis Vicente Villardon
# Dpto. de Estadistica
# Universidad de Salamanca
# Revisado: Noviembre/2013
# This fuction converts a data frame into a binary matrix
# Factors with two levels are converted to 0, 1
# factors with more tan two levels are converted into a binary
# data matrix with a many columns as factor levels
# and continuous variables are splitted using the mean, the median
# or  an user provided value as a cut point



Dataframe2BinaryMatrix <- function(dataf, cuttype="Median", cut =NULL, BinFact=TRUE) {
  if (!is.data.frame(dataf)) dataf=as.data.frame(dataf)
	if (is.data.frame(dataf)) {
		n = dim(dataf)[1]
		p = dim(dataf)[2]
		ColumnNames = colnames(dataf)
		RowNames = rownames(dataf)
		x = NULL
		for (j in 1:p) {
			if (is.factor(dataf[[j]])) 
				if ((length(levels(dataf[[j]])) == 2) & BinFact) {
					x1 = matrix(as.numeric(dataf[[j]]) - 1, n, 1)
					colnames(x1) <- ColumnNames[j]
					# colnames(x1) <- paste(ColumnNames[j], "-", levels(dataf[[j]])[2], sep = "")
				} else {
					x1 = Factor2Binary(dataf[[j]], ColumnNames[j])
				}
      
			if (is.ordered(dataf[[j]])) 
			  if (length(levels(dataf[[j]])) == 2) {
			    x1 = matrix(as.numeric(dataf[[j]]) - 1, n, 1)
			    colnames(x1) <- paste(ColumnNames[j], "-", levels(dataf[[j]])[2], sep = "")
			  } else {
			    x1 = Factor2Binary(dataf[[j]], ColumnNames[j])
			  }

			if (is.integer(dataf[[j]])) {
				x1 = Integer2Binary(dataf[[j]], ColumnNames[j])
			}
			if (is.numeric(dataf[[j]])) {
        if (is.null(cut))
        if (cuttype=="Median") cut1 =median(dataf[[j]])
        else cut1 =mean(dataf[[j]])
        else cut1=cut
				x1 = Numeric2Binary(dataf[[j]], ColumnNames[j], cut=cut1)
			}
			x = cbind(x, x1)
		}

		rownames(x) <- RowNames
		return(x)
	} else stop("Data is not organized as a data frame")
}
