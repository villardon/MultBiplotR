\name{plot.Binary.Logistic.Biplot}
\alias{plot.Binary.Logistic.Biplot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plots the results of a binary Logistic  Biplot
}
Plots the results of a binary Logistic  Biplot
}
\usage{
plot.Binary.Logistic.Biplot(x, F1 = 1, F2 = 2, ShowAxis = FALSE, margin = 0, 
PlotVars = TRUE, PlotInd = TRUE, WhatRows = NULL, WhatCols = NULL, 
LabelRows = TRUE, LabelCols = TRUE, ShowBox = FALSE, RowLabels = NULL, 
ColLabels = NULL, RowColors = NULL, ColColors = NULL, Mode = "s", 
TickLength = 0.01, RowCex = 0.8, ColCex = 0.8, SmartLabels = FALSE, 
MinQualityRows = 0, MinQualityCols = 0, dp = 0, PredPoints = 0, 
SizeQualRows = FALSE, SizeQualCols = FALSE, ColorQualRows = FALSE, 
ColorQualCols = FALSE, PchRows = NULL, PchCols = NULL, PlotClus = FALSE, 
TypeClus = "ch", ClustConf = 1, Significant = TRUE, alpha = 0.05, 
Bonferroni = TRUE, PlotSupVars = TRUE, AbbreviateLabels = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
An object of class Binary.Logistic.Biplot
}
  \item{F1}{
Dimension for the first axis of the representation. Default = 1
}
  \item{F2}{
Dimension for the second axis of the representation. Default = 2
}
  \item{ShowAxis}{
Should the axis of the representation be shown?
}
  \item{margin}{
Margin of the plot as a percentage. It gets some space for the labels.
}
  \item{PlotVars}{
Should the variables be plotted?
}
  \item{PlotInd}{
Should the individuals be plotted?
}
  \item{WhatRows}{
What Rows should be plotted. A binary vector containing which rows (individuals) should be plotted (1) and which should not (0).
}
  \item{WhatCols}{
What Columns should be plotted. A binary vector containing which columns (variables) should be plotted (1) and which should not (0).
}
  \item{LabelRows}{
Should the individuals be labeled?
}
  \item{LabelCols}{
Should the individuals be labeled?
}
  \item{ShowBox}{
Should a box around the points be plotted?
}
  \item{RowLabels}{
A vector of row labels. If NULL the labels contained in the object will be used.
}
  \item{ColLabels}{
A vector of column labels. If NULL the labels contained in the object will be used.
}
  \item{RowColors}{
A vector of alternative row colors. 
}
  \item{ColColors}{
A vector of alternative column colors. 
}
  \item{Mode}{
Mode of the biplot: "p", "a", "b", "h", "ah" and "s".
}
  \item{TickLength}{
Length of the scale ticks for the biplot variables.
}
  \item{RowCex}{
Cex (Size) of the rows (marks and labels). Can be a single common size for all the points or a vector with individual sizes.
}
  \item{ColCex}{
Cex (Size) of the columns (marks and labels). Can be a single common size for all the points or a vector with individual sizes.
}
  \item{SmartLabels}{
Should the labels be placed in a smart way?
}
  \item{MinQualityRows}{
Minimum quality of the rows to be plotted. (Between 0 and 1)
}
  \item{MinQualityCols}{
Minimum quality of the columns to be plotted. (Between 0 and 1)
}
  \item{dp}{
A vector of variable indices to project all the individuals onto each variable of the vector.
}
  \item{PredPoints}{
A vector of row indices to project onto each variable.
}
  \item{SizeQualRows}{
Should the size of the Row points be related to its quality?
}
  \item{SizeQualCols}{
Should the size of the Column points be related to its quality?
}
  \item{ColorQualRows}{
Should the color of the Row points be related to its quality?
}
  \item{ColorQualCols}{
Should the color of the Column points be related to its quality?
}
  \item{PchRows}{
Marks for the rows (numbers). Can be a single common mark for all the points or a vector with individual marks.
}
  \item{PchCols}{
Marks for the columns (numbers). Can be a single common mark for all the points or a vector with individual marks.
}
  \item{PlotClus}{
Should the clusters be plotted?
}
  \item{TypeClus}{
Type of plot for the clusters. ("ch"- Convex Hull, "el"- Ellipse or "st"- Star)
}
  \item{ClustConf}{
Percent of points included in the cluster. only the ClusConf percent of the points nearest to the center will be used to calculate the cluster
}
  \item{Significant}{
Should only the significant variables be plotted?
}
  \item{alpha}{
Signification level.
}
  \item{Bonferroni}{
Should the Bonferroni correction be used?
}
  \item{PlotSupVars}{
Should the Supplementary variables be plotted?
}
  \item{AbbreviateLabels}{
Should labels be abbreviated?
}
  \item{\dots}{
Any other graphical parameter.
}
}
\details{
Plots a biplot for binary data. The Biplot for binary data is taken as the basis of the plot. If there are a mixture of different types of variables (binary, nominal, abundance, ...) are added to the biplot as supplementary parts.

There are several modes for plotting the biplot.
"p".- Points (Rows and Columns are represented by points)

"a" .- Arrows (The traditional representation with points for rows and arrows for columns)

"b" .- The arrows for the columns are extended to both extremes of the plot and labeled outside the plot area.

"h" .- The arrows for the columns are extended to the positive extreme of the plot and labeled outside the plot area.

"ah" .- Same as arrows but labeled outside the plot area.

"s" .- The directions (or biplot axes) have a graded scale for prediction of the original values.
}
\value{
The plot of the biplot.
}
\references{
Vicente-Villardon, J. L., Galindo, M. P. and Blazquez, A. (2006) Logistic Biplots. In Multiple Correspondence Análisis And Related Methods. Grenacre, M & Blasius, J, Eds,  Chapman and Hall, Boca Raton.

Demey, J., Vicente-Villardon, J. L., Galindo, M.P.  AND Zambrano, A. (2008) Identifying Molecular Markers Associated With Classification Of Genotypes Using External Logistic Biplots. Bioinformatics, 24(24): 2832-2838.
}
\author{
Jose Luis Vicente Villardon
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (x, F1 = 1, F2 = 2, ShowAxis = FALSE, margin = 0, PlotVars = TRUE, 
    PlotInd = TRUE, WhatRows = NULL, WhatCols = NULL, LabelRows = TRUE, 
    LabelCols = TRUE, ShowBox = FALSE, RowLabels = NULL, ColLabels = NULL, 
    RowColors = NULL, ColColors = NULL, Mode = "s", TickLength = 0.01, 
    RowCex = 0.8, ColCex = 0.8, SmartLabels = FALSE, MinQualityRows = 0, 
    MinQualityCols = 0, dp = 0, PredPoints = 0, SizeQualRows = FALSE, 
    SizeQualCols = FALSE, ColorQualRows = FALSE, ColorQualCols = FALSE, 
    PchRows = NULL, PchCols = NULL, PlotClus = FALSE, TypeClus = "ch", 
    ClustConf = 1, Significant = TRUE, alpha = 0.05, Bonferroni = TRUE, 
    PlotSupVars = TRUE, ...) 
{
    a = x$RowCoordinates[, c(F1, F2)]
    n = dim(a)[1]
    neje = dim(a)[2]
    p = dim(x$ColumnParameters)[1]
    if (is.null(WhatRows)) 
        WhatRows = matrix(1, n, 1)
    else if (!CheckBinaryVector(WhatRows)) {
        AllRows = matrix(0, n, 1)
        AllRows[WhatRows] = 1
        WhatRows = AllRows
    }
    WhatRows = as.logical(WhatRows)
    RowQualities = (x$RowQualities[, F1] + x$RowQualities[, F2])/100
    if (is.null(WhatCols)) 
        WhatCols = matrix(1, p, 1)
    else if (!CheckBinaryVector(WhatCols)) {
        AllCols = matrix(0, p, 1)
        AllCols[WhatCols] = 1
        WhatCols = AllCols
    }
    WhatCols = as.logical(WhatCols)
    if (Significant) {
        if (Bonferroni) 
            WhatCols = WhatCols & ((x$pvalues * p) < alpha)
        else WhatCols = WhatCols & (x$pvalues < alpha)
    }
    WhatCols = WhatCols & (x$R2 > MinQualityCols)
    WhatCols[which((is.na(WhatCols)))] = FALSE
    if (is.null(RowColors)) 
        RowColors = matrix("blue", n, 1)
    if (length(RowColors) == 1) 
        RowColors = rep(RowColors, n)
    if (is.null(ColColors)) 
        ColColors = matrix("black", p, 1)
    if (is.null(RowLabels)) 
        RowLabels = rownames(x$RowCoordinates)
    if (is.null(ColLabels)) 
        ColLabels = rownames(x$ColumnParameters)
    if (length(RowCex == 1)) 
        RowCex = rep(RowCex, n)
    if (length(ColCex == 1)) 
        ColCex = rep(ColCex, p)
    qlrcols = x$VarInfo$Nagelkerke
    qlrrows = x$RowQualities[, F1] + x$RowQualities[, F2]
    if (SizeQualRows) 
        RowCex = cscale(qlrrows, rescale_pal())
    if (SizeQualCols) 
        ColCex = cscale(qlrcols, rescale_pal())
    if (ColorQualRows) 
        RowColors = cscale(qlrrows, seq_gradient_pal("white", 
            "red"))
    if (ColorQualCols) 
        ColColors = cscale(qlrcols, seq_gradient_pal("white", 
            "blue"))
    if (ShowAxis) {
        xaxt = "s"
        yaxt = "s"
    }
    else {
        xaxt = "n"
        yaxt = "n"
    }
    if ((margin < 0) | (margin > 0.5)) 
        margin = 0
    P = a
    xmin = min(P[, 1])
    xmax = max(P[, 1])
    ymin = min(P[, 2])
    ymax = max(P[, 2])
    xrang = abs(xmax - xmin)
    yrang = abs(ymax - ymin)
    if (xmax < 0) 
        xmax = xmax * (-1)
    print(c(xmin, xmax, ymin, ymax))
    P = rbind(P, c(xmax + (xmax - xmin) * margin, ymax + (ymax - 
        ymin) * margin), c(xmin - (xmax - xmin) * margin, ymin - 
        (ymax - ymin) * margin))
    plot(P[, 1], P[, 2], asp = 1, xaxt = xaxt, yaxt = yaxt, cex = 0, 
        bty = "n", xlab = paste("Dimension", F1), ylab = paste("Dimension", 
            F2), main = x$Biplot, ...)
    if (ShowBox) 
        rect(xmin, ymin, xmax, ymax)
    if (PlotClus) {
        RowColors = PlotBiplotClusters(a, x$Clusters, TypeClus = "ch", 
            ClusterColors = x$ClusterColors, ClusterNames = x$ClusterNames, 
            centers = TRUE, ClustConf = ClustConf, ...)
    }
    if (PlotInd) {
        for (i in 1:n) if (WhatRows[i]) {
            points(a[i, 1], a[i, 2], col = RowColors[i], cex = RowCex[i], 
                pch = PchRows, ...)
            if (LabelRows) 
                text(a[i, 1], a[i, 2], labels = RowLabels[i], 
                  col = RowColors[i], pos = 1, cex = RowCex[i], 
                  ...)
        }
    }
    if (PlotVars) {
        for (i in 1:p) {
            if (WhatCols[i]) {
                PlotBinaryVar(b0 = x$ColumnParameters[i, 1], 
                  bi1 = x$ColumnParameters[i, F1 + 1], bi2 = x$ColumnParameters[i, 
                    F2 + 1], xmin = xmin, xmax = xmax, ymin = ymin, 
                  ymax = ymax, mode = Mode, Color = ColColors[i], 
                  label = ColLabels[i], tl = TickLength, CexPoint = ColCex[i])
            }
        }
    }
    for (idp in dp) {
        if ((idp > 0) & (idp < (p + 1))) 
            if (WhatCols[idp]) {
                g = matrix(c(x$ColumnParameters[idp, F1 + 1], 
                  x$ColumnParameters[idp, F2 + 1]), 2, 1)
                nn = (t(g) \%*\% g)
                scal <- (a \%*\% g)/nn[1, 1]
                Dscal <- diag(as.vector(scal))
                Fpr <- Dscal \%*\% matrix(rep(1, nrow(a)), ncol = 1) \%*\% 
                  t(g)
                nrFpr <- nrow(Fpr)
                dlines(a, Fpr, color = ColColors[idp])
            }
    }
    for (idp in PredPoints) {
        if ((idp > 0) & (idp < (n + 1))) {
            for (j in 1:p) {
                if (WhatCols[j]) {
                  g = matrix(c(x$ColumnParameters[j, F1 + 1], 
                    x$ColumnParameters[j, F2 + 1]), 2, 1)
                  nn = (t(g) \%*\% g)
                  scal <- (a[idp, ] \%*\% g)/nn[1, 1]
                  Fpr <- scal \%*\% t(g)
                  nrFpr <- nrow(Fpr)
                  dlines(matrix(a[idp, ], 1, 2), Fpr, color = ColColors[j])
                }
            }
        }
    }
    if (PlotSupVars) 
        plot.Supplementary.Variables(x, F1 = F1, F2 = F2, xmin = xmin, 
            xmax = xmax, ymin = ymin, ymax = ymax, mode = Mode, 
            ColorSupContVars = "red")
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
