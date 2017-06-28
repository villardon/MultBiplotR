
RidgeOrdinalLogistic <- function(y, x, penalization = 0.1, tol = 1e-04, maxiter = 200, show = FALSE) {
	    if (is.matrix(x)) {
        n <- nrow(x)
    }
    else {
        n <- length(x)
    }
  y=as.numeric(y)
    model=OrdinalLogisticFit(y,x, penalization = penalization, tol = tol,  maxiter = maxiter, show = show)
    null=OrdinalLogisticFit(y,x=NULL, penalization = penalization, tol = tol,  maxiter = maxiter, show = show)
	
	
  model$DevianceNull = null$Deviance
	model$Dif = (model$DevianceNull - model$Deviance)
	model$df = model$nvar
	model$pval = 1 - pchisq(model$Dif, df = model$df)
	model$CoxSnell = 1 - exp(-1 * model$Dif/n)
	model$Nagelkerke = model$CoxSnell/(1 - exp((model$DevianceNull/(-2)))^(2/n))
	model$MacFaden = 1 - (model$Deviance/model$DevianceNull)
	
	class(model) = "OrdinalLogisticRegression"

	return(model)
	
	}
	