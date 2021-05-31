Unfolding <- function(A, ENV = NULL, TransAbund = "Gaussian Columns", offset = 0.5, weight = "All_1", Constrained = FALSE, TransEnv = "Standardize columns", InitConfig = "SVD", 
                      model = "Ratio", condition = "Columns", Algorithm="SMACOF", OptimMethod="CG", r = 2, maxiter = 100, tolerance = 1e-05, lambda = 1, omega = 0, plot = FALSE) {
  # offset is the quantity added to the zeros of the table
  models = c("Absolute", "Ratio", "Interval", "Ordinal")
  conditions = c("Matrix", "Columns")
  InitConfigs = c("SVD", "Random", "Correspondences", "PCAEnvironment")
  TransAbunds = c("None", "Gaussian", "Column Percent", "Gaussian Columns", "Inverse Square Root", "Divide by Column Maximum")
  weights = c("All_1", "0 absences")
  TransEnvs = c("Raw Data", "Column centering", "Standardize columns")
  Algorithms=c("SMACOF", "GD")
  
  SiteNames = rownames(A)
  SpeciesNames = colnames(A)
  
  # Arguments passed as numbers are converted into strings
  if (is.numeric(model)) 
    model = models[model]
  if (is.numeric(condition)) 
    condition = conditions[condition]
  if (is.numeric(InitConfig)) 
    InitConfig = InitConfigs[InitConfig]
  if (is.numeric(TransAbund)) 
    TransAbund = TransAbunds[TransAbund]
  if (is.numeric(weight)) 
    weight = weights[weight]
  if (is.numeric(TransEnv)) 
    TransEnv = TransEnvs[TransEnv]
  if (is.numeric(Algorithm)) 
    Algorithm = Algorithms[Algorithm]
  
  if (is.data.frame(A)) 
    A = as.matrix(A)
  
  rownames(A) <- SiteNames
  
  if (!is.null(ENV)) {
    if (is.data.frame(ENV)) 
      ENV = as.matrix(ENV)
    
    #Initial Weighted averages (before transformation of abundances)
    c = apply(A, 2, sum)
    W1 = diagonal(1/c) %*% t(A) %*% ENV
    rownames(W1) = SpeciesNames
    
    VarNames = colnames(ENV)
    # Useful Parameters of the environmental variables
    Minima_Z = apply(ENV, 2, min)
    Maxima_Z = apply(ENV, 2, max)
    Medians_Z = apply(ENV, 2, quantile, 0.5)
    P25_Z = apply(ENV, 2, quantile, 0.25)
    P75_Z = apply(ENV, 2, quantile, 0.75)
    Means = apply(ENV, 2, mean)
    Deviations = apply(ENV, 2, sd)
    ENV = InitialTransform(ENV, transform = TransEnv)$X
    
    rownames(ENV) = SiteNames
    colnames(ENV) = VarNames
  }
  
  n = dim(A)[1]
  m = dim(A)[2]
  
  # Transforming the abundance/preference values
  P = TransfUnfold(A, TransAbund, offset)
  # Calculating Weights
  
  if (weight == "All_1") {
    W = matrix(1, n, m)
  } else {
    W = matrix(as.numeric(A > 0), n, m)
  }
  
  # Calculating Initial Configurations
  if (InitConfig == "SVD") {
    B = -0.5 * (diag(n) - matrix(1, n, n)/(n)) %*% P^2 %*% (diag(m) - matrix(1, m, m)/(m))
    sdec = svd(B)
    X0 = sdec$u[, 1:r] %*% diag(sqrt(sdec$d[1:r]))
    Y0 = sdec$v[, 1:r] %*% diag(sqrt(sdec$d[1:r]))
  }
  
  if (InitConfig == "Random") {
    X0 = matrix(rnorm(n * r), n, r)
    media = ColMeans(X0)
    X0 = X0 - matrix(1, n, 1) %*% media
    Y0 = matrix(rnorm(m * r), m, r)
    media = ColMeans(Y0)
    Y0 = Y0 - matrix(1, m, 1) %*% media
  }
  
  if (InitConfig == "Correspondences") {
    CAS = CA(P)
    X0 = CAS$RowCoordinates
    Y0 = CAS$ColCoordinates
  }
  
  # Calculating Unfolding coordinates
  if (Algorithm=="SMACOF")
  Unfold = UnfoldingSMACOF(P, W = W, Constrained = Constrained, ENV = ENV, model = model, condition = condition, r = r, maxiter = maxiter, tolerance = tolerance, 
                           lambda = 1, omega = 0, X0 = X0, Y0 = Y0)
  
  if (Algorithm=="GD")
    Unfold = UnfoldingGD(P, W = W, Constrained = Constrained, ENV = ENV, model = model, condition = condition, r = r, maxiter = maxiter, tolerance = tolerance, 
                             lambda = 1, omega = 0, X0 = X0, Y0 = Y0, OptimMethod=OptimMethod)
    
  Unfold$call <- match.call()
  Unfold$Abundances = A
  Unfold$Offset = offset
  if (!is.null(ENV)) {
    # Useful Parameters of the environmental variables
    Unfold$Minima_Z = Minima_Z
    Unfold$Maxima_Z = Maxima_Z
    Unfold$Medians_Z = Medians_Z
    Unfold$P25_Z = P25_Z
    Unfold$P75_Z = P75_Z
    Unfold$Means_Z = Means
    Unfold$Deviations_Z = Deviations
    p = dim(ENV)[2]
    Unfold$Weighted_Averages = W1
  }
  
  Unfold$Initial = InitConfig
  
  if (!is.null(ENV)) {
    # Calculating environmental coordinates
    Unfold$InterSet = cor(ENV, Unfold$X)
    Unfold$Env_Var_Scores = t(ENV) %*% Unfold$X %*% solve(t(Unfold$X) %*% Unfold$X)
    Unfold$nvars = ncol(ENV)
    ESP = Unfold$X %*% t(Unfold$Env_Var_Scores)
    Unfold$Env_fit = sum(ESP^2)/sum(ENV^2)
    Unfold$Var_Fit = matrix(1, p, 1)
    rownames(Unfold$Var_Fit) <- VarNames

    Unfold$QualityVars=cor(ENV, Unfold$X)^2
    for (j in 1:p) Unfold$Var_Fit[j] = sum(ESP[, j]^2)/sum(ENV[, j]^2)
  }
  
  if (plot) 
    plot(Unfold, mode = "s")
  return(Unfold)
}