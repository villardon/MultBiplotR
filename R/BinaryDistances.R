# Autor: Jose Luis Vicente Villardon
# Dpto. de Estadistica
# Universidad de Salamanca
# Revisado: Noviembre/2013

BinaryDistances <- function(x, y=NULL, coefficient= "Simple_Matching", transformation="sqrt(1-S)") {
  if (!is.matrix(x)) stop("Input must be a matrix")
  if (!CheckBinaryMatrix(x)) stop("Input must be a binary matrix (with 0 or 1 values)")
  coefficients = c("Kulezynski", "Russell_and_Rao", "Jaccard", "Simple_Matching", "Anderberg", "Rogers_and_Tanimoto", "Sorensen_Dice_and_Czekanowski", 
                   "Sneath_and_Sokal", "Hamman", "Kulezynski2", "Anderberg2", "Ochiai", "S13", "Pearson_phi", "Yule", "Sorensen", "Dice")
  if (is.numeric(coefficient)) coefficient=coefficients[coefficient]
  if (is.null(y)) y=x
	a = y %*% t(x)
	b = y %*% t(1 - x)
	c = (1 - y) %*% t(x)
	d = (1 - y) %*% t(1 - x)
	switch(coefficient, Kulezynski = {
		sim = a/(b + c)
	}, Russell_and_Rao = {
		sim = a/(a + b + c+d)
	}, Jaccard = {
		sim = a/(a + b + c)
	}, Simple_Matching = {
		sim = (a + d)/(a + b + c + d)
	}, Anderberg = {
		sim = a/(a + 2 * (b + c))
	}, Rogers_and_Tanimoto = {
		sim = (a + d)/(a + 2 * (b + c) + d)
	}, Sorensen_Dice_and_Czekanowski = {
		sim = a/(a + 0.5 * (b + c))
	}, Sneath_and_Sokal = {
		sim = (a + d)/(a + 0.5 * (b + c) + d)
	}, Hamman = {
		sim = (a - (b + c) + d)/(a + b + c + d)
	}, Kulezynski = {
		sim = 0.5 * ((a/(a + b)) + (a/(a + c)))
	}, Anderberg2 = {
		sim = 0.25 * (a/(a + b) + a/(a + c) + d/(c + d) + d/(b + d))
	}, Ochiai = {
		sim = a/sqrt((a + b) * (a + c))
	}, S13 = {
		sim = (a * d)/sqrt((a + b) * (a + c) * (d + b) * (d + c))
	}, Pearson_phi = {
		sim = (a * d - b * c)/sqrt((a + b) * (a + c) * (d + b) * (d + c))
	}, Yule = {
		sim = (a * d - b * c)/(a * d + b * c)
	}, Sorensen = {
	  sim = (2*a)/(2* a  + b + c)
	}, Dice = {
	  sim = (2*a)/(2* a  + b + c)
	})
  
  transformations= c("Identity", "1-S", "sqrt(1-S)", "-log(s)", "1/S-1", "sqrt(2(1-S))", "1-(S+1)/2", "1-abs(S)", "1/(S+1)")
  if (is.numeric(transformation)) transformation=transformations[transformation]
  
  switch(transformation, `Identity` = {
    dis=sim
  }, `Identity` = {
    dis=sim
  }, `1-S` = {
    dis=1-sim
  }, `sqrt(1-S)` = {
    dis = sqrt(1 - sim)
  }, `-log(s)` = {
    dis=-1*log(sim)
  }, `1/S-1` = {
    dis=1/sim -1 
  }, `sqrt(2(1-S))` = {
    dis== sqrt(2*(1 - sim))
  }, `1-(S+1)/2` = {
    dis=1-(sim+1)/2
  }, `1-abs(S)` = {
    dis=1-abs(sim)
  }, `1/(S+1)` = {
    dis=1/(sim)+1
  })


  return(dis)
}

