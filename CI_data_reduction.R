### CI AND DATA REDUCTION ###

## Function reduces dimensionality of input data by SVD, returns U matrix that explains critical amonut of variance
reduce <- function(x, crit=0.95) {
  decomp <- svd(x)
  cum_var <- cumsum(decomp$d^2)/sum(decomp$d^2)
  npc <- length(cum_var[cum_var<crit])+1
  message("explained variance:")
  print(cum_var[npc])
  return(decomp$u[,1:npc, drop=FALSE])
}


## Function inserts pKa which is varied inside CI_ML on right position in pKa vector
## works even if there's only 1 pKa
insert <- function(x1,x2=NULL,pos) {
  x <- numeric((length(x1)+length(x2)))
  x[pos] <- x1
  x[-pos] <- x2
  return(x)
}

## Function calculates lower 95%CI of Cost function, arguments are modified (varied and fixed pKa) for CI_ML
Cost_pKa <- function (pK_var, pK_fix, pos, x, pH, crit, ctot=1) {
  
  ## Varied (pK_var) and fixed (pK_fix) are merged in original pK vector
  pK <- insert(pK_var,pK_fix,pos)
  npk <- length(pK)
  y <- conc_matrix(npk, pH, pK=pK, ctot=ctot)
  y_org <- y$cmat
  pK <- y$pK
  H <- 10^-pH
  K <- 10^-pK
  y_scale <- scale(y_org)
  
  ## J is lower limit of 95%CI of Cost function, calculated from Chi^2 distribution
  J <- length(pH)*Cost_Function(x,y_scale)/qchisq(0.975,length(pH))-crit
  return(J)
}



## Function CI_ML calculates 95% CI of each pKa value calculated by ML_pKa function
## The only argument of a function is the output list from ML_pKa
## 95%CI is calculated from Chi^2 distribution of each point of Cost function
CI_ML <- function(data_ML) {
  x <- data_ML$x
  pH <- data_ML$pH
  pK <- data_ML$pK
  crit <- length(pH)*data_ML$J [length(data_ML$J)]/qchisq(0.025,length(pH))
  CI <- matrix(nrow = 2, ncol=length(pK))
  
  ## tryCatch ensures that if uniroot can't find root of a function, function returns Inf for that CI limit
  ## uniroot looks for lower 95% confidence limit of pKa in range 0-pKa and upper in range pKa-14
  for (i in 1:length(pK)) {
    CI[1,i] <- tryCatch(uniroot(Cost_pKa, pK_fix=pK[-i], pos=i, x=x, pH=pH, crit=crit, interval= c(0, pK[i]))$root, 
                        error= function(e){return(-Inf)})
    CI[2,i] <- tryCatch(uniroot(Cost_pKa, pK_fix=pK[-i], pos=i, x=x, pH=pH, crit=crit, interval= c(pK[i], 14))$root, 
                        error= function(e){return(Inf)})
  }
  
  CI <- data.frame(CI, row.names = c("lower 95% CL", "upper 95% CL"))
  names(CI) <- paste("pKa", seq_along(pK), sep="")
  
  return(CI)
}

