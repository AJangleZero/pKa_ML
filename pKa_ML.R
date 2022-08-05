### MAIN PROGRAM ML pKa ###

ML_pka <- function (x, pH, npK, pK=NA,
                    ctot=1, alpha=0.005, crit=10^-3, n_it=1000,...) {
  if (anyNA(pK)==TRUE) {pK <- rnorm(npK, min(pH)+(max(pH)-min(pH))/(npK+1)*seq(1,npK),(max(pH)-min(pH))/(npK+1)/4)}
  names(pK) <- paste("pKa", seq_along(pK), sep="")
  J_vec <- numeric(n_it)
  J_old <- Inf
  derivation <- numeric(npK)
  message("Initial pKa")
  print(pK)
  
  
  ## Defining matrices used for gradient checking
  {
    der1 <- matrix(nrow=n_it,ncol=npK)
    der2 <- matrix(nrow=n_it,ncol=npK)
  }
  
  
  for (i in 1:n_it) {
    y <- conc_matrix(npK, pH, pK=pK, ctot=ctot)
    y_org <- y$cmat
    pK <- y$pK
    H <- 10^-pH
    K <- 10^-pK
    y_scale <- scale(y_org)
    J <- Cost_Function(x,y_scale)
    J_vec[i] <- J
    
    ## Analytical derivative calculation 
    #for (j in 1:npK) {
    #  derivation [j]<- derivative_composition(dJ_dy(x,y_scale), dy_dyorg(y_org), 
    #    dyorg_dK(H,K), dK_dpK(K[j]), pos_pK = j) 
    #}
    
    #der1 represents analytical derivative, der2 is numerical derivative via gradient check
    {
    #  der1 [i,] <- derivation
      der2 [i,] <- grad_check(x,pH,pK,ctot,eps=10^-12)
    }
    
    ## pKa update
    pK=pK-alpha*der2[i,]
    message("Iteration: ",i)
    print(pK)
    
    
    ## If convergence is achieved, return the result
    if(abs(J_old-J)<crit) {return(list(pK=pK, J=J_vec[1:i], der1=der1[1:i,], der2=der2[1:i,],
                                       pH=pH, ctoc=ctot, x=x))}
    J_old <- J
  }
  message("no convergence achieved after ", n_it, " iterations")
  return(list(pK=pK, J=J_vec, der1=der2[1:i,], der2=der1[1:i,],
              pH=pH, ctoc=ctot, x=x))
}
