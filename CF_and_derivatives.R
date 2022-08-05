### COST FUNCTION AND DERIVATIVES ###
Cost_Function <- function(X,y) {
  q=ncol(y)
  n=nrow(X)
  m=ncol(X)
  Hat=X%*%t(X)
  J <- 0
  for (i in 1:q) {
    for (j in 1:n) {
      J <- J+( (y[j,i]-X[j,]%*%t(X)%*%y[,i])/(1-Hat[j,j]) )^2
    }
  }
  return(J/n/q)
}

dJ_dy <- function (X,y) {
  q=ncol(y)
  n=nrow(X)
  m=ncol(X)
  solution= matrix (1, nrow=n, ncol=q) 
  
  diffy = function(X,y,J_root,n) {
    Hat=X%*%t(X)
    dy=rep(1,n)/(1-diag(Hat))*J_root
    for (i in 1:n) {
      for (j in 1:n) {
        dy[i]=dy[i]-(X[j,]%*%t(X[i,,drop=FALSE]))/(1-Hat[j,j])*J_root[j]
      }
      
    }
    return(2*dy)
  }
  
  J_root = function(X,y,n) {
    y=as.matrix(y)
    J=vector(length=n)
    Hat=X%*%t(X)
    for (j in 1:n) {
      J[j]=(y[j,1]-X[j,]%*%t(X)%*%y[,1])/(1-Hat[j,j]) 
    }
    return(J)
  }
  
  for (i in 1:q) {
    solution[,i]= diffy(X,y[,i],J_root(X,y[,i],n),n)
  }
  return(solution/n/q)
}

dy_dyorg <- function(yorg){
  q=ncol(yorg)
  n=nrow(yorg)
  solution = array (1, dim=c(n,n,q))
  for (j in 1:q) {
    for (i in 1:n) {
      for (k in 1:n) {
        if (i==k) solution[k,i,j]=-(n-1)^(-1)*sd(yorg[,j])^-3 *(yorg[i,j]-mean(yorg[,j]))^2 + sd(yorg[,j])^-1 *(1-1/n)
        if (i!=k) solution[k,i,j]=-(n-1)^(-1)*sd(yorg[,j])^-3 *(yorg[i,j]-mean(yorg[,j]))*(yorg[k,j]-mean(yorg[,j])) + sd(yorg[,j])^-1 *(-1/n)
      }
    }
  }
  return(solution)
}

dyorg_dK <- function(H, K) {
  n <- length(H)
  npk <- length(K)
  K <- matrix(rep(K,n), byrow=TRUE, nrow=n)
  derivatives <- array(dim=c(n,(npk+1),npk))
  brojnik <- matrix(nrow=n,ncol=(npk+1))
  nazivnik <- vector("numeric", n)
  for (i in 1:(npk+1)) {
    fac <- H^(i-1)
    if (i!=npk+1) {
      for (j in 1:(npk+1-i)) {
        fac <- fac*K[,j]
      }
    }
    brojnik[,(npk+2-i)] <- fac 
    nazivnik <- nazivnik+fac
  }
  der_naz <- array(dim=c(n,(npk+1),npk))
  for (i in 1:npk) {
    der_naz [,2:(npk+1),i] <- brojnik [,2:(npk+1)]
  }
  for (i in 1:npk) {
    for (j in 1:(npk+1)) {
      if (i>=j) {
        der_naz [,j,i] <- 0
      }
    }
  }
  for (i in 1:(npk+1)){
    der_naz[,i,] <- der_naz[,i,]/K
  }
  for (i in 1:(npk+1)) {
    for (j in 1:npk) {
      derivatives [,i,j] <- der_naz[,i,j]*nazivnik^-1-brojnik[,i]*nazivnik^-2*rowSums(der_naz[,,j])
    }
  }
  return(derivatives)
}

dK_dpK <- function(K) {
  solution = -log(10)*K
}

derivative_composition <- function (..., pos_pK) {
  derivatives <- list(...)
  q=ncol(derivatives[[1]])
  n=nrow(derivatives[[1]])
  derivative= matrix(1,nrow=n,ncol=q)
  for (i in 1:q) {
    derivative [,i] = derivatives[[1]] [,i] %*% derivatives[[2]] [,,i]
  }
  temp_derivative=0
  for (i in 1:q) {
    temp_derivative <- temp_derivative+t(derivatives[[3]] [,i,pos_pK])%*%derivative[,i]
  }
  derivative <- temp_derivative
  if (length(derivatives)>3) {
    for (i in 4:length(derivatives)) {
      derivative <- derivative*derivatives[[i]]
    }
  }
  return (derivative)
}

grad_check <- function(x,pH,pK,ctot=1,eps=0.004) {
  n <- nrow(x)
  npk <- length(pK)
  checks <- numeric(npk)
  for (i in 1:npk) {
    pK1 <- pK
    pK2 <- pK
    pK1[i] <- pK[i]-eps
    pK2[i] <- pK[i]+eps
    y1 <- conc_matrix(npk, pH, pK=pK1, ctot=ctot)
    y2 <- conc_matrix(npk, pH, pK=pK2, ctot=ctot)
    y_org1 <- y1$cmat
    y_org2 <- y2$cmat
    y_scale1 <- scale(y_org1)
    y_scale2 <- scale(y_org2)
    checks [i] <- (Cost_Function(x,y_scale2)-Cost_Function(x,y_scale1))/2/eps
  }
  return(checks)
}
