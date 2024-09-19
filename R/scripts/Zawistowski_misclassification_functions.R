#### R CODE FOR MISCLASSIFICATION-ADJUSTED ROC #####

#### FUNCTIONS FOR MISCLASSIFICATION ROC ####

## Implementation of the misclassification-adjusted logistic regression model
## described in Neuhaus (Biometrika 1999)
Misclassify_logistic_IWLS <- function(X, y, gamma0, gamma1, maxIter=300, tol=1E-6){
  # X is the covariate matrix
  # y is the response vector
  # maxIter is the maximum number of iterations
  # tol is a convergence criterion
  # gamma0 = P(Y=1 | T=0)
  # gamma1 = P(Y=0 | T=1)
  
  b <- bLast <- glm(y~X, family="binomial")$coef
  X <- cbind(1, X) # add constant
  
  w=0
  it <- 1 # iteration index
  while (it <= maxIter){
    eta <- X %*% b
    mu <- (1-gamma0-gamma1)*(1/(1 + exp(-eta)))+gamma0
    nu <- as.vector(mu*(1 - mu))
    w <- nu
    z <- log((mu-gamma0)/(1-gamma1-mu)) + (y - mu)/nu
    b <- lsfit(X, z, w, intercept=FALSE)$coef
    if (max(abs(b - bLast)/(abs(bLast) + 0.01*tol)) < tol) {break}
    bLast <- b
    it <- it + 1 # increment index
  }
  if(it > maxIter) {warning("maximum iterations exceeded"); list(coefficients=0, var=0, iterations=it)}
  else{ 
    Vb <- solve(t(X) %*% diag(w) %*% X)
    list(coefficients=b, var=Vb, iterations=it)
  }
}

# function to create misclassified version of simulated data
misclassify <- function(v,gamma0,gamma1) {
  p00 = 1-gamma0
  p11 = 1-gamma1
  miss_v=rep(-9, length(v)) 
  if(length(p00)==1) {p00 = rep(p00, length(v)); p11 = rep(p11, length(v))}
  for(i in 1:length(v)) {
    new_v=v[i]
    if(v[i]==0) {if(runif(1) > p00[i]) {new_v=1}}
    else{{if(runif(1) > p11[i]) {new_v=0}}}
    miss_v[i]=new_v
  }
  return(miss_v)
}

# Reverse misclassification function.
# Adaptation (Javier Mancilla-Galindo) of the 'misclassify' function, to 
# use reverse misclassification parameters instead of gamma0 and gamma1: 
misclassify_reverse <- function(v, rev_g0, rev_g1) {
  miss_v = rep(-9, length(v))
  for(i in 1:length(v)) {
    new_v = v[i]
    if(v[i] == 1) {  # If the CXR result is positive (1)
      # Misclassify based on the probability of HRCT+ given CXR+
      if(runif(1) > rev_g0) { new_v = 0 }  # Misclassify to HRCT-
    } else {  # If the CXR result is negative (0)
      # Misclassify based on the probability of HRCT- given CXR-
      if(runif(1) > rev_g1) { new_v = 1 }  # Misclassify to HRCT+
    }
    miss_v[i] = new_v  # 
  }
  
  return(miss_v)
}


# this expects the data columns and the betas to be in the same order, 
# following the intercept beta
predict <- function(data, beta_v) {
  n_rows = NULL
  if(is.vector(data)) {n_rows=length(data)}
  else{n_rows=dim(data)[1]}
  xb=apply(as.matrix(cbind(rep(1,n_rows),data)) %*% diag(beta_v), 1, sum)
  pred = exp(xb) / (1+exp(xb))
}

# Standard ROC Analysis
ROC <- function(pheno, score) {
  tp=NULL
  fp=NULL
  auc=0 
  n.pts=100
  cut = seq(max(score), min(score), length=n.pts)
  for(c in cut) { 
    tp = c(tp, sum(score>c & pheno==1, na.rm=T)/sum(pheno==1 & !is.na(score), na.rm=T))
    fp = c(fp, sum(score>c & pheno==0, na.rm=T)/sum(pheno==0 & !is.na(score), na.rm=T))
    if(c<max(score)) {auc = auc + ((tail(fp,n=1)-tail(fp,n=2)[1])*tail(tp,n=2)[1]) + .5*((tail(fp,n=1)-tail(fp,n=2)[1])*((tail(tp,n=1)-tail(tp,n=2)[1])))}
  }
  list(tp=tp, fp=fp, auc=auc)
}

# Misclassification-adjusted ROC procedure
mis_ROC <- function(y, score, gamma0, gamma1) {
  
  # Compute conidtional predictive probability
  Pr_T = ((gamma1-y*(2*gamma1-1)) * score) / ((1-y) + ((-1)^(1-y))*((1-gamma1-gamma0)*score+gamma0))
  
  tp=NULL
  fp=NULL
  auc=0 
  n.pts=100
  cut = seq(max(score), min(score), length=n.pts)
  t=cbind(Pr_T, 1-Pr_T, score)
  
  for(c in cut) { 
    tp = c(tp, sum(subset(t, score>c)[,1]) /sum(Pr_T))
    fp = c(fp, sum(subset(t, score>c)[,2]) /sum(1-Pr_T))
    if(c<max(score)) {auc = auc + ((tail(fp,n=1)-tail(fp,n=2)[1])*tail(tp,n=2)[1]) + .5*((tail(fp,n=1)-tail(fp,n=2)[1])*((tail(tp,n=1)-tail(tp,n=2)[1])))}
  }
  list(tp=tp, fp=fp, auc=auc, Pr_T=Pr_T)
}

### Return the probability from a logitistic model
logit.pred = function(beta,x) {
  xb= beta[1] + x %*% beta[-1]  
  return(exp(xb)/(1+exp(xb)))
}

#### END OF FUNCTIONS #######################
