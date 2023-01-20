#############################################################################################
##
## AUXILIARY ROUTINES FOR THE STATISTICAL COMPUTING COURSE
##
## AUTHOR: DAVID ROSSELL
##
#############################################################################################

# LIST OF FUNCTIONS
# - coefSummary: summarize the estimated coefficients of a lm/glm object
# - lmResPerm.onevar: residual permutation test for a single variable in linear regression
# - lmResPerm: residual permutation test for all variables in linear regression
# - flogregvec, flogreg: (negative) logistic regression log-likelihood
# - fplogreg: gradient and hessian of flogreg
# - fpoisregvec, fpoisreg: (negative) Poisson regression log-likelihood
# - fppoisreg: gradient and hessian of Poisson




# coefSummary: summarize the estimated coefficients of a lm/glm object
# Input
# - lmfit: fitted model of class 'lm' or 'glm'
# - level: confidence level for the confidence intervals
# - digits: number of digits to report for estimated coefficients (P-values always reported to 5 digits)
# Output: tibble with one row per coefficient, indicating its point estimate, confidence interval and P-value
coefSummary= function(lmfit, level=0.95, digits=3) {
  require(tidyverse)
  if (!('lm' %in% class(lmfit))) stop('lmfit must be of class lm')
  b= round(coef(lmfit), digits)
  ci= round(confint(lmfit, level=level), digits)
  ci= paste('(',ci[,1],',',ci[,2],')',sep='')
  pval= round(summary(lmfit)$coef[,4],5)
  pval[pval < 0.00001]= '<0.00001'
  ans= tibble(names(b), b, ci, pval)
  names(ans)= c('Parameter','Estimate','Conf. Int.','P-value')
  return(ans)
}



#Residual permutation test for a single variable in linear regression
#Input
# - y: outcome
# - x1: covariate of interest
# - x2: other covariates (excluding the intercept)
# - B: number of residual permutations to consider
# Output: estimated coefficients for x1 in the B permutations
lmResPerm.onevar= function(y, x1, x2, B=5000) {
  if (!is.matrix(x2)) x2= as.matrix(x2)
  fit= lm(x1 ~ x2)
  x1hat= predict(fit)
  e= residuals(fit)
  bperm= double(B)
  for (b in 1:B) {
    eperm= e[sample(1:length(e), size=length(e), replace=FALSE)]
    x1tilde= x1hat + eperm
    fitperm= lm(y ~ x1tilde + x2)
    bperm[b]= coef(fitperm)['x1tilde']    
  }
  return(bperm)  
}


#Residual permutation test for all variables in linear regression
#Input
# - y: outcome
# - x: covariates (excluding the intercept)
# - B: number of residual permutations
#Output: a list with two elements
# - bperm: matrix with B rows and ncol(x) columns. Column j contains the permuted least-squares estimates for covariate j
# - pvalue: a vector with ncol(x) entries with the permutation-based P-value for each covariate
lmResPerm= function(y, x, B=5000) {
  if (!is.matrix(x)) x= as.matrix(x)
  p= ncol(x)
  bobs= coef(lm(y ~ x))[-1]  #exclude the intercept
  bperm= matrix(NA, nrow=B, ncol=p)
  pvalue= double(p)
  colnames(bperm)= names(pvalue)= colnames(x)
  for (j in 1:p) {
    bperm[,j]= lmResPerm.onevar(y, x1=x[,j], x2=x[,-j,drop=FALSE], B=B)
    pvalue[j]= mean(abs(bperm[,j]) > abs(bobs[j]))
  }
  ans= list(bperm=bperm, pvalue=pvalue)  
  return(ans)
}




###########################################################################################
# AUXILIARY FUNCTIONS FOR LOGISTIC REGRESSION
###########################################################################################

logit= function(z) log(z/(1-z))  #logit function
expit= function(z) 1/(1+exp(-z)) #inverse of logit function

#Negative logistic regression log-likelihood, for multiple beta values
flogregvec= function(beta, y, ytX, n, X=X, logscale=TRUE) {
  if (missing(ytX)) ytX = matrix(y,nrow=1) %*% X
  if (is.vector(beta)) beta= matrix(beta,ncol=1)
  apply(beta, 1, function(z) flogreg(z, ytX=ytX, n=n, X=X, logscale=logscale))
}

#Negative logistic regression log-likelihood
flogreg= function(beta, y, X, ytX, logscale=TRUE) {
  if (missing(ytX)) ytX = matrix(y,nrow=1) %*% X
  if (any(beta != 0)) {
    Xbeta= as.vector(X %*% matrix(beta,ncol=1))
    ans= -sum(ytX * beta) + sum(log(1+exp(Xbeta)))
  } else {
    n= length(y)
    ans= n * log(2)
  }
  if (!logscale) ans= exp(ans)
  return(ans)
}


#Gradient and Hessian of negative logistic regression log-likelihood
# If beta==0, only ytX, colSumX, XtX are used
# If beta!=0, only Xbeta and X are used
fplogreg= function(beta, y, ytX, Xbeta, X, colsumX, XtX) {
  if (missing(ytX)) ytX = matrix(y,nrow=1) %*% X
  if (any(beta != 0)) {
    if (missing(Xbeta)) Xbeta= as.vector(X %*% beta)
    prob= 1.0/(1.0+exp(-Xbeta))
    g= -ytX + colSums(X * prob)
    H= t(X) %*% (prob*(1-prob) * X)
  } else {
    g= -ytX + 0.5 * colsumX
    H= 0.25 * XtX
  }
  return(list(g=g,H=H))
}



###########################################################################################
# AUXILIARY FUNCTIONS FOR POISSON REGRESSION
###########################################################################################


#Negative logistic regression log-likelihood, for multiple beta values
fpoisregvec= function(beta, ytX, n, X=X, sumlfactorialy, logscale=TRUE) {
  if (is.vector(beta)) beta= matrix(beta,ncol=1)
  apply(beta, 1, function(z) fpoisreg(z, ytX=ytX, n=n, X=X, sumlfactorialy=sumlfactorialy, logscale=logscale))
}

#Negative logistic regression log-likelihood
fpoisreg= function(beta, ytX, Xbeta, n, X, sumlfactorialy, logscale=TRUE) {
  if (any(beta != 0)) {
    if (missing(Xbeta)) Xbeta= as.vector(X %*% matrix(beta,ncol=1))
    ans= -sum(ytX * beta) + sum(exp(Xbeta)) + sumlfactorialy
  } else {
    if (missing(n)) stop("If beta==0, n must be specified")
    ans= n + sumlfactorialy
  }
  if (!logscale) ans= exp(ans)
  return(ans)
}

#Gradient and Hessian of negative logistic regression log-likelihood
# If beta==0, only ytX, colSumX, XtX are used
# If beta!=0, only Xbeta and X are used
fppoisreg= function(beta, ytX, Xbeta, X, colsumX, XtX) {
  if (any(beta != 0)) {
    if (missing(Xbeta)) Xbeta= as.vector(X %*% beta)
    mu= exp(Xbeta)
    g= -ytX + colSums(X * mu)
    H= t(X) %*% (mu * X)
  } else {
    g= -ytX + colsumX
    H= XtX
  }
  return(list(g=g,H=H))
}


###########################################################################################
## BOOTSTRAP FOR GLMs
###########################################################################################

#Bootstrap confidence intervals for a GLM
# Input
# - data: dataset
# - formula: GLM formula of the type formula(y ~ x)
# - family: GLM family
# - level: confidence level
# - R: number of bootstrap samples
# Output: a tibble with a confidence interval for each regression parameter
bootGLM= function(data, formula, family, level=0.95, R=2000) {
  require(tidyverse)
  require(boot)
  fit= glm(formula, data=data, family=family)
  bootfit= boot(data, mleGLM, R=R, formula=formula, family=family)
  ans= bootCI(fit, bootfit, level=level)
  return(ans)
}

#Obtain MLE under a GLM
# Input
# - data: dataset
# - indices: MLE is obtained for data[indices,]
# - formula: GLM formula of the type formula(y ~ x)
# - family: GLM family
# Output: MLE
mleGLM= function(data, indices=1:nrow(data), formula, family) {
  if (missing(family)) stop("family must be specified")
  fit= glm(formula, data=data[indices,], family=family)
  return(coef(fit))
}

#Extract confidence intervals from object returned by "boot"
# Input
# - fit: object returned by glm, containing the MLE for the original observed data
# - bootfit: object returned by "boot"
# - level: confidence level
# Output: a tibble with a confidence interval for each regression parameter
bootCI= function(fit, bootfit, level=0.95) {
  bhat= coef(fit)
  probs= c((1-level)/2, 1 - (1-level)/2)
  bhat.boot= map_df(as_tibble(bootfit$t), quantile, probs=probs)
  bhat.boot= cbind(bhat, bhat.boot)
  return(bhat.boot)
}

