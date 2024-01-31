#############################################################################################
##
## AUXILIARY ROUTINES FOR THE STATISTICAL COMPUTING COURSE
##
## AUTHOR: DAVID ROSSELL
##
#############################################################################################

# LIST OF FUNCTIONS
# - coefSummary: summarize the estimated coefficients of a lm/glm object

# PERMUTATION TESTS
# - lmResPerm.onevar: residual permutation test for a single variable in linear regression
# - lmResPerm: residual permutation test for all variables in linear regression

# AUXILIARY FUNCTIONS FOR LOGISTIC REGRESSION
# - flogregvec, flogreg: (negative) logistic regression log-likelihood
# - fplogreg: gradient and hessian of flogreg

# AUXILIARY FUNCTIONS FOR POISSON REGRESSION
# - fpoisregvec, fpoisreg: (negative) Poisson regression log-likelihood
# - fppoisreg: gradient and hessian of Poisson

# BOOTSTRAP FOR GLMs
# - bootGLM: Bootstrap confidence intervals for a GLM
# - mleGLM: obtain MLE under a GLM
# - bootCI: Extract confidence intervals from object returned by "boot"
# - bootstrapAnova: BOOTSTRAP PARAMETRIC TEST FOR NESTED RANDOM EFFECTS MODELS

# LOSS FUNCTIONS FOR CROSS-VALIDATION
# - cost_loglik_logistic: logistic regression log-likelihood loss
# - cost_misclass: proportion of miss-classified observations given predicted class probabilities








# coefSummary: summarize the estimated coefficients of a lm/glm object
# Input
# - lmfit: fitted model of class 'lm' or 'glm'
# - level: confidence level for the confidence intervals
# - digits: number of digits to report for estimated coefficients (P-values always reported to 5 digits)
# - transform: optional argument. If specified, it should be a function that will be applied to the estimated parameters. 
#   For example, set transform= function(x) exp(x) to return exp(coef(lmfit))
# Output: tibble with one row per coefficient, indicating its point estimate, confidence interval and P-value
coefSummary= function(lmfit, level=0.95, digits=3, transform) {
  require(tidyverse)
  if (!inherits(lmfit, "lm")) stop('lmfit must be of class lm')
  b= coef(lmfit)
  ci= confint(lmfit, level=level)
  if (!missing(transform)) {
    b= transform(b)
    ci= transform(ci)
  }
  b= round(b, digits); ci= round(ci, digits)
  ci= paste('(',ci[,1],',',ci[,2],')',sep='')
  pval= round(summary(lmfit)$coef[,4],5)
  pval[pval < 0.00001]= '<0.00001'
  ans= tibble(names(b), b, ci, pval)
  names(ans)= c('Parameter','Estimate','Conf. Int.','P-value')
  return(ans)
}



###########################################################################################
# PERMUTATION TESTS
###########################################################################################

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
    eperm= sample(e, size=length(e), replace=FALSE)
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
flogreg= function(beta, y, X) {
  ytX = matrix(y,nrow=1) %*% X
  Xbeta= as.vector(X %*% matrix(beta,ncol=1))
  ans= -sum(ytX * beta) + sum(log(1+exp(Xbeta)))
  return(ans)
}



#Gradient and Hessian of negative logistic regression log-likelihood
fplogreg= function(beta, y, X) {
  ytX = matrix(y,nrow=1) %*% X
  Xbeta= as.vector(X %*% beta)
  prob= expit(-Xbeta)
  g= -ytX + colSums(X * prob)
  H= t(X) %*% (prob*(1-prob) * X)
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


# BOOTSTRAP PARAMETRIC TEST FOR NESTED RANDOM EFFECTS MODELS
# Input
# - modelA: larger model, fitted by lmer
# - model0: null model, i.e. maller model nested within modelA, fitted by lmer
# - B: number of boostrap samples
# Output: likelihood ratio test, with Bootstrap-based P-value
# NOTE: Adapted from https://github.com/proback/BeyondMLR
bootstrapAnova = function(modelA, model0, B=1000){
  require(lme4)
  oneBootstrap = function(model0, modelA){ #LRT test statistic for 1 simulated dataset
    d = drop(simulate(model0))             #simulate data under model0
    m2 = lme4::refit(modelA, newresp=d)    #fit modelA to sim data
    m1 = lme4::refit(model0, newresp=d)    #fit model0 to sim data
    return(anova(m2,m1)$Chisq[2])          #LRT test statistic
  }  
  suppressMessages( nulldist <- replicate(B, oneBootstrap(model0, modelA)) )
  ret = anova(modelA, model0)
  ret$"Pr(>Chisq)"[2] = mean(nulldist > ret$Chisq[2])
  names(ret)[8] = "Pr_boot(>Chisq)"
  attr(ret, "heading") = c(attr(ret, "heading")[1], 
                            paste("Parametric bootstrap with", B,"samples."),
                            attr(ret, "heading")[-1])
  attr(ret, "nulldist") = nulldist
  return(ret)
}


###########################################################################################
## LOSS FUNCTIONS FOR CROSS-VALIDATION
###########################################################################################

# cost_loglik_logistic: logistic regression log-likelihood loss
# Input
# - yobs: observed class (must be binary)
# - ypred: predicted probability for yobs=1
# Output: logistic regression logistic regression log-likelihood
cost_loglik_logistic= function(yobs, ypred) {
  loglik= dbinom(yobs, size=1, prob=ypred, log=TRUE)
  return(sum(-loglik))
}

# cost_misclass: proportion of miss-classified observations given predicted class probabilities
# Input
# - yobs: observed class (must be binary)
# - ypred: predicted probability for yobs=1
# - threshold: if ypred > threshold then yobs=1 is predicted, else yobs=0 is predicted
# Output: proportion of miss-classified observations in yobs
cost_misclass= function(yobs, ypred, threshold=0.5) {
  err1= (ypred > threshold) & (yobs==0)
  err2= (ypred < threshold) & (yobs==1)
  ans= sum(err1 | err2) / length(yobs)
  return(ans)
}
