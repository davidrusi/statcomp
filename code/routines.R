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


# coefSummary: summarize the estimated coefficients of a lm/glm object
# Input
# - lmfit: fitted model of class 'lm' or 'glm'
# - level: confidence level for the confidence intervals
# - digits: number of digits to report for estimated coefficients (P-values always reported to 5 digits)
# Output: tibble with one row per coefficient, indicating its point estimate, confidence interval and P-value
coefSummary= function(lmfit, level=0.95, digits=3) {
  if (class(lmfit) != 'lm') stop('lmfit must be of class lm')
  b= round(coef(lmfit), digits)
  ci= round(confint(lmfit, level=level), digits)
  ci= paste('(',ci[,1],',',ci[,2],')',sep='')
  pval= round(summary(lmfit)$coef[,'Pr(>|t|)'],5)
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
