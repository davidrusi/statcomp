#############################################################################################
##
## AUXILIARY ROUTINES FOR THE STATISTICAL COMPUTING COURSE
##
## AUTHOR: DAVID ROSSELL
##
#############################################################################################

# LIST OF FUNCTIONS
# - coefSummary: summarize the estimated coefficients of a lm/glm object



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