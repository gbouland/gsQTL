#' Regresses out covariates
#' @param expr A matrix with features on the rows and individual on the columns
#' @param samplesheet data.frame of meta data, columns corresponding to the rows 
#' in expression
#' @param name names of the features
#' @param cov vector of strings corresponding to columns in the samplesheet that need
#' to be regressed out in the same order as expression. Rows in expression are
#' are columns in sample 
#' @export
#' @importFrom stats lm
#' @examples \dontrun{
#' library(gsQTL)
#' }
getResiduals <- function(expr,samplesheet,name,cov){
  tmpData <- samplesheet
  tmpData[,"toResid"] <- expr[name,]
  form <- sprintf("%s ~ %s","toResid",paste0(cov, collapse = " + "))
  resid <- lm(formula = form, data = tmpData)$residuals
  return(unname(resid))
}