#' Runs linear model
#' @param expression A matrix with features on the rows and individual on the columns
#' @param samplesheet data.frame of meta data, columns corresponding to the rows 
#' in expression
#' @param predictor predictor variable for the linear model. Should be a column 
#' in samplesheet
#' @param cov vector of strings corresponding to columns in the samplesheet that need
#' to be regressed out in the same order as expression. Rows in expression are
#' are columns in sample
#' @param padj type of multiple testing correction. See p.adjust
#' @param type string linear or LR for logistic regression
#' @export
#' @importFrom stats lm glm p.adjust
#' @examples \dontrun{
#' library(gsQTL)
#' }
#' @return Matrix of gene set by ID
runLinearModel <- function(expression, samplesheet, predictor, cov = NULL, padj, type = "linear"){
  outcomes <- rownames(expression)
  rownames(expression) <- paste0(make.names(rownames(expression)),"_expr")
  tmpData <- cbind(samplesheet,t(expression))
  tmpData <- tmpData[!is.na(tmpData[,predictor]),]
  res <- lapply(rownames(expression),function(x){
    #if(match(x,rownames(expression))%%500 ==0){message(match(x,rownames(expression)))}
    
    if(!is.null((cov))){
      form <- sprintf("%s ~ %s + %s",x,predictor, paste0(cov, collapse = " + "))
    }
    else{
      form <- sprintf("%s ~ %s",x,predictor)
    }
    
    if(type == "LR"){
      form <- sprintf("%s ~ %s + %s",predictor,x, paste0(cov, collapse = " + "))
      res <- glm(formula = form, data = tmpData,family = "binomial")
      summ <- summary(res)
    }else{
      res <- lm(formula = form,data = tmpData)
      summ <- summary(res)
    }
    out <- data.frame(summ$coefficients[2,1],
                      summ$coefficients[2,2],
                      summ$coefficients[2,3],
                      summ$coefficients[2,4])
    
    colnames(out) <- switch(type,
                            "LR" = c("Estimate","stdError","z","p"),
                            "linear" = c("Estimate","stdError","t","p"))
    
    
    return(out)
  }) |> do.call(what = "rbind")
  
  res$name <- outcomes
  res$p_adj <- p.adjust(res$p,method = padj)
  res <- res[order(res$p),]
  return(res)
}