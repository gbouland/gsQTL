#' Creates a gene set matrix from a gene expression matrix'
#' @param pathwayInfo A data.frame with two columns: pathname and symbol
#' where pathname represent the gene set/pathway and symbol represents the genes
#' belonging to that pathway
#' @param expression A matrix with features on the rows and individual on the columns
#' @param minGenes Minimal number of genes a gene set must consists of. This is
#' the minimal number after filtering the expression matrix based on features 
#' that overlap with pathwayInfo
#' @param minExplained minimal requirement of variance explained of PC1
#' @param regressOut Boolean stating whether variables need to be regressed out
#' @param samplesheet data.frame of meta data, columns corresponding to the rows 
#' in expression 
#' @param cov vector of strings corresponding to columns in the samplesheet that need
#' to be regressed out in the same order as expression. Rows in expression are
#' are columns in sample 
#' @param linear Boolean indicating whether linear or kernel PCA needs to be used
#' @param sigma inverse kernel width for the Radial Basis kernel
#' @export
#' @import kernlab
#' @importFrom stats na.omit cor prcomp
#' @examples \dontrun{
#' library(gsQTL)
#' gsQTLs
#' }
#' @return Matrix of gene set by ID
makePathwayExpression <- function(pathwayInfo, expression, minGenes, minExplained, regressOut = F, samplesheet = NULL, cov = NULL, linear = TRUE, sigma = 10){
  pathways <- unique(pathwayInfo$pathname)
  ID_path <- pathwayInfo[pathwayInfo$symbol %in% rownames(expression),]
  if(regressOut){
    genes_test <- unique(ID_path$symbol)
    resids <- lapply(genes_test,function(gene){
      getResiduals(expression,samplesheet,gene,cov)
    }) |> do.call(what = 'rbind')
    rownames(resids) <- genes_test
    colnames(resids) <- colnames(expression)
    expression <- resids
  }
  path_expr <- lapply(pathways,function(x){
    if(match(x,pathways)%%1000 ==0){
      message(sprintf("%s of %s",match(x,pathways),length(pathways)))
    }
    genes <- ID_path[ID_path$pathname == x, "symbol"]
    genes <- unique(genes)
    if(length(genes) <= minGenes){
      return(NA)
    }
    else{
      geneSums <- rowSums(expression[genes,])
      genes <- names(geneSums)[geneSums!=0]
      if(linear){
        PCs <- prcomp(scale(t(expression[genes,])))
        PC1 <- PCs$x[,1]
        rotation <- PCs$rotation[,1]
        eigenValues <- PCs$sdev
      }else{
        PCs <- kpca(t(expression[genes,]),
                    kernel = "rbfdot",
                    kpar = list(sigma = sigma),
                    features = 0)
        PC1 <- PCs@pcv[,1]
        rotation <- cor(PC1,t(expression[genes,]))#PCs@rotated[,1] 
        eigenValues <- PCs@eig
      }
      if(eigenValues[1] / sum(eigenValues) >= minExplained){
        sumLoad <- sum(rotation)
        if(sumLoad<0){
          out <- PC1 * -1
        }else{
          out <- PC1
        }
        return(out)
      }else{
        return(NA)
      }
    }
  })|> do.call(what = "rbind")
  rownames(path_expr) <- pathways
  path_expr <- na.omit(path_expr)
  colnames(path_expr) <- colnames(expression)
  return(path_expr)
}