normTmmSE <- function(SE, useNormFactors=TRUE, priorCount=0.25){
  
  colDataSlots <- c('norm.factors', 'lib.size')
  assaySlots <- c('cpm', 'logcpm')
  
  if(any(colDataSlots %in%  names(colData(SE)))){
    msg <- colDataSlots[which(colDataSlots %in% names(colData(SE)))]
    warning(paste(c("The following colData columns will be recomputed:", msg), collapse=" "))
  }
  
  if(any(assaySlots %in%  names(assays(SE)))){
    msg <- assaySlots[which(assaySlots %in% names(assays(SE)))]
    warning(paste(c("The following SE assays will be recomputed:", msg), collapse=" "))
  }
  
  norm <- edgeR::calcNormFactors(SE, method='TMM')  #Returns DGElist object
  
  SE$norm.factors <- norm$samples$norm.factors
  SE$lib.size <- norm$samples$lib.size
  
  # Compute cpm/logcpm on DGElist so it can consider norm.factors and add to assays
  assays(SE)$cpm <- edgeR::cpm(norm, normalized.lib.sizes=useNormFactors)
  assays(SE)$logcpm <- edgeR::cpm(norm, log=TRUE, prior.count=priorCount, normalized.lib.sizes=useNormFactors)
  
  return(SE)
}