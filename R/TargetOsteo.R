options('DCC' = '/data/CCRBioinfo/projects/TargetOsteosarcoma/OtherData/DCC')

#' Get TARGET clinical data (only 89 discovery samples) as a data frame
#' 
#' @param DCC The string location of the DCC directory
#' 
#' @export
#' @examples
#' clindf = getClinical()
#' head(clindf)
getClinical = function(DCC = options('DCC')) {
  tmp = read.xls(file.path(DCC,'clinical/TARGET_OS_Discovery_Toronto_COG_89_FCCs_ClinicalData_7_29_2015_harmonized.xlsx'))
  rownames(tmp) = tmp$TARGET.USI
  return(tmp)
}

#' Get TARGET miRNA data (only 89 discovery samples) as an ExpressionSet
#' 
#' This function loads the miRNA data and does sample matching with 
#' the clinical data (loaded using \code{\link{getClinical()}}).
#' 
#' @param DCC The string location of the DCC directory
#' 
#' @import Biobase
#' 
#' @export
#' 
#' @examples
#' mirna = getmiRNA()
#' mirna
getmiRNA = function(DCC = options('DCC')) {
  clindf = getClinical(DCC)
  tmp    = read.delim(file.path(DCC,'miRNA_pcr/L2/miRNA dCt values.txt'),check.names = FALSE)
  colnames(tmp)[1]='miRNA'
  rownames(tmp) = tmp[,1]
  colnames(tmp) = sub('-01A','',colnames(tmp))
  tmpm = as.matrix(tmp[,-1])
  rownames(tmpm) = tmp[,1]
  matchingsamples = intersect(colnames(tmpm),rownames(clindf))
  e = ExpressionSet(assayData = tmpm[,matchingsamples],
                    phenoData = AnnotatedDataFrame(clindf[matchingsamples,]))
  fData(e) = data.frame(miRNA = tmp[,1])
  featureNames(e) = rownames(tmpm)
  return(e)
}
