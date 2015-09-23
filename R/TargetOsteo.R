options('DCC' = '/data/CCRBioinfo/projects/TargetOsteosarcoma/OtherData/DCC')

#' Get TARGET clinical data (only 89 discovery samples) as a data frame. The
#' actual file used is "TARGET_OS_Discovery_Toronto_COG_89_FCCs_ClinicalData_7_29_2015_harmonized.xlsx".
#'
#' @param DCC The string location of the DCC directory
#'
#' @importFrom gdata read.xls
#'
#' @export
#' @examples
#' clindf = getClinical()
#' head(clindf)
getClinical = function(DCC = options('DCC')) {
  tmp = gdata::read.xls(file.path(DCC,'clinical/TARGET_OS_Discovery_Toronto_COG_89_FCCs_ClinicalData_7_29_2015_harmonized.xlsx'))
  rownames(tmp) = tmp$TARGET.USI
  return(tmp)
}

#' Get TARGET miRNA data (only discovery samples) as an ExpressionSet
#'
#' This function loads the miRNA data and does sample matching with
#' the clinical data (loaded using \code{\link{getClinical()}}). The
#' actual file name used is "miRNA dCt values.txt". No further normalization
#' is performed.
#'
#' @param DCC The string location of the DCC directory
#'
#' @import Biobase
#' @importFrom gdata read.xls
#'
#' @export
#'
#' @examples
#' mirna = getmiRNA()
#' mirna
getmiRNA = function(DCC = options('DCC'), hsaOnly=TRUE) {
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
  if(hsaOnly) {
    e = e[grepl('hsa',fData(e)[,1]),]
  }
  return(e)
}

#' Get TARGET methylation data (only discovery samples) as SummarizedExperiment
#'
#' This is a specialized function to read TARGET methylation data from the
#' DCC repository structure and join it up with clinical data.  Use getAnnotation(obj)
#' to get associated feature annotation.
#'
#' @param DCC The string location of the DCC directory
#'
#' @value A \link{\code{GenomicMethylSet-class}} (subclass of SummarizedExperiment)
#'
#' @import minfi
#' @export
#'
#' @example
#' \donotrun {
#' methdat = getMethylation()
#' }
getMethylation = function(DCC=options('DCC'),normalization=preprocessIllumina) {
  clin      = getClinical()
  methpath  = file.path(DCC, 'methylation_array')
  idatpath  = file.path(methpath, 'L1')
  sdrfpath  = file.path(methpath, 'METADATA')
  sdrf      = read.delim(file.path(sdrfpath, 'MAGE-TAB_TARGET_OS_Meth_Illumina_20140327.sdrf.txt'))
  sdrf      = sdrf[sdrf$Label=='Cy3', ]
  rownames(sdrf) = make.unique(as.character(sdrf[,1]))
  olapsamps = intersect(clin$TARGET.USI, sdrf$Source.Name)
  sdrf      = sdrf[olapsamps, ]
  clin      = clin[olapsamps, ]
  basenames = sdrf$Array.Data.File
  names(basenames) = sdrf$Source.Name
  dat = normalization(read.450k(file.path(idatpath,basenames)))
  gset = mapToGenome(dat)
  colData(gset) = DataFrame(clin)
  return(gset)
}


