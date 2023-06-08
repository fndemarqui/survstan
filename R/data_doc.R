#' IRESSA Pan-Asia Study (IPASS) data set
#'
#' @name ipass
#' @docType data
#' @author Fabio N. Demarqui \email{fndemarqui@est.ufmg.br}
#' @keywords datasets
#' @description Reconstructed IPASS clinical trial data reported in Argyropoulos and Unruh (2015). Although reconstructed, this data set preserves all features exhibited in references with full access to the observations from this clinical trial. The data base is related to the period of March 2006 to April 2008. The main purpose of the study is to compare the drug gefitinib against carboplatin/paclitaxel doublet chemotherapy as first line treatment, in terms of progression free survival (in months), to be applied to selected non-small-cell lung cancer (NSCLC) patients.
#' @format A data frame with 1217 rows and 3 variables:
#' \itemize{
#'   \item time: progression free survival (in months)
#'   \item status: failure indicator (1 - failure; 0 - otherwise)
#'   \item arm: (1 - gefitinib; 0 - carboplatin/paclitaxel doublet chemotherapy)
#' }
#' @references  Argyropoulos, C. and Unruh, M. L. (2015). Analysis of time to event outcomes in randomized controlled trials by generalized additive models. PLOS One 10, 1-33.
#'
NULL


#' Gastric cancer data set
#'
#' @name gastric
#' @docType data
#' @author Fabio N. Demarqui \email{fndemarqui@est.ufmg.br}
#' @keywords datasets
#' @description Data set from a clinical trial conducted by the Gastrointestinal Tumor Study Group (GTSG) in 1982. The data set refers to the survival times of patients with locally nonresectable gastric cancer. Patients were either treated with chemotherapy combined with radiation or chemotherapy alone.
#' @format A data frame with 90 rows and 3 variables:
#' \itemize{
#'   \item time: survival times (in days)
#'   \item status: failure indicator (1 - failure; 0 - otherwise)
#'   \item trt: treatments (1 - chemotherapy + radiation; 0 - chemotherapy alone)
#' }
#' @references  Gastrointestinal Tumor Study Group. (1982) A Comparison of Combination Chemotherapy and Combined Modality Therapy for Locally Advanced Gastric Carcinoma. Cancer 49:1771-7.
#'
NULL
