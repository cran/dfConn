#' Sample median output of Multivariate Linear Process Bootstrap (MLPB) 
#'
#' A dataset containing summarized by calculating the median for dynamic functional connectivity estimates of pairwise region of interest (ROI) comparisions.
#' 
#' @docType data
#' 
#' @usage data(MLPB_output_median)
#' 
#' @keywords datasets
#' 
#' @format A list with 3 matrices, each matrix has 5 rows and 731 columns, where rows represent subjects-level multivariate linear process bootstrapping (MLPB) estimates and columns represent median of MLPB estimates at each window.
#' \describe{
#' 
#'   \item{median_1_2}{matrix with 5 rows and 731 columns, containing median of multivariate linear process bootstrapping results between region 1 and 2 in subject-level.}
#'   \item{median_1_3}{matrix with 5 rows and 731 columns, containing median of multivariate linear process bootstrapping results between region 1 and 3 in subject-level.}
#'   \item{median_2_3}{matrix with 5 rows and 731 columns, containing median of multivariate linear process bootstrapping results between region 2 and 3 in subject-level.}
#' }
"MLPB_output_median"


#' Sample confidence interval output of dynamic functional connectivity (dFC) model 
#'
#' A sample confidence interval estimates (between pairwise regions of interest 1 and 2) coming from dynamic functional connectivity model
#' 
#' @docType data
#' 
#' @usage data(DynModel_results)
#' 
#' @keywords datasets
#' 
#' @format A list with 5 objects:
#' \describe{
#' \item{params}{Store the parameters for the linear mixed model}
#' \item{parallel}{Store the information about parallel computing environment and parameters}
#' \item{output_by_row}{Store the model results information for data between each pair of regions}
#' \item{modelDyn_results}{A combined dataframe with all the model results informations}
#' \item{est_CI}{A list of matrices, contains information of confidence band estimate, see description below.}
#' In each matrix of \code{est_CI}:
#' 
#'   \item{row 1}{Base line estimate, difference between condition 1 and condition 2}
#'   \item{row 2}{Dynamic functional connecvity estimates for condition 1}
#'   \item{row 3}{Dynamic functional connecvity estimates for condition 2}
#'   \item{row 4}{Lower bound of 95\% confidence interval of condition-difference estimates}
#'   \item{row 5}{Upper bound of 95\% confidence interval of condition-difference estimates}
#'   \item{row 6}{Lower bound of 95\% confidence interval of condition-1 estimates}
#'   \item{row 7}{Upper bound of 95\% confidence interval of condition-1 estimates}
#'   \item{row 8}{Lower bound of 95\% confidence interval of condition-2 estimates}
#'   \item{row 9}{Upper bound of 95\% confidence interval of condition-2 estimates}
#' }
"DynModel_results"


#' Sample fMRI time series data matrices
#'
#' Sample functional magnetic resonance imaging (fMRI) data contains a list of data matrices with consisting time series data of two subjects, with 6 scans for each subject in total.
#' 
#' @docType data
#' 
#' @usage data(fMRI_dataList)
#' 
#' @keywords datasets
#'
#' @format A list consisting of two matrices, each with 278 rows and 750 columns:
#' \describe{
#'   \item{subject1}{sample time series data matrix for subject 1.}
#'   \item{subject2}{sample time series data matrix for subject 2.}
#'   
#'   
#'   For each data matrix in the list:
#'   \item{row}{index of regions of interest .}
#'   \item{column}{fMRI time series data at each time point.}
#' }
"fMRI_dataList"

#' A shrinked sample subject functional magnetic resonance imaging (fMRI) time series data
#'
#' Contains a list of data matrices with consisting time series data of only one subject, and only 1 scan is preserved. This shrinked dataset is useful for user to test the programs and functions.
#' 
#' @docType data
#' 
#' @usage data(fMRI_dataList_shrinked)
#' 
#' @keywords datasets
#' 
#' @format A list which contains a matrix with 2 rows and 19 columns:
#' \describe{
#'   \item{row}{index of regions of interest (ROIs).}
#'   \item{column}{fMRI time series data at each time point.}
#' }
"fMRI_dataList_shrinked"

#' Sample dynamic functional connectivity estimates summary tables
#'
#' Contains a list of data frames with comparisons and corresponding significant non-zero estimates coverage.
#' 
#' @docType data
#' 
#' @usage data(coverage.tbl.list)
#' 
#' @keywords datasets
#' 
#' @format A list which contains 9 data frames:
#' \describe{
#'   \item{coverage.zero.cond.diff}{Percentage of estimates for condition difference that are within 95\% confidence interval with a center of zero.}
#'   \item{coverage.zero.cond.1}{Percentage of estimates for condition 1 that are within 95\% confidence interval with a center of zero.}
#'   \item{coverage.zero.cond.2}{Percentage of estimates for condition 2 that are below 95\% confidence interval with a center of zero.}
#'   \item{coverage.below.zero.cond.2}{Percentage of estimates for condition 2 that are below 95\% confidence interval with a center of zero.}
#'   \item{coverage.below.zero.cond.diff}{Percentage of estimates for condition difference that are below 95\% confidence interval with a center of zero.}
#'   \item{coverage.below.zero.cond.1}{Percentage of estimates for condition 1 that are below 95\% confidence interval with a center of zero.}
#'   \item{coverage.above.zero.cond.diff}{Percentage of estimates for condition difference that are above 95\% confidence interval with a center of zero.}
#'   \item{coverage.above.zero.cond.1}{Percentage of estimates for condition 1 that are above 95\% confidence interval with a center of zero.}
#'   \item{coverage.above.zero.cond.2}{Percentage of estimates for condition 1 that are above 95\% confidence interval with a center of zero.}
#' }
"coverage.tbl.list"






