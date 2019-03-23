#' @description Summarize the data from dynamic functional connectivity linear mixed modeling.
#' @title Obtain the non-zero coverage for both two conditions and the condition-based difference
#' @param dataList.CI list, a list of data matrices which include the confidence interval estimates obtained from linear mixed models.
#' @param save_res logical, whether to save the result or not, if true, an output directory should be provided
#' @param output_dir character, directory for output files
#' @return A list which contains 9 data frames: 
#'  \describe{
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
#' 
#' @examples 
#' # Summarize the output of linear mixed effect model for dynamic functional connectivity
#' data(DynModel_results)
#' data_summary_coverage(DynModel_results$est_CI)
#' \dontshow{
#' rm(list = c("DynModel_results"))
#' gc()
#' }
#' @family summary
#' @export

data_summary_coverage <- function(dataList.CI, save_res = FALSE, output_dir = NULL) {
    

    
    names.CI.files <- dataList.CI
    if (isTRUE(save_res)){
      if (!is.null(output_dir)){
        output_dir <- paste(output_dir, "/coverage_above_below_zero/", sep = "")
        dir.create(output_dir, showWarnings = F, recursive = T)
      }else{
        stop("Output directory being set as NULL")
      }
    }

    comparisonsName <- names(dataList.CI)
    
    coverage <- table_summarize.coverage(dataList.CI, comparisonsName)
    
    
    ## ======================= Write csv file into output directory ===================================
    if (isTRUE(save_res)) {
        utils::write.csv(x = coverage$coverage.zero.cond.diff, file = paste(output_dir, "coverage_zero_cond_diff.csv", sep = ""), row.names = F)
        
        utils::write.csv(x = coverage$coverage.zero.cond.1, file = paste(output_dir, "coverage_zero_cond_1.csv", sep = ""), row.names = F)
        
        utils::write.csv(x = coverage$coverage.zero.cond.2, file = paste(output_dir, "coverage_zero_cond_2.csv", sep = ""), row.names = F)
        
        utils::write.csv(x = coverage$coverage.above.zero.cond.2, file = paste(output_dir, "cov_above_zero_cond_2.csv", sep = ""), row.names = F)
        
        utils::write.csv(x = coverage$coverage.above.zero.cond.1, file = paste(output_dir, "cov_above_zero_cond_1.csv", sep = ""), row.names = F)
        
        utils::write.csv(x = coverage$coverage.above.zero.cond.diff, file = paste(output_dir, "cov_above_zero_cond_diff.csv", sep = ""), row.names = F)
        
        utils::write.csv(x = coverage$coverage.below.zero.cond.2, file = paste(output_dir, "cov_below_zero_cond_2.csv", sep = ""), row.names = F)
        
        utils::write.csv(x = coverage$coverage.below.zero.cond.1, file = paste(output_dir, "cov_below_zero_cond_1.csv", sep = ""), row.names = F)
        
        utils::write.csv(x = coverage$coverage.below.zero.cond.diff, file = paste(output_dir, "cov_below_zero_cond_diff.csv", sep = ""), row.names = F)
    }
    
    return(coverage)
    
}
