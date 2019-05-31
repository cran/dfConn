#' @description Visualize non-zero coverage of dynamic functional connectivity modeling estimates.
#' @title Non-zero coverage correlation matrix visualization
#' @param coverage.tbl.list A list with a list of zero-coverage table, can be obtained from \code{data_summary_coverage}.
#' @param roi integer vector, specify the list of regions of interests.
#' @param save_fig logical, whether to save the figures locally.
#' @param output_dir character string, output directory of figures to save.
#' @param model_results data frame, output of either \code{lmmConn} or \code{lmmDyn}.
#' @param ntps.per.scan integer, number of timepoints per scan.
#' @param alpha alpha level of proportion-testing, default is 0.05.
#' @param thres.cond1 numeric, threshold for cond1 flavor proportion-test, default is 0.5.
#' @param thres.cond2 numeric, threshold for cond2 flavor proportion-test, default is 0.5.
#' @param thres.diff numeric, threshold for proportion-test of difference between condition 1 and condition 2 .
#' @examples 
#' # load sample summary of the linear mixed model output
#' data(coverage.tbl.list)
#' data(DynModel_results)
#' nzc_vis(coverage.tbl.list, DynModel_results$modelDyn_results,
#'         roi=c(54,191,235), DynModel_results$params$ntp_per_scan, save_fig=FALSE)
#' \dontshow{
#' rm(list=c('DynModel_results','coverage.tbl.list'))
#' gc()}       
#' @export
nzc_vis <- function(coverage.tbl.list, model_results, roi, ntps.per.scan, save_fig = FALSE, output_dir = NULL, alpha = 0.05, thres.cond1 = 0.5, thres.cond2 = 0.5, thres.diff = 0.1) {
    

  # Precheck
  if (isTRUE(save_fig)) {
    if (!is.null(output_dir)) {
            dir.create(output_dir, recursive = T, showWarnings = F)
    }
  }
    
    roi <- unique(roi)  # Remove duplicate roi
    
    mat_plot_object <- list()
    mat_plot_object[["ntp_per_scan"]] <- ntps.per.scan
    mat_plot_object[["save"]] <- save_fig
    mat_plot_object[["roi"]] <- roi
    mat_plot_object[["roi_names2plot"]] <- roi_names[roi, ]$Complete_names
    mat_plot_object[["thres.cond1"]] <- thres.cond1
    mat_plot_object[["thres.cond2"]] <- thres.cond2
    mat_plot_object[["thres.diff"]] <- thres.diff
    
    
    ################################################################# 
    lp.roi <- length(roi)
    
    
    mat_plot_object[["model_results"]] <- model_results
    mat_plot_object[["tbl.list"]] <- coverage.tbl.list
    mat_plot_object[["alpha"]] <- alpha
    
    
    cat("========================================================\n")
    cat("       Mat object attributes successfuly loaded.\n")
    
    
    
    cat("========================================================\n")
    cat("                   Generating Matrices.\n")
    
    ## Obtain value matrix for each direction of lag
    mat_list <- matrix_construction(mat_plot_object)
    
    cat("======================================================== \n")
    cat("          Successfully Generating Matrices. \n")
    
    value.map.roi.diff <- mat_list$diff
    value.map.roi.cond1 <- mat_list$cond1
    value.map.roi.cond2 <- mat_list$cond2
    
    
    cat("======================================================== \n")
    cat("                 Generating Figures \n")
    
    # print(value.map.roi.diff.neg)
    
    
    diag(value.map.roi.diff) <- -0.01
    
    tab1 <- rep(1, lp.roi)
    tab2 <- cumsum(tab1)
    tab2 <- c(0, (tab2))
    tab2.tick <- seq(0, lp.roi, 0.5)
    
    rois_yeo_names <- matplot_roisnames_setup(mat_plot_object[["roi_names2plot"]])
    
    label.x <- rois_yeo_names
    label.y <- rois_yeo_names
    
    
    
    ####################################### cond1 ##############################
    
    together.cond1.diff <- matrix(NA, lp.roi, lp.roi)
    
    for (i in c(1:(lp.roi - 1))) {
        for (j in c((i + 1):(lp.roi))) {
            
            together.cond1.diff[i, j] <- as.numeric(value.map.roi.cond1[i, j])
            together.cond1.diff[j, i] <- as.numeric(value.map.roi.diff[j, i])
            
        }
    }
    
    diag(together.cond1.diff) <- -0.01
    plot_non_zero_coverage(mat = together.cond1.diff, save_fig = save_fig, file = file.path(output_dir, "nzc_cond1_diff_combined_sig_diff.png", sep = ""), 
                           label.x = label.x, label.y = label.y, 
                           title = "Non-zero coverage for condition 1 and condition difference", 
                           height = 17, width = 17, type = "combined", res = 700, mar = c(14, 14, 3.5, 6.5))
    
    
    cat("Condition 1 vs difference matrix plot generated. \n")
    ############################################### cond2 ##############################
    
    together.cond1.cond2 <- matrix(NA, lp.roi, lp.roi)
    
    for (i in c(1:(lp.roi - 1))) {
        for (j in c((i + 1):(lp.roi))) {
            
            together.cond1.cond2[i, j] <- as.numeric(value.map.roi.cond1[i, j])
            together.cond1.cond2[j, i] <- as.numeric(value.map.roi.cond2[j, i])
            
        }
    }
    
    diag(together.cond1.cond2) <- -0.01
    
    tab1 <- rep(1, lp.roi)
    tab2 <- cumsum(tab1)
    tab2 <- c(0, (tab2))
    tab2.tick <- seq(0, lp.roi, 0.5)
    
    ############################################################### 
    
    diag(together.cond1.cond2) <- -0.01
    
    title.cond1.cond2 <- sprintf("Non-zero coverage for condition1 and condition2 (p=%f)", thres.cond1)
    
    plot_non_zero_coverage(mat = together.cond1.cond2, save_fig = save_fig,
                           file = paste(output_dir, "nzc_cond1_cond2_sig_diff.png", sep = ""),
                           label.x = label.x, label.y = label.y, title = title.cond1.cond2, 
                           height = 17, width = 17, type = "combined", res = 700, mar = c(14, 14, 3.5, 6.5))
    diag(together.cond1.cond2) <- 0
    cat("Condition 1 vs Condition 2 matrix plot generated. \n")
}
