#' @description Visualize the dFC time dependence  modeled in dynamic functional connectivity modelling along with its confidence band.
#' @title 95\% Coverage time dependence plots for dFCM
#' @param dataList.CI The path of confidence interval data
#' @param scan.length integer, usually set as number of total timepoints per scan
#' @param rois Numerical array, regions of interests to generate plots
#' @param save_fig logical, whether or not the plots should be saved as files
#' @param output_dir character, output directory of PNG-format figures
#' @param cores integer, indicate the number of cores to use in parallel execution
#' @examples 
#' data(DynModel_results) # Load sample model result object
#' scan.length <- DynModel_results$params$ntp_per_scan
#' plot_generate_timeDependence(DynModel_results$est_CI[1], scan.length, rois=c(114,134))
#' \dontshow{
#' rm(list = c('scan.length', 'DynModel_results'))
#' gc()
#' }
#' @export
#'

plot_generate_timeDependence <- function(dataList.CI, scan.length, rois, save_fig = FALSE, output_dir = NULL, cores = 0) {
    
    
    if (is.null(rois)) {
        stop(paste("Error: No ROIs input, nothing to generate."))
    }
    
    if(isTRUE(save_fig)){
      if (!is.null(output_dir)) {
        if (!dir.exists(output_dir)) 
          dir.create(output_dir, recursive = T, showWarnings = F)
      }
    }

    
    
    # Parse roi from comparison names
    comparisonNames <- names(dataList.CI)
    
    roi1 <- as.integer(unlist(lapply(strsplit(names(dataList.CI), "_"), function(x) return(x[2]))))
    roi2 <- as.integer(unlist(lapply(strsplit(names(dataList.CI), "_"), function(x) return(x[3]))))
    
    main_filenames <- comparisonNames
    
    names.roi.plots <- roi_names[rois, 6]
    
    Namesnames <- cbind(rois, names.roi.plots)
    
    num.of.plots <- length(main_filenames)
    
    parallel <- ifelse(cores < 1, FALSE, TRUE)
    
    if (parallel) {
        # Save system initial status - pre-setup for parallel computing
        print("=========== Loading Parallel Computing Packages ==========")
        "%dopar%" <- foreach::"%dopar%"
        doParallel::registerDoParallel(cores = cores)
        foreach::foreach(i = 1:num.of.plots, .packages = c("ggplot2", "doParallel", "dfConn"), .export = c("corr.scale", "plot_timeDependence_dFCM")) %dopar% {
            
            data.p <- corr.scale(x = dataList.CI[[i]])
            # print('Plotting.')
            g.list <- list()
            g.list[[i]] <- plot_timeDependence_dFCM(data.p = data.p, roi1 = roi1[i], roi2 = roi2[i], n_timepoints = scan.length)
            
            print(sprintf("finish plot %d", i))
            if (isTRUE(save_fig)) {
                ggplot2::ggsave(file.path(output_dir, comparisonNames[i]), plot = g.list[[i]], device = NULL, path = NULL, width = 13, height = 11, scale = 1, units = "in", 
                  dpi = 800, limitsize = TRUE)
            }
            
            
        }
        
    } else {
        g.list <- list()
        for (i in 1:num.of.plots) {
            data.p <- corr.scale(x = dataList.CI[[i]])
            
            # ggplot
            g <- plot_timeDependence_dFCM(data.p = data.p, roi1 = roi1[i], roi2 = roi2[i], n_timepoints = scan.length)
            g.list[[i]] <- g
            # Save ggplot in file
            if (isTRUE(save_fig)) {
                ggplot2::ggsave(file.path(output_dir, comparisonNames[i]), plot = g, device = NULL, path = NULL, width = 13, height = 11, scale = 1, units = "in", dpi = 800, 
                  limitsize = TRUE)
            }
            
            
        }
    }
    
    return(g.list)
    
}
