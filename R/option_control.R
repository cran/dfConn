#' @description Options control for running multivariate linear Multivariate Linear Process Bootstrap (MLPB) method.
#' @title Options control - Multivariate linear process bootstrap
#' @param input_list, list of matrices, subject-specific time series data.
#' @param output_dir character string, output directory for bootstrap samples.
#' @param rois vector of integers, specify the list of regions of interests
#' @param boot_rep integer, bootstrapping repetition times
#' @param subset vector of character, which subject to run, specify the index in \code{dataList}, for example, c(1, 2, ...) 
#' @param window_size integer, window size for sliding window technique
#' @param timepoints integer, number of timepoints in total
#' @param number_of_intervals integer, number of intervals in sliding window technique
#' @param n_boot integer, number of bootstrap sample to be generated in MLPB
#' @param cores integer, number of cores to register for parallel execution, if set to 1 a sequential job will be run
#' @param save_file_suffix character, suffix of output files, treated as a labels.
#' @return A list of options 
#' 
#' @examples 
#' # Summarize the output of linear mixed effect model for dynamic functional connectivity
#' data(fMRI_dataList)
#' op <- options_boot(input_list = fMRI_dataList,
#'                    output_dir = NULL,
#'                    rois = c(1,2,3),
#'                    window_size = 10, 
#'                    save_file_suffix = "",
#'                    number_of_interval = 1,
#'                    boot_rep = 250,
#'                    n_boot = 1,
#'                    cores = 1)
#'                    
#' \dontshow{
#' rm(list = c("fMRI_dataList","op"))
#' gc()
#' }
#' @export

options_boot <- function(input_list,
                         output_dir = NULL,
                         rois = NULL,
                         subset = NULL,
                         timepoints = NULL, 
                         window_size = 20, 
                         number_of_intervals = 1,
                         boot_rep = 250,
                         n_boot = 1,
                         cores = 1,
                         save_file_suffix = ""){
  
  everything_ok <- TRUE
  if(length(input_list) == 0){
    stop("No data containing in input_list")
  }else{
    # Check each element in the input list
    each_timepoints <- unlist(lapply(input_list, ncol))
    
    if(length(unique(each_timepoints))!=1){
      cat("Number of timepoints is not the same across the data list, please check.\n")
      everything_ok <- FALSE
    }
    
    i <- 0
    for (mat in input_list){
      if (is.null(mat)){
        cat("The ", i, "th", " data matrix in input_list is NULL. Please check again. \n")
        everything_ok <- FALSE
      }
      i <- i+1
    }
    if (is.null(timepoints)) timepoints <- ncol(input_list[[1]])
  }
  
  # Check output_dir
  if(is.null(output_dir)){
    cat("output directory is set to NULL, will not save output.\n")
    save_output <- F
  }else{
    save_output <- T
  }
  
  
  # Check rois
  if (is.null(rois)){
    rois <- 1:min(unlist(lapply(input_list, nrow)))
  }
  
  if (any(is.na(rois))){
    cat("rois contains NA, please check again.\n")
    everything_ok <- FALSE
  }else{
    if(any(rois<0)){
      cat("rois contains negative value, please check again.\n")
      everything_ok <- FALSE
    }
    max_roi <- max(rois)
    
    if(max_roi > max(unlist(lapply(input_list, nrow)))){
      cat("rois contains value larger than total number of regions in input_list, please check again.\n")
      everything_ok <- FALSE
    }
  }
  
  # check subject subset
  if (!is.null(subset)){
    if (any(!(subset %in% names(input_list)))){
      cat("subset contains subject(s) that is/are not in the input_list, please check.\n")
      everything_ok <- FALSE
    }
  }
  
  # check save_file_suffix
  if (length(save_file_suffix)>230){
    cat("save_file_suffix is too long, please shorten the file name.\n")
    everything_ok <- FALSE
  }
  
  
  # check windows size
  if (window_size > max(unlist(lapply(input_list, ncol)))){
    cat("window_size is larger than number of timepoint, please check.\n ")
    everything_ok <- FALSE
  }
  
  # check other arguments
  if (!is.numeric(number_of_intervals)){
    cat("number_of_intervals should be number, please check.\n ")
    everything_ok <- FALSE
  }
  
  if (!is.numeric(boot_rep)){
    cat("boot_rep should be number, please check.\n ")
    everything_ok <- FALSE
  }
  
  if (!is.numeric(n_boot)){
    cat("n_boot should be number, please check.\n ")
    everything_ok <- FALSE
  }
  
  if (!is.numeric(cores)){
    cat("cores should be number, please check.\n ")
    everything_ok <- FALSE
  }
  
  if(cores>1) parallel_run <- TRUE
  else parallel_run <- FALSE
  
  if(isTRUE(everything_ok)){
    return(list(dataList = input_list, 
                output_dir = output_dir,
                rois = rois, 
                subset = subset, 
                save_file_suffix = save_file_suffix,
                save_output = save_output,
                window_size = window_size, 
                number_of_intervals = number_of_intervals,
                boot_rep = boot_rep,
                n_boot = n_boot,
                cores = cores,
                parallel_run = parallel_run,
                timepoints = timepoints))
  }else{
    return()
  }
  
}


################# Linear Mixed Model Control #####################

#' @description Options control for running linear mixed model.
#' @title Options control - Linear mixed model
#' @param effective_tp integer, effective scan time points.
#' @param output_dir string, directory for storing output files.
#' @param cores integer, number of cores to register while running parallel jobs.
#' @param seed integer, random seed.
#' @param subjects character vector, names of subjects.
#' @param num.scan integer, number of scan.
#' @param ntps.per.scan integer, number of timepoints per scan.
#' @param ngrid integer, number of grids.
#' @param ci_level numeric, level of confidence interval, with default 0.975.
#' @param numIntKnots integer, number of knots for O’Sullivan penalized splines, default is 40
#' @return A list of options 
#' @examples 
#' data("MLPB_output_median")
#' 
#' subjects <- c('subject1', 'subject2', 'subject3', 'subject4', 'subject5')
#' 
#' # In our demo data, each subject has a scan with a total of 750 time points
#' time.points <- c(1:105, 126:230, 251:355,
#'                  376:480, 501:605, 626:730) 
#'                  
#'                  
#' num.scan <- 6 # Each subject has 6 scans
#' ntps.per.scan <- 105 # Each scan has 105 time points
#' 
#' op <- options_lme(effective_tp = time.points, 
#'                  ntps.per.scan = ntps.per.scan,
#'                  subjects = subjects, 
#'                  num.scan = num.scan, 
#'                  cores = 5)
#' \dontshow{
#' rm(list = c("MLPB_output_median","op"))
#' gc()
#' }
#' @export

options_lme <- function(effective_tp, 
                        num.scan, 
                        ntps.per.scan,
                        output_dir = NULL,
                        subjects = NULL,
                        ci_level = 0.975,
                        numIntKnots = 40,
                        cores = 1,
                        ngrid = 201,
                        seed = 1114){
  
  
  everything_ok <- TRUE
  
  # Check output_dir
  if(is.null(output_dir)){
    cat("output directory is set to NULL, will not save output.\n")
    save_output <- F
  }else{
    save_output <- T
  }

  # check scan and effective time points
  if (ntps.per.scan * num.scan != length(effective_tp)){
    cat("Total number of effective time points should be equal to timepoint-per-scan * num.scan.\n ")
    everything_ok <- FALSE
  }
  

  if (!is.numeric(cores)){
    cat("cores should be number, please check.\n ")
    everything_ok <- FALSE
  }
  
  if(cores>1) parallel_run <- TRUE
  else parallel_run <- FALSE
  
  if(isTRUE(everything_ok)){
    return(list(output_dir = output_dir,
                effective_tp = effective_tp, 
                num.scan = num.scan, 
                ntps.per.scan = ntps.per.scan,
                subjects = subjects,
                ci_level = ci_level,
                cores = cores,
                ngrid = ngrid,
                numIntKnots = numIntKnots,
                seed = seed,
                parallel_run = parallel_run,
                save_output = save_output))
  }else{
    return()
  }
  
}
