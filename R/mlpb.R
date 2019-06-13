#' @description Multivariate Linear Process Bootstrap (MLPB) method to assess the uncertainty in dynamic functional connectivity (dFC) by providing its confidence band.
#' @title  Multivariate Linear Process Bootstrap Method
#' @param boot_options a list of option being specify by \code{option_control} function. 
#' @importFrom stats cor
#' @importFrom stats mad
#' @importFrom stats median
#' @importFrom parallel makePSOCKcluster
#' @importFrom gtools combinations
#' @importFrom Rcpp evalCpp
#' @importFrom utils combn
#' @import Rcpp
#' @details The \code{dataList} parameter is a list of matrices which contains time series data of each region of interest (ROI). Output directory is required here because the results to be generate is massive. 
#' @references Kudela et al. (2017) NeuroImage 10.1016/j.neuroimage.2017.01.056
#' (\href{https://www.ncbi.nlm.nih.gov/pubmed/28132931}{PubMed})
#' @examples
#' \dontshow{
#' # Load sample data
#' data(fMRI_dataList_shrinked)
#' boot_option <- options_boot(input_list = fMRI_dataList, output_dir = tempdir(),
#'  rois = c(1,2), timepoints = 5, window_size = 5)
#' MLPB_boot(boot_option)
#' rm(list=c('fMRI_dataList_shrinked','boot_option'))
#' gc()
#' }
#' \donttest{
#' 
#' # Load sample data
#' 
#' data(fMRI_dataList)
#' boot_option <- options_boot(input_list = fMRI_dataList, output_dir = tempdir(),
#'  rois = c(1,2,3), timepoints = 5, window_size = 5)
#' MLPB_boot(boot_option)
#'
#' }
#' @export
#' @keywords bootstarpping
#' 

MLPB_boot <- function(boot_options) {
    
  
    # Read parameters
    output_dir <- boot_options[["output_dir"]]
    dataList <- boot_options[["dataList"]]
    subset.subject <- boot_options[["subset"]]
    save_file_suffix <- boot_options[["save_file_suffix"]]
    window_size <- boot_options[["window_size"]]
    number_of_intervals <- boot_options[["number_of_intervals"]] 
    boot_rep <- boot_options[["boot_rep"]] 
    n_boot <- boot_options[["n_boot"]] 
    cores <- boot_options[["cores"]] 
    parallel_run <- boot_options[["parallel_run"]] 
    timepoints <- boot_options[["timepoints"]] 
    rois <- boot_options[["rois"]]
    
    if (is.null(output_dir)){
      stop("Output directory is not set, exiting...")
    }
    if (!dir.exists(output_dir)) 
        dir.create(output_dir, recursive = T)
    ########################################################### 
    
    ################################################## 
    
    
    ## To get the names of all the files
    nsubjects <- length(dataList)
    data.names <- vector()
    rozmiar <- matrix(NA, nsubjects, 2)
    names.data.sub <- matrix(NA, nsubjects, 1)
    tjj.name <- vector()
    
    if (is.null(subset.subject)) 
        subset.subject <- 1:nsubjects
    else
        subset.subject <- which(subset.subject %in% names(dataList))
    
    for (jj in subset.subject) {
        regions <- rois
        tjj <- dataList[[jj]]
        tjj <- data.matrix(tjj)
        tjj <- tjj[regions, ]
        throw.NAN1 <- list()
        
        if (is.null(dim(tjj))) {
            tjj <- t(as.matrix(tjj, 1, timepoints))
        }
        
        for (ipk in 1:timepoints) {
            throw.NAN1[[ipk]] <- which(!is.finite(tjj[, ipk]))
        }
        
        
        throw.NAN <- if (length(unique(unlist(throw.NAN1))) != 0) {
            unique(unlist(throw.NAN1))
        } else {
            0
        }
        
        
        if (length(throw.NAN) != 1) {
            tjj1 <- tjj[-throw.NAN, ]
        } else {
            if (throw.NAN != 0) {
                tjj1 <- tjj[-throw.NAN, ]
            } else {
                tjj1 <- tjj
            }
        }
        
        
        tjj.name[jj] <- names(dataList)[jj]
        dim.data <- dim(tjj)
        all.data <- matrix(NA, dim.data[1], dim.data[2])
        all.data <- tjj
        rozmiar[jj, ] <- dim(tjj1)
        names.data.sub[jj, 1] <- tjj.name[jj]
      
        path.out <- file.path(output_dir, paste(tjj.name[jj], save_file_suffix, "/", sep = ""))
        
        
        
        
        ############################################ 
        comparison <- t(utils::combn(rois[!(rois %in% throw.NAN)], 2)) # pairwise time series to be analysis

        
        ############# change window size 20, mad instead of sd
        n.sig <- timepoints  #number of time points
        up_limit <- n.sig - window_size + 1

        boot.cjj <- matrix(NA, nrow = boot_rep, ncol = (up_limit))
        
        
        coverage.zero <- vector()
        coverage.static <- vector()
        est.sd <- vector()
        sd.boot <- matrix()
        
        
        cat(paste("Bootstrapping on ", tjj.name[jj], "\n"))
        
        if(isTRUE(parallel_run)){
          cl <- parallel::makePSOCKcluster(2)
          doParallel::registerDoParallel(cl, cores = cores)
        }

        ############################### parallel
        
        '%dopar%' <- foreach::'%dopar%'
        
        i = 0
        foreach::foreach(i = 1:nrow(comparison)) %dopar% # for (i in 1:nrow(comparison))
        {
            
            boot.cjj <- matrix(NA, nrow = boot_rep, ncol = (up_limit))
            coverage.zero <- vector()
            coverage.static <- vector()
            est.sd <- vector()
            sd.boot <- matrix()
            
            ik <- comparison[i, 1]
            ij <- comparison[i, 2]
            
            kk <- ik
            ll <- ij
            name.of.subfile <- paste(tjj.name[jj], "_", as.character(regions[kk]), "_", as.character(regions[ll]), save_file_suffix, 
                "_medians.RData", sep = "")
            cat(sprintf("ROI %d and ROI %d\n", regions[kk], regions[ll]))
            
            xjj <- all.data[kk, ]
            yjj <- all.data[ll, ]
            
            Xjj <- rbind(xjj, yjj)
            
            
            
            r <- stats::cor(xjj, yjj)
            set.seed(101 + jj)
            
            
            boot.cjj <- bootcjj(Xjj, MLPB3, up_limit, window_size, boot_rep, n_boot, n.sig)
            
            
            
            
            n1 <- dim(boot.cjj)[1]
            n2 <- dim(boot.cjj)[2]
            
            
            boot.ks.y1 <- 0.5 * log((1 + boot.cjj)/(1 - boot.cjj))
            sd.boot <- apply(boot.ks.y1, 2, mad)
            
            k11 <- apply(boot.ks.y1, 2, function(x) quantile(x, 0.975, na.rm = T))
            k22 <- apply(boot.ks.y1, 2, function(x) quantile(x, 0.5, na.rm = T))
            k33 <- apply(boot.ks.y1, 2, function(x) quantile(x, 0.025, na.rm = T))
            
            
            est.sd <- c(k11, k22, k33, sd.boot)
            
            
            coverage.zero <- sum((0 <= k11) & (k33 <= 0))
            
            coverage.static <- sum((r <= k11) & (k33 <= r))
            
            
            to_save <- list(i, coverage.zero, coverage.static, est.sd, boot.cjj)
            
            median_to_save <- apply(to_save[[5]], 2, median)
            
            
            
            #rd_dir <- file.path(path.out, "Rdata/")
            
            #dir.create(rd_dir, recursive = TRUE, showWarnings = FALSE)
            dir.create(path = paste(output_dir, "/result/", save_file_suffix, .Platform$file.sep, tjj.name[jj], save_file_suffix, 
                sep = ""), showWarnings = FALSE, recursive = TRUE)
            # name.of.subfile <- paste(tjj.name[jj], "_", as.character(regions[kk]), "_", as.character(regions[ll]), save_file_suffix, 
            #     ".RData", sep = "")
            #location.of.subfile <- file.path(rd_dir, name.of.subfile)
            # save(median_to_save, file= location.of.subfile)
            name.of.csvfile <- paste(tjj.name[jj], "_", as.character(regions[kk]), "_", as.character(regions[ll]), save_file_suffix, 
                ".csv", sep = "")
            
            DF <- data.table::data.table(to_save[[5]])
            data.table::fwrite(DF, paste(paste(output_dir, "/result", save_file_suffix, .Platform$file.sep, tjj.name[jj], save_file_suffix, 
                sep = ""), .Platform$file.sep, name.of.csvfile, sep = ""), row.names = F)
            #data.table::fwrite(DF, location.of.subfile)
            
            invisible({
                rm(list = c("coverage.zero", "coverage.static", "est.sd", "boot.cjj"))
                gc()
            })
            
            
        }
        
        
    }
    if(parallel_run){
      doParallel::stopImplicitCluster()
    }
}
