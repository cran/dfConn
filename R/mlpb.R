#' @description Multivariate Linear Process Bootstrap (MLPB) method to assess the uncertainty in dynamic functional connectivity (dFC) by providing its confidence band.
#' @title  Multivariate Linear Process Bootstrap Method
#' @param dataList list, subject-specific data.
#' @param output_dir character, output directory for bootstrap samples.
#' @param rois vector of integers, specify the list of regions of interests
#' @param boot_rep integer, bootstrapping repetition times
#' @param subset.subject vector of character, which subject to run, specify the index in \code{dataList}, for example, c(1, 2, ...) 
#' @param window_size integer, window size for sliding window technique
#' @param timepoints integer, number of timepoints in total
#' @param number_of_intervals integer, number of intervals in sliding window technique
#' @param n_boot integer, number of bootstrap sample to be generated in MLPB
#' @param cores integer, number of cores to register for parallel execution, if set to 1 a sequential job will be run
#' @param save_file_suffix character, suffix of output files, treated as a labels.
#' @importFrom stats cor
#' @importFrom stats mad
#' @importFrom stats median
#' @importFrom parallel makePSOCKcluster
#' @importFrom gtools combinations
#' @importFrom Rcpp evalCpp
#' @import Rcpp
#' @details The \code{dataList} parameter is a list of matrices which contains time series data of each region of interest (ROI). Output directory is required here because the results to be generate is massive. 
#' @references Kudela et al. (2017) NeuroImage 10.1016/j.neuroimage.2017.01.056
#' (\href{https://www.ncbi.nlm.nih.gov/pubmed/28132931}{PubMed})
#' @examples
#' \dontshow{
#' # Load sample data
#' data(fMRI_dataList_shrinked)
#' MLPB_boot(fMRI_dataList_shrinked, output_dir=tempdir(), rois = 1:2, timepoints = 5, window_size=5)
#' rm(list=c('fMRI_dataList_shrinked'))
#' gc()
#' }
#' \donttest{
#' 
#' # Load sample data
#' 
#' data(fMRI_dataList)
#' 
#' MLPB_boot(fMRI_dataList, output_dir = tempdir(), 
#'           rois = c(54,191,235), 
#'           timepoints = 750)
#' }
#' 
#' @export
#' @keywords bootstarpping
#' 

MLPB_boot <- function(dataList, output_dir, rois, timepoints, subset.subject = NULL, save_file_suffix = "", window_size = 20, 
    number_of_intervals = 1, boot_rep = 250, n_boot = 1, cores = 1) {
    
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
    
    for (jj in subset.subject) {
        regions <- rois
        tjj <- dataList[[jj]]
        tjj <- data.matrix(tjj)
        tjj <- tjj[regions, ]
        wywal.NAN1 <- list()
        
        if (is.null(dim(tjj))) {
            tjj <- t(as.matrix(tjj, 1, timepoints))
        }
        
        for (ipk in 1:timepoints) {
            wywal.NAN1[[ipk]] <- which(!is.finite(tjj[, ipk]))
        }
        
        
        wywal.NAN <- if (length(unique(unlist(wywal.NAN1))) != 0) {
            unique(unlist(wywal.NAN1))
        } else {
            0
        }
        
        
        if (length(wywal.NAN) != 1) {
            tjj1 <- tjj[-wywal.NAN, ]
        } else {
            if (wywal.NAN != 0) {
                tjj1 <- tjj[-wywal.NAN, ]
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
        
        
        
        a.rozmiar <- rozmiar[jj, 1]
        num.comp <- a.rozmiar * (a.rozmiar + 1)/2 - a.rozmiar
        comparison <- matrix(NA, num.comp, 2)
        
        
        path.out <- file.path(output_dir, paste(tjj.name[jj], save_file_suffix, "/", sep = ""))
        
        
        
        
        ############################################ 
        
        go.through <- c(1:length(rois) - 1)
        go.through <- go.through[!(go.through %in% wywal.NAN)]
        
        iii <- 1
        for (kk in go.through) {
            
            
            doll <- max(2 - kk + 1, kk + 1)
            go.through2 <- c(doll:(max(go.through) + 1))
            go.through2 <- go.through2[!(go.through2 %in% wywal.NAN)]
            
            for (ll in go.through2) {
                
                comparison[iii, 1] <- kk
                comparison[iii, 2] <- ll
                iii <- iii + 1
            }
        }
        
        
        ############# change window size 20, mad instead of sd
        n.sig <- timepoints  #number of time points
        wind <- "w20"
        up_limit <- n.sig - window_size + 1
        
        upper_limit <- up_limit
        
        
        cov.zero <- cov.static <- vector()
        boot.cjj <- matrix(NA, nrow = boot_rep, ncol = (up_limit))
        
        
        coverage.zero <- vector()
        coverage.static <- vector()
        est.sd <- vector()
        sd.boot <- matrix()
        
        
        cat(paste("Bootstrapping on ", tjj.name[jj], "\n"))
        
        if(cores > 1){
          cl <- parallel::makePSOCKcluster(2)
          doParallel::registerDoParallel(cl, cores = cores)
        }

        ############################### parallel
        
        '%dopar%' <- foreach::'%dopar%'
        
        # Check files exists or not, returns a list of files that need to be generated file2run <- function(target_folder,
        # comparison, tjj.name, jj){ n <- nrow(comparison) file.existed <- list.files(target_folder) res <- c() for (i in 1:n){ ik
        # <- comparison[i,1] ij <- comparison[i,2] kk<-ik ll<-ij name.of.subfile <-
        # paste(tjj.name[jj],'_',as.character(regions[kk]),'_', as.character(regions[ll]),save_file_suffix,'_medians.RData' ,sep='')
        # if (!name.of.subfile %in% file.existed){ res <- c(res, i) } } return(res) }
        
        # res <- file2run(path.out, comparison, tjj.name, jj) doParallel::registerDoParallel(cores = 3) print(comparison)
        
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
            
            
            
            rd_dir <- file.path(path.out, "Rdata/")
            
            dir.create(rd_dir, recursive = TRUE, showWarnings = FALSE)
            dir.create(path = paste(output_dir, "/result", save_file_suffix, .Platform$file.sep, tjj.name[jj], save_file_suffix, 
                sep = ""), showWarnings = FALSE, recursive = TRUE)
            name.of.subfile <- paste(tjj.name[jj], "_", as.character(regions[kk]), "_", as.character(regions[ll]), save_file_suffix, 
                ".RData", sep = "")
            location.of.subfile <- file.path(rd_dir, name.of.subfile)
            # save(median_to_save, file= location.of.subfile)
            name.of.csvfile <- paste(tjj.name[jj], "_", as.character(regions[kk]), "_", as.character(regions[ll]), save_file_suffix, 
                ".csv", sep = "")
            
            DF <- data.table::data.table(to_save[[5]])
            data.table::fwrite(DF, paste(paste(output_dir, "/result", save_file_suffix, .Platform$file.sep, tjj.name[jj], save_file_suffix, 
                sep = ""), .Platform$file.sep, name.of.csvfile, sep = ""), row.names = F)
            data.table::fwrite(DF, location.of.subfile)
            
            invisible({
                rm(list = c("coverage.zero", "coverage.static", "est.sd", "boot.cjj"))
                gc()
            })
            
            
        }
        
        
    }
    if(cores>1){
      doParallel::stopImplicitCluster(cl)
    }
}
