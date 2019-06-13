#' @description Linear mixed effect modeling for static functional connectivity 
#' 
#' @title Semiparametric modelling - Static
#' @param dataList list, list of data matrices .
#' @param op_lme list, options constructed by function \code{options_lme}, see \code{options_lme}.
#' @importFrom utils write.csv
#' @importFrom parallel makePSOCKcluster
#' @details 
#' \code{lmmDyn} and \code{lmmConn}, we need to summarize the bootstrapping by taking the mean or median of the dFC estimate at each time point.
#' @return An object of list containing all the information of the static functional connectivity linear mixed model. It stores the information of model parameters, input data, and model results. 
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
#' @examples 
#' 
#' # Assuming user has run MLPB_boot() and has summarized the bootstrapping results 
#' # by calculating the median or mean.
#' 
#' ## Example with 5 subjects bootstrap-based static functional connectivity estimates data included
#' ## Each subject's data has 731 time points in total, which includes 6 scans and 105 effective 
#' ## time point for each scan. 
#' ## 3 regions of interest (ROIs) comparision pairs are selected
#' 
#' data(MLPB_output_median)
#' subjects <- c('subject1', 'subject2', 'subject3', 'subject4', 'subject5')
#' 
#' # In our demo data, each subject has a scan with a total of 750 time points
#' time.points <- c(1:105, 126:230, 251:355,
#'                  376:480, 501:605, 626:730) 
#'                  
#' num.scan <- 6 # Each subject has 6 scans
#' ntps.per.scan <- 105 # Each scan has 105 time points
#' 
#' op <- options_lme(effective_tp = time.points, 
#'                  ntps.per.scan = ntps.per.scan,
#'                  subjects = subjects, 
#'                  num.scan = num.scan)
#'                  
#' resConn <- lmmConn(MLPB_output_median, op)
#' \dontshow{
#' rm(list = c('subjects', 'MLPB_output_median', 'time.points', 'num.scan', 'ntps.per.scan','resConn'))
#' gc()}
#' @export
#'




lmmConn <- function(dataList, op_lme) {
    
  
  # Load input arguments from options list
  subjects <- op_lme$subjects
  eff_time_points <- op_lme$effective_tp
  num.scan <- op_lme$num.scan
  ntps.per.scan <- op_lme$ntps.per.scan
  output_dir <- op_lme$output_dir
  ci_level <- op_lme$ci_level
  cores <- op_lme$cores
  ngrid <- op_lme$ngrid
  seed <- op_lme$seed
  parallel_run <- op_lme$parallel_run
  save_output <- op_lme$save_output
  
    requireNamespace("foreach")
    num.subjects <- length(subjects)
    

    if (!is.null(seed)) 
        set.seed(seed)
    
    files.to.run.analysis <- dataList
    
    if (isTRUE(save_output)) {
        if (!is.null(output_dir)) {
          if (!dir.exists(output_dir)) {
            dir.create(output_dir)
          }
          
            path1 <- file.path(output_dir, "RData/")
            path2 <- file.path(output_dir, "modelConn", "lme_models_results/")
            path3 <- file.path(output_dir, "modelConn", "lme_CI_output", "fit0_REML/")
            path4 <- file.path(output_dir, "modelConn", "lme_models_results_by_one", "fit0_REML/")
            path5 <- file.path(output_dir, "modelConn", "lme_models_results_by_one", "rerun_4", "fit0_REML/")
            
            
            dir.create(path1, recursive = TRUE, showWarnings = FALSE)
            dir.create(path2, recursive = TRUE, showWarnings = FALSE)
            dir.create(path3, recursive = TRUE, showWarnings = FALSE)
            dir.create(path4, recursive = TRUE, showWarnings = FALSE)
            dir.create(path5, recursive = TRUE, showWarnings = FALSE)
        }
    }
    
    
    output_obj <- list()
    # save parameters
    output_obj$params <- list()
    output_obj$params$dataList <- dataList
    output_obj$params$subjects <- subjects
    output_obj$params$effective_tps <- eff_time_points
    output_obj$params$num_scans <- num.scan
    output_obj$params$ntp_per_scan <- ntps.per.scan
    output_obj$params$ci_level <- ci_level
    output_obj$params$ngrid <- ngrid
    
    
    
    # Parallel computing
    output_obj$parallel <- list()
    output_obj$parallel$is_parallel <- ifelse(cores > 1, T, F)
    output_obj$parallel$cores <- cores
    
    if(cores>1){
      cl <- parallel::makePSOCKcluster(2)
      doParallel::registerDoParallel(cl, cores = cores)
    }

    
    "%dopar%" <- foreach::"%dopar%"
    ijk <- NA
    # print(files.to.run.analysis)
    output.m0.list <- foreach::foreach(ijk = 1:length(files.to.run.analysis), .combine = "comb", .packages = "doParallel", .init = list(list(), list())) %dopar% {
        # for(ijk in 1:length(files.to.run.analysis)){
        
        # load data
        
        data.data1 <- files.to.run.analysis[[ijk]]
        series_not_NA <- which(is.finite(data.data1[, 1]))
        
        if (length(series_not_NA) == num.subjects) {
            num_NA_series <- 0
        } else {
            ind.sub <- c(1:num.subjects)
            num_NA_series <- ind.sub[!(c(1:num.subjects) %in% series_not_NA)]
        }
        
        
        minus.subj <- num.subjects - length(series_not_NA)
        
        
        
        if (minus.subj != 0) {
            data.data1 <- data.data1[-c(num_NA_series), ]
        } else {
            data.data1 <- data.data1
        }
        
        num.of.sub <- num.subjects - minus.subj
        
        data.data1 <- data.data1[, eff_time_points]
        name.rdata.save <- names(files.to.run.analysis)[ijk]
        name.rdata.save1 <- name.rdata.save
        
        # dummy data
        data.data2 <- data.data1
        
        if (num.of.sub != num.subjects) {
          subjects <- subjects[-c(num_NA_series)]
        } else {
          subjects <- subjects
        }
        
        ###### data for 6 scans
        full.ntps <- ntps.per.scan * num.scan
        subject <- as.vector(t(matrix(rep(subjects, full.ntps), num.of.sub, full.ntps)))
        long.scans.est <- as.vector(t(data.data1))
        scan <- c()
        for (i in 1:num.scan) {
            scan <- c(scan, rep(sprintf("s%d", i), ntps.per.scan))
        }
        scan <- rep(scan, num.of.sub)
        
        # Assuming that each condition is equally distributed
        group <- as.vector(matrix(rep(c(rep("Condition1", 0.5 * full.ntps), rep("Condition2", 0.5 * full.ntps)), num.of.sub), full.ntps, num.of.sub))
        
        time <- as.vector(matrix(rep(rep(c(1:ntps.per.scan), num.scan), num.of.sub), full.ntps, num.of.sub))
        
        scan.long <- data.frame(cbind(subject, long.scans.est, time, group, scan))
        # print(scan.long)
        colnames(scan.long) <- c("subject", "Corr.est", "time", "group", "scan")
        
        scan.long <- transform(scan.long, Corr.est = as.numeric(as.character(scan.long$Corr.est)), time = as.numeric(as.character(time)))
        
        group <- factor(scan.long$group)
        
        
        # Create relevant variables:
        y <- scan.long$Corr.est
        x <- scan.long$time
        idnumOrig <- scan.long$subject
        typeIsB <- (scan.long$group == "Condition2") * 1
        scan <- scan.long$scan
        flavor <- scan.long$group
        
        # Create new (ordered) ID numbers:
        
        idnum <- rep(NA, length(idnumOrig))
        uqID <- unique(idnumOrig)
        for (i in 1:length(uqID)) idnum[idnumOrig == uqID[i]] <- i
        
        
        numObs <- length(x)
        numGrp <- length(unique(idnum))
        
        
        Xbase <- cbind(rep(1, numObs))
        XBvA <- typeIsB * Xbase
        X <- cbind(Xbase, XBvA)
        colnames(X) <- c("Intercept", "Intercept difference")
        
        
        Zscan <- kronecker(diag(num.scan * num.of.sub), rep(1, ntps.per.scan))
        subj.scan <- paste(idnumOrig, scan, sep = "")
        
        # Fit using linear mixed model software:
        
        dummyId <- factor(rep(1, numObs))
        Zblock <- list(idnum = nlme::pdIdent(~1), idnum = nlme::pdIdent(~-1 + scan))
        
        try({
            fit.0 <- nlme::lme(y ~ -1 + X, random = Zblock, control = nlme::lmeControl(opt = "optim"), method = "REML")
            
            output.0 <- summary(fit.0)
            
            n1.table <- dim(output.0$tTable)[1]
            n2.table <- dim(output.0$tTable)[2]
            size.fixed <- n1.table * n2.table
            output.0.row <- matrix(NA, 1, (size.fixed + 3 + 2 + 1 + 1))  # 3=(AIC, BIC, logLIK); 2=(b_i0, a_is); 1=name_of_comaprison, 1=number of people with all the data
            
            for (ijkm in (1:n1.table)) {
                output.0.row[1, ((1 + (ijkm - 1) * n2.table):(ijkm * n2.table))] <- output.0$tTable[ijkm, ]
            }
            
            # fit.0$sigma*exp(unlist(fit.0$modelStruct))
            
            output.0.row[size.fixed + 1] <- output.0$BIC
            output.0.row[size.fixed + 2] <- output.0$AIC
            output.0.row[size.fixed + 3] <- output.0$logLik
            
            output.0.row[size.fixed + 4] <- fit.0$sigma
            output.0.row[size.fixed + 5] <- unlist(fit.0$modelStruct)[2]
            output.0.row[size.fixed + 6] <- unlist(fit.0$modelStruct)[1]
            
            output.0.row[size.fixed + 7] <- name.rdata.save1
            output.0.row[size.fixed + 8] <- minus.subj
            
            if (isTRUE(save_output)) {
                utils::write.csv(output.0.row, paste(path4, "outputConn_REML_", name.rdata.save1, ".csv", sep = ""), row.names = FALSE)
            }
            
            
            ##### 
            sig.epsHat <- fit.0$sigma
            sig.uHat <- sig.epsHat * exp(unlist(fit.0$modelStruct))
            
            lamVal <- (sig.epsHat/sig.uHat)^2
            
            # Set up plotting grids:
            Xg <- cbind(rep(1, ngrid))
            
            betaHat <- as.vector(fit.0$coef$fixed)
            
            fhatBaseg <- Xg %*% betaHat[1]
            Contg <- Xg %*% betaHat[2]
            fhatGatg <- fhatBaseg + Contg
            
            
            # Obtain penalty matrix:
            
            sigsq.epsHat <- fit.0$sigma^2
            sig.uHat <- as.numeric(sqrt(sigsq.epsHat) * exp(unlist(fit.0$modelStruct)))
            
            SigmaHat <- sig.uHat[2]^2
            sigsq.Scan <- sig.uHat[1]^2
            
            DmatLinSbj <- sigsq.epsHat * kronecker(diag(numGrp), solve(SigmaHat))  #random slope and intercept
            
            DmatScanSbj <- sigsq.epsHat * kronecker(diag(num.of.sub * num.scan), solve(sigsq.Scan))  # random for scan
            dimVec <- c(2, nrow(DmatLinSbj), nrow(DmatScanSbj))
            
            lamMat <- matrix(0, sum(dimVec), sum(dimVec))
            csdV <- cumsum(dimVec)
            lamMat[(csdV[1] + 1):csdV[2], (csdV[1] + 1):csdV[2]] <- DmatLinSbj
            lamMat[(csdV[2] + 1):csdV[3], (csdV[2] + 1):csdV[3]] <- DmatScanSbj
            
            
            # Obtain C matrix:
            
            uqID <- unique(idnum)
            Cmat <- cbind(X)
            for (iSbj in 1:numGrp) {
                newCols <- matrix(0, numObs, 1)
                indsCurr <- (1:numObs)[idnum == uqID[iSbj]]
                newCols[indsCurr, ] <- X[indsCurr, 1]
                Cmat <- cbind(Cmat, newCols)
            }
            
            Cmat <- cbind(Cmat, Zscan)
            
        })
        
        # print(lamMat)
        try({
            # print(Cmat)
            CTC <- crossprod(Cmat)
            fullCovMat <- solve(CTC + lamMat)
            S.1 <- tcrossprod(fullCovMat, Cmat)
            diag.S <- matrix(NA, dim(Cmat)[1], 1)
            for (i in (1:dim(Cmat)[1])) {
                diag.S[i, 1] <- sum(Cmat[i, ] * S.1[, i])
            }
            
            
            output.0.row[(length(output.0.row) + 1)] <- sum(diag(diag.S))
            
            if (isTRUE(save_output)) {
                utils::write.csv(output.0.row, file.path(path5, paste("outputConn_REML_", name.rdata.save1, ".csv", sep = "")), row.names = FALSE)
            }
            
            # Find subset of covariance corresponding to the contrast curve:
            
            contInds <- c(2)
            contCovMat <- fullCovMat[contInds, contInds]
            
            # Obtain approximate pointwise 95% confidence limits:
            
            Cg <- cbind(Xg)
            sdg <- sqrt(sigsq.epsHat) * sqrt(diag(Cg %*% contCovMat %*% t(Cg)))
            lowerg <- Contg - qnorm(ci_level) * sdg
            upperg <- Contg + qnorm(ci_level) * sdg
            
            
            # Find subset of covariance corresponding to the contrast curve:
            
            contInds <- c(2)
            contCovMat <- fullCovMat[contInds, contInds]
            
            Condition1Inds <- c(1)
            Condition1CovMat <- fullCovMat[Condition1Inds, Condition1Inds]
            
            
            GatInds <- c(1, 2)
            GatCovMat <- fullCovMat[GatInds, GatInds]
            
            sdg.B <- sqrt(sigsq.epsHat) * sqrt(diag(Cg %*% Condition1CovMat %*% t(Cg)))
            lowerg.Condition1 <- fhatBaseg - qnorm(ci_level) * sdg.B
            upperg.Condition1 <- fhatBaseg + qnorm(ci_level) * sdg.B
            
            
            Cg.G <- cbind(Xg, Xg)
            sdg.G <- sqrt(sigsq.epsHat) * sqrt(diag(Cg.G %*% GatCovMat %*% t(Cg.G)))
            lowerg.Gat <- fhatGatg - qnorm(ci_level) * sdg.G
            upperg.Gat <- fhatGatg + qnorm(ci_level) * sdg.G
            
            
            est.CI <- matrix(NA, 9, ngrid)
            est.CI[1, ] <- fhatBaseg
            est.CI[2, ] <- Contg
            est.CI[3, ] <- fhatGatg
            est.CI[4, ] <- lowerg
            est.CI[5, ] <- upperg
            
            
            est.CI[6, ] <- lowerg.Condition1
            est.CI[7, ] <- upperg.Condition1
            
            est.CI[8, ] <- lowerg.Gat
            est.CI[9, ] <- upperg.Gat
            
            if (isTRUE(save_output)) {
                utils::write.csv(est.CI, file.path(path3, paste("modelConn_REML_lme", "_", name.rdata.save1, ".csv", sep = "")), row.names = FALSE)
                save(output.0, file = file.path(path1, paste("output_fit0_REML", "_", name.rdata.save1, ".RData", sep = "")))
            }
            
            
            ##### print(name.rdata.save1)
            output.0.row <- as.data.frame(matrix(output.0.row, 1, length(output.0.row)))
            est.CI <- as.data.frame(est.CI)
            
            list(est.CI, output.0.row)
        })
    }
    output_obj$output_by_row <- output.m0.list[[1]]
    output_obj$modelConn_results <- output.m0.list[[2]]
    print(output.m0.list[[2]])
    names(output_obj$modelConn_results) <- names(files.to.run.analysis)
    output.m0 <- data.table::rbindlist(output.m0.list[[2]])
    
    if (is.null(dim(output.m0))) 
        output.m0 <- as.data.frame(matrix(output.m0, 1, length(output.m0)))
    
    colnames(output.m0) <- c("XIntercept_Value", "XIntercept_Std.Error", "XIntercept_DF", "XIntercept_t-value", "XIntercept_p-value", "Xintercept_difference_Value", "Xintercept_difference_Std.Error", 
        "Xintercept_difference_DF", "Xintercept_difference_t-value", "Xintercept_difference_p-value", "AIC", "BIC", "logLik", "sigma_eps", "sigma_b0_unstruct", "sigma_a_unstruct", 
        "comparison", "missing_series", "edf")
    
    if (isTRUE(save_output)) {
        write.csv(output.m0, file.path(path2, paste("modelConn_REML_lme_all.csv", sep = "")), row.names = FALSE)
    }
    #class(output_obj) <- "dFClmm"
    if (cores>1){
      doParallel::stopImplicitCluster()
    }

    return(output_obj)
}

