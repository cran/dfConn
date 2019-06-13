#' @description Linear mixed effect model for dynamic functional connectivity.
#' @title Semiparametric modelling - dynamic functional connectivity model
#' @param dataList list, a list of data matrices.
#' @param op_lme list, options constructed by function \code{options_lme}, see \code{options_lme}.
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
#' @importFrom data.table rbindlist
#' @importFrom stats qnorm
#' @importFrom stats quantile
#' @importFrom parallel makePSOCKcluster
#' @details The output of \code{MLPB_boot} is a complete bootstrapping result for the dynamic functional connectivity (dFC) estimates of each pair of region of interest, however, for further analysis, like \code{lmmDyn} and \code{lmmConn}, we need to summarize the bootstrapping by taking the mean or median of the dFC estimate at each time point.
#' @examples 
#' \dontshow{
#' # Example for testing
#' # Assuming user has run \code{MLPB_boot} and has summarize the bootstrapping result
#' 
#' data(MLPB_output_median)
#' dataList <- list(median_1_2 = MLPB_output_median[[1]])
#' subjects <- c('subject1', 'subject2', 'subject3', 'subject4', 'subject5')
#' 
#' # In our demo data, each subject has a scan with a total of 90 time points
#' time.points <- c(1:40,51:90) 
#'                  
#' num.scan <- 2 # Each subject has 2 scans
#' 
#' op <- options_lme(effective_tp = time.points, 
#'                  ntps.per.scan = 40,
#'                  subjects = subjects, 
#'                  num.scan = 2)
#'                  
#' resDyn <- lmmDyn(dataList, op)
#' rm(list = c('subjects', 'MLPB_output_median', 'time.points', 'num.scan', 'ntps.per.scan', 'resDyn'))
#' gc()
#' }
#' \donttest{
#' # Assuming user has run MLPB_boot() and has summarize the bootstrapping result
#' 
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
#'                  
#' resDyn <- lmmDyn(MLPB_output_median, op)
#' }
#' @export
#'


lmmDyn <- function(dataList, op_lme) {
    
    
    # Load input arguments from options list
    subjects <- op_lme$subjects
    eff_time_points <- op_lme$effective_tp
    num.scan <- op_lme$num.scan
    ntps.per.scan <- op_lme$ntps.per.scan
    output_dir <- op_lme$output_dir
    ci_level <- op_lme$ci_level
    cores <- op_lme$cores
    ngrid <- op_lme$ngrid
    numIntKnots <- op_lme$numIntKnots
    seed <- op_lme$seed
    parallel_run <- op_lme$parallel_run
    save_output <- op_lme$save_output
    
    requireNamespace("foreach")
    requireNamespace("stats")
    num.subjects <- length(subjects)
    
    # Precheck
    if (!is.null(seed)) 
        set.seed(seed)
    files.to.run.analysis <- dataList
    
    
    ## Create directory for results
    if (isTRUE(save_output)) {
        if (!is.null(output_dir)) {
          if (!dir.exists(output_dir)) {
            dir.create(output_dir, showWarnings = FALSE)
          }
          
            path1 <- file.path(output_dir, "RData/")
            path2 <- file.path(output_dir, "modelDyn", "lme_models_results/")
            path3 <- file.path(output_dir, "modelDyn", "lme_CI_output", "fit6_REML/")
            path4 <- file.path(output_dir, "modelDyn", "lme_models_results_by_one", "fit6_REML/")
            path5 <- file.path(output_dir, "modelDyn", "lme_models_results_by_one", "rerun_4", "fit6_REML/")
            
            
            dir.create(path1, recursive = T, showWarnings = F)
            dir.create(path2, recursive = T, showWarnings = F)
            dir.create(path3, recursive = T, showWarnings = F)
            dir.create(path4, recursive = T, showWarnings = F)
            dir.create(path5, recursive = T, showWarnings = F)
        } else {
            warning("No output directory being set, the result will not be saved.")
        }
    }
    
    
    
    output_obj <- list() # Empty object initialization
    
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
    output_obj$parallel$is_parallel <- parallel_run
    output_obj$parallel$cores <- cores
    
    # if cores > 1, we consider it to be run in parallel, set up parallel cluster
    if(parallel_run){
      cl <- parallel::makePSOCKcluster(2)
      doParallel::registerDoParallel(cl, cores = cores)
    }
    
    "%dopar%" <- foreach::"%dopar%"
    ijk <- NA
    
    
    output.m6.REML.list <- foreach::foreach(ijk = 1:length(files.to.run.analysis),
                                            .combine = "comb",
                                            .packages = c("doParallel"),
                                            .export = c("ZOSull"),
                                            .init = list(list(), list())) %dopar% {

        data_org <- files.to.run.analysis[[ijk]]
        series_not_NA <- which(is.finite(data_org[, 1]))
        
        if (length(series_not_NA) == num.subjects) {
            missing_to_remove <- 0
        } else {
            
            ind.sub <- c(1:num.subjects)
            missing_to_remove <- ind.sub[!(c(1:num.subjects) %in% series_not_NA)]
        }
        
        
        minus.subj <- num.subjects - length(series_not_NA)
        
        
        
        if (minus.subj != 0) {
          data_org <- data_org[-c(missing_to_remove), ]
        } else {
          data_org <- data_org
        }
        
        num.of.sub <- num.subjects - minus.subj

        data_org <- data_org[, eff_time_points] # Extract the effective dFC trajectories 
        
        name.rdata.save <- names(files.to.run.analysis)[ijk]
        name.rdata.save1 <- name.rdata.save
        
        # dummy data
        
        data.data2 <- data_org
      
        
        if (num.of.sub != num.subjects) {
          # Remove missing scans
          subjects <- subjects[-c(missing_to_remove)]
        } else {
          # No missing scan, do nothing
          subjects <- subjects
        }
        
        
        ###### data for 6 scan
        full.ntps <- ntps.per.scan * num.scan
        subject <- as.vector(t(matrix(rep(subjects, full.ntps), num.of.sub, full.ntps)))
        long.scans.est <- as.vector(t(data_org))
        scan <- c()
        for (i in 1:num.scan) {
            scan <- c(scan, rep(sprintf("s%d", i), ntps.per.scan))
        }
        scan <- rep(scan, num.of.sub)
        group <- as.vector(matrix(rep(c(rep("ConditionA", 0.5 * full.ntps), rep("ConditionB", 0.5 * full.ntps)), num.of.sub), ntps.per.scan * num.scan, num.of.sub))
        time <- as.vector(matrix(rep(rep(c(1:ntps.per.scan), num.scan), num.of.sub), ntps.per.scan * num.scan, num.of.sub))
        
        scan.long <- data.frame(cbind(subject, long.scans.est, time, group, scan))
        colnames(scan.long) <- c("subject", "Corr.est", "time", "group", "scan")
        
        scan.long <- transform(scan.long, Corr.est = as.numeric(as.character(scan.long$Corr.est)), time = as.numeric(as.character(time)))
        group <- factor(scan.long$group)
        
        
        # Create relevant variables:
        y <- scan.long$Corr.est
        x <- scan.long$time
        idnumOrig <- scan.long$subject
        typeIsB <- (scan.long$group == "ConditionB") * 1
        scan <- scan.long$scan
        flavor <- scan.long$group
        
        # Create new (ordered) ID numbers:
        
        idnum <- rep(NA, length(idnumOrig))
        uqID <- unique(idnumOrig)
        for (i in 1:length(uqID)) idnum[idnumOrig == uqID[i]] <- i
        
        
        numObs <- length(x)
        numGrp <- length(unique(idnum))
        
        
        
        
        Xbase <- cbind(rep(1, numObs), x)
        XBvA <- typeIsB * Xbase
        X <- cbind(Xbase, XBvA)
        colnames(X) <- c("Intercept", "time", "Intercept difference", "Slope difference")
        
        
        intKnots <- quantile(unique(x), seq(0, 1, length = numIntKnots + 2))[-c(1, numIntKnots + 2)]
        
        range.x <- c(1.01 * min(x) - 0.01 * max(x), 1.01 * max(x) - 0.01 * min(x))
        
        Zbase <- ZOSull(x, intKnots = intKnots, range.x = range.x)
        ZBvA <- typeIsB * Zbase
        
        Zscan <- kronecker(diag(num.scan * num.of.sub), rep(1, ntps.per.scan))
        
        # Set up Z matrix for group specific curves:
        
        subj.scan <- paste(idnumOrig, scan, sep = "")
        
        # Let up dimension variables:
        
        ncZ <- ncol(Zbase)
        
        # Fit using linear mixed model software:
        
        dummyId <- factor(rep(1, numObs))
        Zblock <- list(dummyId = nlme::pdIdent(~-1 + Zbase), dummyId = nlme::pdIdent(~-1 + ZBvA), idnum = nlme::pdSymm(~1), idnum = nlme::pdIdent(~-1 + scan))
        
        
        
        try({
            fit.6.REML <- nlme::lme(y ~ -1 + X, random = Zblock, method = "REML")
            output.6.REML <- summary(fit.6.REML)
            
            
            
            n1.table <- dim(output.6.REML$tTable)[1]
            n2.table <- dim(output.6.REML$tTable)[2]
            size.fixed <- n1.table * n2.table
            output.6.REML.row <- as.data.frame(matrix(NA, 1, (size.fixed + 3 + 4 + 1 + 1)))
            
            for (ijkm in (1:n1.table)) {
                output.6.REML.row[1, ((1 + (ijkm - 1) * n2.table):(ijkm * n2.table))] <- output.6.REML$tTable[ijkm, ]
            }
            
            
            
            output.6.REML.row[1, size.fixed + 1] <- output.6.REML$AIC
            output.6.REML.row[1, size.fixed + 2] <- output.6.REML$BIC
            output.6.REML.row[1, size.fixed + 3] <- output.6.REML$logLik
            
            output.6.REML.row[1, size.fixed + 4] <- fit.6.REML$sigma
            output.6.REML.row[1, size.fixed + 5] <- unlist(fit.6.REML$modelStruct)[4]
            output.6.REML.row[1, size.fixed + 6] <- unlist(fit.6.REML$modelStruct)[3]
            output.6.REML.row[1, size.fixed + 7] <- unlist(fit.6.REML$modelStruct)[2]
            output.6.REML.row[1, size.fixed + 8] <- unlist(fit.6.REML$modelStruct)[1]
            output.6.REML.row[1, size.fixed + 9] <- name.rdata.save1
            output.6.REML.row[1, size.fixed + 10] <- minus.subj
            
            if (isTRUE(save_output)) 
                write.csv(output.6.REML.row, paste(path4, "outpu6_REML_", name.rdata.save1, ".csv", sep = ""), row.names = FALSE)
            
            ####################################################################################################### 
            sig.epsHat <- fit.6.REML$sigma
            sig.uHat <- sig.epsHat * exp(unlist(fit.6.REML$modelStruct))
            
            lamVal <- (sig.epsHat/sig.uHat)^2
            
            # Set up plotting grids:
            
            xg <- seq(range.x[1], range.x[2], length = ngrid)
            Xg <- cbind(rep(1, ngrid), xg)
            Zg <- ZOSull(xg, intKnots = intKnots, range.x = range.x)
            betaHat <- as.vector(fit.6.REML$coef$fixed)
            uHatBase <- as.vector(fit.6.REML$coef$random[[1]])
            uHatCont <- as.vector(fit.6.REML$coef$random[[2]])
            
            fhatBaseg <- Xg %*% betaHat[1:2] + Zg %*% uHatBase
            Contg <- Xg %*% betaHat[3:4] + Zg %*% uHatCont
            fhatGatg <- fhatBaseg + Contg
            
            
            # Obtain penalty matrix:
            
            sigsq.epsHat <- fit.6.REML$sigma^2
            sig.uHat <- as.numeric(sqrt(sigsq.epsHat) * exp(unlist(fit.6.REML$modelStruct)))
            
            sigsq.Gbl <- sig.uHat[4]^2
            sigsq.Gblcont <- sig.uHat[3]^2
            SigmaHat <- sig.uHat[2]^2
            sigsq.Scan <- sig.uHat[1]^2
            
            
            DmatGbl <- (sigsq.epsHat/sigsq.Gbl) * diag(ncZ)  # u_k
            DmatGblcont <- (sigsq.epsHat/sigsq.Gblcont) * diag(ncZ)  # w_k
            DmatLinSbj <- sigsq.epsHat * kronecker(diag(numGrp), solve(SigmaHat))  #random slope and intercept
            
            DmatScanSbj <- sigsq.epsHat * kronecker(diag(num.of.sub * num.scan), solve(sigsq.Scan))  # random for scan
            dim(DmatScanSbj)
            
            dimVec <- c(4, nrow(DmatGbl), nrow(DmatGblcont), nrow(DmatLinSbj), nrow(DmatScanSbj))
            
            lamMat <- matrix(0, sum(dimVec), sum(dimVec))
            csdV <- cumsum(dimVec)
            lamMat[(csdV[1] + 1):csdV[2], (csdV[1] + 1):csdV[2]] <- DmatGbl
            lamMat[(csdV[2] + 1):csdV[3], (csdV[2] + 1):csdV[3]] <- DmatGblcont
            lamMat[(csdV[3] + 1):csdV[4], (csdV[3] + 1):csdV[4]] <- DmatLinSbj
            lamMat[(csdV[4] + 1):csdV[5], (csdV[4] + 1):csdV[5]] <- DmatScanSbj
            
            # Obtain C matrix:
            
            uqID <- unique(idnum)
            Cmat <- cbind(X, Zbase, ZBvA)
            for (iSbj in 1:numGrp) {
                newCols <- matrix(0, numObs, 1)
                indsCurr <- (1:numObs)[idnum == uqID[iSbj]]
                newCols[indsCurr, ] <- X[indsCurr, 1]
                Cmat <- cbind(Cmat, newCols)
            }
            
            Cmat <- cbind(Cmat, Zscan)
            
        })
        
        try({
            CTC <- crossprod(Cmat)
            fullCovMat <- qr.solve(CTC + lamMat)
            
            # calculate edf according to formula
            S.1 <- tcrossprod(fullCovMat, Cmat)
            diag.S <- matrix(NA, dim(Cmat)[1], 1)
            for (i in (1:dim(Cmat)[1])) {
                diag.S[i, 1] <- sum(Cmat[i, ] * S.1[, i])
            }
            
            
            output.6.REML.row[(length(output.6.REML.row) + 1)] <- sum(diag.S)
            
            
            
            # Find subset of covariance corresponding to the contrast curve:
            
            contInds <- c(3, 4, (ncZ + 5):(2 * ncZ + 4))
            contCovMat <- fullCovMat[contInds, contInds]
            
            # Obtain approximate pointwise 95% confidence limits:
            
            Cg <- cbind(Xg, Zg)
            sdg <- sqrt(sigsq.epsHat) * sqrt(diag(Cg %*% contCovMat %*% t(Cg)))
            lowerg <- Contg - qnorm(ci_level) * sdg
            upperg <- Contg + qnorm(ci_level) * sdg
            
            
            
            # Find subset of covariance corresponding to the contrast curve:
            
            contInds <- c(3, 4, (ncZ + 5):(2 * ncZ + 4))
            contCovMat <- fullCovMat[contInds, contInds]
            
            Condition.A.Inds <- c(1, 2, (csdV[1] + 1):(csdV[2]))
            Condition.A.CovMat <- fullCovMat[Condition.A.Inds, Condition.A.Inds]
            
            
            Condition.B.Inds <- c(1, 2, 3, 4, (csdV[1] + 1):(csdV[3]))
            Condition.B.CovMat <- fullCovMat[Condition.B.Inds, Condition.B.Inds]
            
            sdg.B <- sqrt(sigsq.epsHat) * sqrt(diag(Cg %*% Condition.A.CovMat %*% t(Cg)))
            lowerg.Condition.A <- fhatBaseg - qnorm(ci_level) * sdg.B # Condition A estimate confidence interval lower bound
            upperg.Condition.A <- fhatBaseg + qnorm(ci_level) * sdg.B # Condition A estimate confidence interval upper bound
            
            
            Cg.G <- cbind(Xg, Xg, Zg, Zg)
            sdg.G <- sqrt(sigsq.epsHat) * sqrt(diag(Cg.G %*% Condition.B.CovMat %*% t(Cg.G)))
            lowerg.Condition.B <- fhatGatg - qnorm(ci_level) * sdg.G # Condition B estimate confidence interval lower bound
            upperg.Condition.B <- fhatGatg + qnorm(ci_level) * sdg.G # Condition B estimate confidence interval upper bound
            
            
            est.CI <- matrix(NA, 9, ngrid)
            est.CI[1, ] <- fhatBaseg 
            est.CI[2, ] <- Contg
            est.CI[3, ] <- fhatGatg
            est.CI[4, ] <- lowerg
            est.CI[5, ] <- upperg
            
            est.CI[6, ] <- lowerg.Condition.A
            est.CI[7, ] <- upperg.Condition.A
            
            est.CI[8, ] <- lowerg.Condition.B
            est.CI[9, ] <- upperg.Condition.B
            
            if (isTRUE(save_output)) {
                write.csv(est.CI, paste(path3, "modelDyn_REML_lme", "_", name.rdata.save1, ".csv", sep = ""), row.names = FALSE)
                save(output.6.REML.row, file = file.path(path1, paste("output_fit6_REML", "_", name.rdata.save1, ".RData", sep = "")))
                write.csv(output.6.REML.row, paste(path5, "outpu6_REML_", name.rdata.save1, ".csv", sep = ""), row.names = FALSE)
                
            }
            # print(name.rdata.save1)
            est.CI <- as.data.frame(est.CI)
            
            list(est.CI, output.6.REML.row)
        })
        
    }
    # print(output.m6.REML) Combine output.ml.REML data frame
    output_obj$output_by_row <- output.m6.REML.list[[2]]
    output_obj$modelDyn_results <- output.m6.REML.list[[2]]
    output_obj$est_CI <- output.m6.REML.list[[1]]
    names(output_obj$output_by_row) <- names(files.to.run.analysis)
    names(output_obj$est_CI) <- names(files.to.run.analysis)
    output.m6.REML <- data.table::rbindlist(output.m6.REML.list[[2]])
    
    if (is.null(dim(output.m6.REML))) 
        output.m6.REML <- as.data.frame(matrix(output.m6.REML, 1, length(output.m6.REML)))
    
    colnames(output.m6.REML) <- c("XIntercept_Value", "XIntercept_Std.Error", "XIntercept_DF", "XIntercept_t-value", "XIntercept_p-value", "Xtime_Value", "Xtime_Std.Error", 
        "Xtime_DF", "Xtime_t-value", "Xtime_p-value", "Xintercept_difference_Value", "Xintercept_difference_Std.Error", "Xintercept_difference_DF", "Xintercept_difference_t-value", 
        "Xintercept_difference_p-value", "Xslope_difference_Value", "Xslope_difference_Std.Error", "Xslope_difference_DF", "Xslope_difference_t-value", "Xslope_difference_p-value", 
        "AIC", "BIC", "logLik", "sigma_eps", "sigma_u_unstruct", "sigma_w_unstruct", "sigma_b0_unstruct", "sigma_a_unstruct", "comparison", "missing_series", "edf")
    
    output_obj$modelDyn_results <- as.data.frame(output.m6.REML)
    
    if (isTRUE(save_output)) {
        write.csv(output.m6.REML, paste(path2, "modelDyn_REML_lme_all.csv", sep = ""), row.names = FALSE)
    }
  
    # Clean up cluster for parallel computing
    if(parallel_run){
      doParallel::stopImplicitCluster()
    }

    return(output_obj)
    
}



