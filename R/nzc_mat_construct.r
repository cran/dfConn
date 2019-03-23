# description For creation of Non-zero coverage matrix, including cond1, cond2 and difference between two flavors.  title Non-zero coverage Correlation Matrix param
# mat_plot_object matrix-to-plot object param direction c('pos', 'neg') is it a positive lag (t > 0) or a negative lag (t < 0) export keywords Non-zero coverage, matrix
# return return a list of correlation matrices

matrix_construction <- function(mat_plot_object) {
    
    requireNamespace("stats", "utils")
    
    data.plot <- mat_plot_object[["model_results"]]
    tbl.list <- mat_plot_object[["tbl.list"]]
    
    
    roi_names <- mat_plot_object[["roi_names"]]
    roi <- mat_plot_object[["roi"]]
    alpha <- mat_plot_object[["alpha"]]
    roi_char <- as.character(roi)
    cond1.thres <- mat_plot_object[["thres.cond1"]]
    cond2.thres <- mat_plot_object[["thres.cond2"]]
    thres.diff <- mat_plot_object[["thres.diff"]]
    ntps.per.scan <- mat_plot_object[["ntp_per_scan"]]
    cat("Mat object attributes loaded.\n")
    # print(model_results_all)
    
    
    # Read data
    data.order <- yeo.order
    names.our.data <- colnames(data.plot)
    
    ######################################## 
    part1 <- vector()
    part2 <- vector()
    part3 <- vector()
    part4 <- vector()
    
    for (i in c(1:nrow(data.plot))) {
        
        name1 <- unlist(strsplit(data.plot$comparison[i], "_"))
        
        part1[i] <- unlist(strsplit(name1, "_"))[1]
        part2[i] <- unlist(strsplit(name1, "_"))[2]
        part3[i] <- unlist(strsplit(name1, "_"))[3]
        part4[i] <- unlist(strsplit(name1, "_"))[4]
        
    }
    
    data.plot$loc1 <- part2
    data.plot$loc2 <- part3
    data.plot$lp <- part4
    
    ######################################## 
    voi1 <- data.plot[, "missing_series"]
    voi <- voi1
    
    # Smaller Matrix should be 7 by 7
    value.miss <- matrix(NA, 278, 278)
    
    for (i in c(1:choose(length(roi), 2))) {
        

        kk <- as.numeric(data.plot$loc1[i])
        ll <- as.numeric(data.plot$loc2[i])
        
        value.miss[kk, ll] <- as.numeric(voi[i])
        value.miss[ll, kk] <- as.numeric(voi[i])
        
    }
    
    diag(value.miss) <- 0
    cat("Mat object data cleaned.\n")
    
    # yeo.order
    dane.plot.miss <- value.miss[, order(data.order[, 2])]
    dane.plot.miss.yeo <- dane.plot.miss[order(data.order[, 2]), ]
    
    
    ####################################### 
    coverage.zero.diff <- tbl.list$coverage.zero.cond.diff
    coverage.zero.cond1 <- tbl.list$coverage.zero.cond.1
    coverage.zero.cond2 <- tbl.list$coverage.zero.cond.2
    
    coverage.zero.diff[, 2] <- 1 - coverage.zero.diff[, 2]
    coverage.zero.cond1[, 2] <- 1 - coverage.zero.cond1[, 2]
    coverage.zero.cond2[, 2] <- 1 - coverage.zero.cond2[, 2]
    
    
    loc_diff_roiPair <- extract_roi_pair(coverage.zero.diff[, 1])
    loc1.diff <- loc_diff_roiPair$roi1
    loc2.diff <- loc_diff_roiPair$roi2
    #print()
    # lp.diff <- unlist(strsplit(coverage.zero.diff[, 1], '_'))[seq(3, (length(coverage.zero.diff[, 1]) * 5), 5)]
    
    loc_cond1_roiPair <- extract_roi_pair(coverage.zero.cond1[, 1])
    loc1.cond1 <- loc_cond1_roiPair$roi1
    loc2.cond1 <- loc_cond1_roiPair$roi2
    # lp.cond1 <- unlist(strsplit(coverage.zero.cond1[, 1], '_'))[seq(3, (length(coverage.zero.cond1[, 1]) * 5), 5)]
    
    loc_cond2_roiPair <- extract_roi_pair(coverage.zero.cond2[, 1])
    loc1.cond2 <- loc_cond1_roiPair$roi1
    loc2.cond2 <- loc_cond1_roiPair$roi2
    # lp.cond2 <- unlist(strsplit(coverage.zero.cond2[, 1], '_'))[seq(3, (length(coverage.zero.cond2[, 1]) * 5), 5)]
    
    
    numberr <- length(roi)
    cat("Non-zero coverage information loaded.\n")
    
    ############################# difference ####################
    
    ommm <- 0
    
    
    coverage.zero.diff[is.na(coverage.zero.diff[, 2]), 2] <- 0
    
    prop.res.chisq <- matrix(NA, nrow(coverage.zero.diff), 1)
    prop.res.chisq.pvalue.diff <- matrix(NA, nrow(coverage.zero.diff), 1)
    val.prop.test <- vector()
    
    for (imk in 1:nrow(coverage.zero.diff)) {
        
        val.prop.test[imk] <- round((coverage.zero.diff[imk, 2]) * ntps.per.scan, 0)
        prop.res <- stats::prop.test(val.prop.test[imk], ntps.per.scan, alternative = "greater", p = thres.diff)
        prop.res.chisq[imk, 1] <- prop.res$statistic
        prop.res.chisq.pvalue.diff[imk, 1] <- prop.res$p.value
    }
    
    coverage.nonzero.diff.chsq <- cbind(coverage.zero.diff[, 1], (coverage.zero.diff[, 2]), prop.res.chisq, loc1.diff, loc2.diff)
    
    colnames(coverage.nonzero.diff.chsq) <- c("comparison", "Significant coverage for Diff", "Chi-Square", "loc1.diff", "loc2.diff")
    
    
    # cat('dummy matrix \n')
    
    dummy.m <- matrix(NA, 278, 278)
    
    for (i in 1:278) {
        for (j in i:278) {
            dummy.m[i, j] <- paste(i, "_", j, sep = "")
            dummy.m[j, i] <- paste(i, "_", j, sep = "")
        }
    }
    
    dummy.m1 <- dummy.m[, order(data.order[, 2])]
    dummy.m2 <- dummy.m1[order(data.order[, 2]), ]
    
    Cerebellum.loc.id <- dummy.m2[249:278, ]
    
    Cerebellum.loc.id <- unique(as.vector(Cerebellum.loc.id))
    
    ######################################## 
    loc1_loc2.diff <- paste(coverage.nonzero.diff.chsq[, "loc1.diff"], "_", coverage.nonzero.diff.chsq[, "loc2.diff"], sep = "")
    
    ######################################## 
    
    cerebellum <- data.order[order(data.order[, 2]), ][249:278, ]
    
    p.val.diff.without.cerebellum <- prop.res.chisq.pvalue.diff[!(loc1_loc2.diff %in% Cerebellum.loc.id)]
    
    ################ Choose alpha level ##############
    order.p.values.diff <- p.val.diff.without.cerebellum[order(p.val.diff.without.cerebellum)]
    m <- length(order.p.values.diff)
    k <- seq(1, m, 1)
    
    # BH dependent correction https://en.wikipedia.org/wiki/False_discovery_rate
    c_m <- sum(1/c(1:m))
    
    # c_m <- log(m)+ 0.5772156649+ 1/(2*m)
    BH.corr.dep <- (k/(m * c_m)) * alpha
    max.BH.dep.diff <- max(order.p.values.diff[order.p.values.diff <= BH.corr.dep])
    
    significant.BH.dep <- rep(2, m)
    significant.BH.dep[order.p.values.diff <= BH.corr.dep] <- 1
    significant.BH.dep[order.p.values.diff > BH.corr.dep] <- 0
    # table(significant.BH.dep)
    
    row.num.BH.dep <- which(prop.res.chisq.pvalue.diff == max.BH.dep.diff)
    chisq.BH.dep <- abs(prop.res.chisq[row.num.BH.dep])[1]
    
    
    ################ 
    
    coverage.nonzero.diff.chsq[as.numeric(prop.res.chisq.pvalue.diff) > max.BH.dep.diff, 2] <- 0
    
    ################ 
    
    coverage.zero.diff.all <- coverage.nonzero.diff.chsq
    
    loc.diff <- coverage.zero.diff.all[(coverage.zero.diff.all[, 4] %in% roi_char), ]
    coverage.zero.roi.diff <- loc.diff[(loc.diff[, 5] %in% roi_char), ]
    
    
    ############################################################################## select all that significant ##########################
    
    sel.sig.pics <- coverage.zero.roi.diff[coverage.zero.roi.diff[, 2] != "0", ]
    if (is.null(sel.sig.pics)) {
        stop("No significant comparisons, exiting...")
    }else if(is.vector(sel.sig.pics)){
      temp_mat <- as.data.frame(matrix(NA, 1, 5))
      temp_mat[1,] <- sel.sig.pics
      sel.sig.pics <- temp_mat
      names(sel.sig.pics) <- c("comparison", "Significant coverage for Diff", "Chi-Square", "loc1.diff", "loc2.diff")
    }
    nrow_sel.sig.pics <- nrow(sel.sig.pics)
    roi1.temp <- as.integer(sel.sig.pics[, 4])
    roi2.temp <- as.integer(sel.sig.pics[, 5])
    
    
    
    roi.temp <- unique(c(as.numeric(roi1.temp), as.numeric(roi2.temp)))
    
    
    
    
    ################################################################# 
    lp.roi <- length(roi)
    value.map.roi.diff <- matrix(NA, lp.roi, lp.roi)
    statistic.value.map.roi.diff <- matrix(NA, lp.roi, lp.roi)
    #print(coverage.zero.roi.diff)
    for (i in c(1:(lp.roi * (lp.roi + 1)/2 - lp.roi))) {
      print(i)
        kk <- which(roi == coverage.zero.roi.diff[i, 4])
        ll <- which(roi == coverage.zero.roi.diff[i, 5])
        
        value.map.roi.diff[kk, ll] <- as.numeric(coverage.zero.roi.diff[i, 2])
        value.map.roi.diff[ll, kk] <- as.numeric(coverage.zero.roi.diff[i, 2])
        statistic.value.map.roi.diff[kk, ll] <- as.numeric(coverage.zero.roi.diff[i, 3])
        statistic.value.map.roi.diff[ll, kk] <- as.numeric(coverage.zero.roi.diff[i, 3])
    }
    
    diag(value.map.roi.diff) <- -0.01
    diag(statistic.value.map.roi.diff) <- 0
    
    
    ############################################################### 
    
    
    diag(value.map.roi.diff) <- 0
    
    
    
    ####################################### cond1 ##############################
    
    
    
    ommm <- 0
    
    coverage.zero.cond1[is.na(coverage.zero.cond1[, 2]), 2] <- 0
    
    
    prop.res.chisq <- matrix(NA, nrow(coverage.zero.cond1), 1)
    prop.res.chisq.pvalue.cond1 <- matrix(NA, nrow(coverage.zero.cond1), 1)
    val.prop.test <- vector()
    
    for (imk in 1:nrow(coverage.zero.cond1)) {
        
        
        val.prop.test[imk] <- round((coverage.zero.cond1[imk, 2]) * 105, 0)
        prop.res <- stats::prop.test(val.prop.test[imk], 105, alternative = "greater", p = cond1.thres)
        prop.res.chisq[imk, 1] <- prop.res$statistic
        prop.res.chisq.pvalue.cond1[imk, 1] <- prop.res$p.value
    }
    
    coverage.nonzero.cond1.chsq <- cbind(coverage.zero.cond1[, 1], (coverage.zero.cond1[, 2]), prop.res.chisq, loc1.cond1, loc2.cond1)
    colnames(coverage.nonzero.cond1.chsq) <- c("comparison", "Significant coverage for Diff", "Chi-Square", "loc1.cond1", "loc2.cond1")
    ######################################################### 
    
    ######################################## 
    loc1_loc2.cond1 <- paste(coverage.nonzero.cond1.chsq[, "loc1.cond1"], "_", coverage.nonzero.cond1.chsq[, "loc2.cond1"], sep = "")
    
    ######################################## 
    
    p.val.cond1.without.cerebellum <- prop.res.chisq.pvalue.cond1[!(loc1_loc2.cond1 %in% Cerebellum.loc.id)]
    
    ################ Choose alpha level alpha <- 0.05
    order.p.values.cond1 <- p.val.cond1.without.cerebellum[order(p.val.cond1.without.cerebellum)]
    
    m <- length(order.p.values.cond1)
    k <- seq(1, m, 1)
    # BH dependent correction https://en.wikipedia.org/wiki/False_discovery_rate
    c_m <- sum(1/c(1:m))
    
    BH.corr.dep <- (k/(m * c_m)) * alpha
    max.BH.dep.cond1 <- max(order.p.values.cond1[order.p.values.cond1 <= BH.corr.dep])
    
    significant.BH.dep <- rep(2, m)
    significant.BH.dep[order.p.values.cond1 <= BH.corr.dep] <- 1
    significant.BH.dep[order.p.values.cond1 > BH.corr.dep] <- 0
    table(significant.BH.dep)
    
    row.num.BH.dep <- which(prop.res.chisq.pvalue.cond1 == max.BH.dep.cond1)
    chisq.BH.dep <- abs(prop.res.chisq[row.num.BH.dep])[1]
    
    
    ################ 
    
    coverage.nonzero.cond1.chsq[as.numeric(prop.res.chisq.pvalue.cond1) > max.BH.dep.cond1, 2] <- 0
    
    ################## 
    
    coverage.zero.cond1.all <- coverage.nonzero.cond1.chsq
    
    
    loc.cond1 <- coverage.zero.cond1.all[(coverage.zero.cond1.all[, 4] %in% roi_char), ]
    coverage.zero.roi.cond1 <- loc.cond1[(loc.cond1[, 5] %in% roi_char), ]
    
    lp.roi <- length(roi)
    value.map.roi.cond1 <- matrix(NA, lp.roi, lp.roi)
    statistic.value.map.roi.cond1 <- matrix(NA, lp.roi, lp.roi)
    
    for (i in c(1:(lp.roi * (lp.roi + 1)/2 - lp.roi))) {
        
        kk <- which(roi == coverage.zero.roi.cond1[i, 4])
        ll <- which(roi == coverage.zero.roi.cond1[i, 5])
        
        value.map.roi.cond1[kk, ll] <- as.numeric(coverage.zero.roi.cond1[i, 2])
        value.map.roi.cond1[ll, kk] <- as.numeric(coverage.zero.roi.cond1[i, 2])
        statistic.value.map.roi.cond1[kk, ll] <- as.numeric(coverage.zero.roi.cond1[i, 3])
        statistic.value.map.roi.cond1[ll, kk] <- as.numeric(coverage.zero.roi.cond1[i, 3])
    }
    
    diag(value.map.roi.cond1) <- -0.01
    diag(statistic.value.map.roi.cond1) <- 0
    
    together.cond1.diff <- matrix(NA, lp.roi, lp.roi)
    
    for (i in c(1:(lp.roi - 1))) {
        for (j in c((i + 1):(lp.roi))) {
            
            together.cond1.diff[i, j] <- as.numeric(value.map.roi.cond1[i, j])
            together.cond1.diff[j, i] <- as.numeric(value.map.roi.diff[j, i])
            
        }
    }
    
    diag(together.cond1.diff) <- 0
    
    ########################################################################################### cond2 ############################################################
    
    ommm <- 0
    
    coverage.zero.cond2[is.na(coverage.zero.cond2[, 2]), 2] <- 0
    
    prop.res.chisq <- matrix(NA, nrow(coverage.zero.cond2), 1)
    prop.res.chisq.pvalue.cond2 <- matrix(NA, nrow(coverage.zero.cond2), 1)
    val.prop.test <- vector()
    
    for (imk in 1:nrow(coverage.zero.cond2)) {
        
        # imk <- 38502
        
        val.prop.test[imk] <- round((coverage.zero.cond2[imk, 2]) * 105, 0)
        prop.res <- stats::prop.test(val.prop.test[imk], 105, alternative = "greater", p = cond2.thres)
        prop.res.chisq[imk, 1] <- prop.res$statistic
        prop.res.chisq.pvalue.cond2[imk, 1] <- prop.res$p.value
    }
    
    coverage.nonzero.cond2.chsq <- cbind(coverage.zero.cond2[, 1], (coverage.zero.cond2[, 2]), prop.res.chisq, loc1.cond2, loc2.cond2)
    colnames(coverage.nonzero.cond2.chsq) <- c("comparison", "Significant coverage for Diff", "Chi-Square", "loc1.cond2", "loc2.cond2")
    
    ######################################################### 
    loc1_loc2.cond2 <- paste(coverage.nonzero.cond2.chsq[, "loc1.cond2"], "_", coverage.nonzero.cond2.chsq[, "loc2.cond2"], sep = "")
    
    ######################################## 
    
    p.val.cond2.without.cerebellum <- prop.res.chisq.pvalue.cond2[!(loc1_loc2.cond2 %in% Cerebellum.loc.id)]
    
    ################ Choose alpha level
    alpha <- 0.05
    order.p.values.cond2 <- p.val.cond2.without.cerebellum[order(p.val.cond2.without.cerebellum)]
    
    m <- length(order.p.values.cond2)
    k <- seq(1, m, 1)
    # BH dependent correction https://en.wikipedia.org/wiki/False_discovery_rate
    c_m <- sum(1/c(1:m))
    
    BH.corr.dep <- (k/(m * c_m)) * alpha
    max.BH.dep.cond2 <- max(order.p.values.cond2[order.p.values.cond2 <= BH.corr.dep])
    
    significant.BH.dep <- rep(2, m)
    significant.BH.dep[order.p.values.cond2 <= BH.corr.dep] <- 1
    significant.BH.dep[order.p.values.cond2 > BH.corr.dep] <- 0
    table(significant.BH.dep)
    
    row.num.BH.dep <- which(prop.res.chisq.pvalue.cond2 == max.BH.dep.cond2)
    chisq.BH.dep <- abs(prop.res.chisq[row.num.BH.dep])[1]
    
    
    ################ 
    coverage.nonzero.cond2.chsq[as.numeric(prop.res.chisq.pvalue.cond2) > max.BH.dep.cond2, 2] <- 0
    
    ################## 
    
    coverage.zero.cond2.all <- coverage.nonzero.cond2.chsq
    
    
    loc.cond2 <- coverage.zero.cond2.all[(coverage.zero.cond2.all[, 4] %in% roi_char), ]
    coverage.zero.roi.cond2 <- loc.cond2[(loc.cond2[, 5] %in% roi_char), ]
    
    lp.roi <- length(roi)
    value.map.roi.cond2 <- matrix(NA, lp.roi, lp.roi)
    statistic.value.map.roi.cond2 <- matrix(NA, lp.roi, lp.roi)
    
    for (i in c(1:(lp.roi * (lp.roi + 1)/2 - lp.roi))) {
        
        kk <- which(roi == coverage.zero.roi.cond2[i, 4])
        ll <- which(roi == coverage.zero.roi.cond2[i, 5])
        
        value.map.roi.cond2[kk, ll] <- as.numeric(coverage.zero.roi.cond2[i, 2])
        value.map.roi.cond2[ll, kk] <- as.numeric(coverage.zero.roi.cond2[i, 2])
        statistic.value.map.roi.cond2[kk, ll] <- as.numeric(coverage.zero.roi.cond2[i, 3])
        statistic.value.map.roi.cond2[ll, kk] <- as.numeric(coverage.zero.roi.cond2[i, 3])
    }
    
    diag(value.map.roi.cond2) <- -0.01
    diag(statistic.value.map.roi.cond2) <- 0
    
    together.cond1.cond2 <- matrix(NA, lp.roi, lp.roi)
    
    for (i in c(1:(lp.roi - 1))) {
        for (j in c((i + 1):(lp.roi))) {
            
            together.cond1.cond2[i, j] <- as.numeric(value.map.roi.cond1[i, j])
            together.cond1.cond2[j, i] <- as.numeric(value.map.roi.cond2[j, i])
            
        }
    }
    
    diag(together.cond1.cond2) <- -0.01
    
    matrixList <- list()
    
    matrixList$diff <- value.map.roi.diff
    matrixList$cond1 <- value.map.roi.cond1
    matrixList$cond2 <- value.map.roi.cond2
    matrixList$diff_cond1 <- together.cond1.diff
    matrixList$cond1_cond2 <- together.cond1.cond2
    
    return(matrixList)
    
}
