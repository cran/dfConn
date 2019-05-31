# Summarise table Summarize results (coverage) of confidence interval data param data.CI.path path to the param comparisonsName A character string, the indicator of group
# comparisons, typically can be extracted from the filename.  return A list of matrix family summarize_table

table_summarize.coverage <- function(dataList.CI, comparisonsName) {
    
    coverage.zero.cond.diff <- data.frame(ROI_COMPARISON = comparisonsName, ZERO_COVERAGE = rep(0, length(comparisonsName)))
    coverage.zero.cond.1 <- data.frame(ROI_COMPARISON = comparisonsName, ZERO_COVERAGE = rep(0, length(comparisonsName)))
    coverage.zero.cond.2 <- data.frame(ROI_COMPARISON = comparisonsName, ZERO_COVERAGE = rep(0, length(comparisonsName)))
    
    # Loop over comparisons and fill the table with average coverage
    for (i in 1:length(comparisonsName)) {
        
        data.p <- dataList.CI[[i]]
        ntp <- ncol(data.p)
        coverage.zero.cond.diff[i, 2] <- sum((data.p[4, ] < 0) & (0 < data.p[5, ]))/ntp
        
        coverage.zero.cond.1[i, 2] <- sum((data.p[6, ] < 0) & (0 < data.p[7, ]))/ntp
        
        coverage.zero.cond.2[i, 2] <- sum((data.p[8, ] < 0) & (0 < data.p[9, ]))/ntp
        
        
    }
    
    coverage.below.zero.cond.2 <- data.frame(ROI_COMPARISON = comparisonsName,
                                             ZERO_COVERAGE = rep(0, length(comparisonsName)))
    
    coverage.below.zero.cond.diff <- data.frame(ROI_COMPARISON = comparisonsName,
                                                ZERO_COVERAGE = rep(0, length(comparisonsName)))
    
    coverage.below.zero.cond.1 <- data.frame(ROI_COMPARISON = comparisonsName,
                                             ZERO_COVERAGE = rep(0, length(comparisonsName)))
    
    coverage.above.zero.cond.2 <- data.frame(ROI_COMPARISON = comparisonsName,
                                             ZERO_COVERAGE = rep(0, length(comparisonsName)))
    
    coverage.above.zero.cond.diff <- data.frame(ROI_COMPARISON = comparisonsName,
                                                ZERO_COVERAGE = rep(0, length(comparisonsName)))
    
    coverage.above.zero.cond.1 <- data.frame(ROI_COMPARISON = comparisonsName,
                                             ZERO_COVERAGE = rep(0, length(comparisonsName)))
    
    for (i in 1:length(comparisonsName)) {
        
        data.p <- dataList.CI[[i]]
        ntp <- ncol(data.p)
        
        coverage.below.zero.cond.diff[i, 2] <- sum((data.p[4, ] < 0) & (data.p[5, ] < 0))/ntp
        
        coverage.below.zero.cond.1[i, 2] <- sum((data.p[6, ] < 0) & (data.p[7, ] < 0))/ntp
        
        coverage.below.zero.cond.2[i, 2] <- sum((data.p[8, ] < 0) & (data.p[9, ] < 0))/ntp
        
        coverage.above.zero.cond.diff[i, 2] <- sum((data.p[4, ] > 0) & (data.p[5, ] > 0))/ntp
        
        coverage.above.zero.cond.1[i, 2] <- sum((data.p[6, ] > 0) & (data.p[7, ] > 0))/ntp
        
        coverage.above.zero.cond.2[i, 2] <- sum((data.p[8, ] > 0) & (data.p[9, ] > 0))/ntp
    }
    
    # return a list containg three tables, coverage.zero.cond.diff, coverage.zero.cond.1, coverage.zero.cond.2, coverage.below.zero.cond.2, coverage.below.zero.cond.1,
    # coverage.below.zero.cond.diff, coverage.above.zero.cond.diff, coverage.above.zero.cond.1, coverage.above.zero.cond.2
    
    colnames(coverage.zero.cond.diff) <- c("ROI_COMPARISON", "ZERO-COVERAGE(condition_diff)")
    colnames(coverage.zero.cond.1) <- c("ROI_COMPARISON", "ZERO-COVERAGE(condition1)")
    colnames(coverage.zero.cond.2) <- c("ROI_COMPARISON", "ZERO-COVERAGE(condition2)")
    
    
    colnames(coverage.below.zero.cond.2) <- c("ROI_COMPARISON", "BELOW-ZERO-COVERAGE(condition2)")
    colnames(coverage.below.zero.cond.1) <- c("ROI_COMPARISON", "BELOW-ZERO-COVERAGE(condition1)")
    colnames(coverage.below.zero.cond.diff) <- c("ROI_COMPARISON", "BELOW-ZERO-COVERAGE(condition_diff)")
    
    colnames(coverage.above.zero.cond.diff) <- c("ROI_COMPARISON", "ABOVE-ZERO-COVERAGE(condition_diff)")
    colnames(coverage.above.zero.cond.1) <- c("ROI_COMPARISON", "ABOVE-ZERO-COVERAGE(condition1)")
    colnames(coverage.above.zero.cond.2) <- c("ROI_COMPARISON", "ABOVE-ZERO-COVERAGE(condition2)")
    
    list(coverage.zero.cond.diff = coverage.zero.cond.diff, coverage.zero.cond.1 = coverage.zero.cond.1, coverage.zero.cond.2 = coverage.zero.cond.2, coverage.below.zero.cond.2 = coverage.below.zero.cond.2, 
        coverage.below.zero.cond.diff = coverage.below.zero.cond.diff, coverage.below.zero.cond.1 = coverage.below.zero.cond.1, coverage.above.zero.cond.diff = coverage.above.zero.cond.diff, 
        coverage.above.zero.cond.1 = coverage.above.zero.cond.1, coverage.above.zero.cond.2 = coverage.above.zero.cond.2)
    
}
