# Correlation Scale This function gives the scaled-correlation param x An numeric vector keywords internal return A numeric vector with the correlation scaled


corr.scale <- function(x) {
    (exp(2 * x) - 1)/(exp(2 * x) + 1)
}




# For creation of O'Sullivan-type Z matrices.  title O'Sullivan-type Z matrices.  description For creation of O'Sullivan-type Z matrices.  param x is a function argument
# param range.x defines the range of an argument x param intKnots are the internal knots param drv is the derivative keywords internal return Z matrix references Kindle,
# E.  'Semiparametric Regression (Cambridge Series in Statistical and Probabilistic Mathematics).' (2003).  importFrom stats prop.test qnorm quantile

ZOSull <- function(x, range.x, intKnots, drv = 0) {
    if (drv > 2) 
        stop("splines not smooth enough for more than 2 derivatives")
    
    requireNamespace("splines", "stats")
    
    # Set defaults for `range.x' and `intKnots'
    
    if (missing(range.x)) 
        range.x <- c(1.05 * min(x) - 0.05 * max(x), 1.05 * max(x) - 0.05 * min(x))
    
    if (missing(intKnots)) {
        numIntKnots <- min(length(unique(x)), 35)
        intKnots <- stats::quantile(unique(x), seq(0, 1, length = (numIntKnots + 2))[-c(1, (numIntKnots + 2))])
    }
    numIntKnots <- length(intKnots)
    
    # Obtain the penalty matrix.
    
    allKnots <- c(rep(range.x[1], 4), intKnots, rep(range.x[2], 4))
    K <- length(intKnots)
    L <- 3 * (K + 8)
    xtilde <- (rep(allKnots, each = 3)[-c(1, (L - 1), L)] + rep(allKnots, each = 3)[-c(1, 2, L)])/2
    wts <- rep(diff(allKnots), each = 3) * rep(c(1, 4, 1)/6, K + 7)
    Bdd <- splines::spline.des(allKnots, xtilde, derivs = rep(2, length(xtilde)), outer.ok = TRUE)$design
    Omega <- t(Bdd * wts) %*% Bdd
    
    # Use the spectral decomposition of Omega to obtain Z.
    
    eigOmega <- eigen(Omega)
    indsZ <- 1:(numIntKnots + 2)
    UZ <- eigOmega$vectors[, indsZ]
    LZ <- t(t(UZ)/sqrt(eigOmega$values[indsZ]))
    
    # Perform stability check.
    
    ind <- (numIntKnots + 3):(numIntKnots + 4)
    UX <- eigOmega$vectors[, ind]
    L <- cbind(UX, LZ)
    stabCheck <- t(crossprod(L, t(crossprod(L, Omega))))
    if (sum(stabCheck^2) > 1.0001 * (numIntKnots + 2)) 
        print("WARNING: NUMERICAL INSTABILITY ARISING\\
          FROM SPECTRAL DECOMPOSITION")
    
    # Obtain B and post-multiply by LZ matrix to get Z.
    
    B <- splines::spline.des(allKnots, x, derivs = rep(drv, length(x)), outer.ok = TRUE)$design
    
    Z <- B %*% LZ
    
    # Add the `range.x' and 'intKnots' as attributes of the return object.
    
    attr(Z, "range.x") <- range.x
    attr(Z, "intKnots") <- intKnots
    
    # Return Z matrix with 2 attributes.
    
    return(Z)
}



# Name parsing for matrix plots.  title Matrix plots label names parsing.  description For generating correct axis label for matrix plot.  param roiNames list of roi
# names to plot param delay need to label as delay or not?  keywords internal return Label vectors

matplot_roisnames_setup <- function(roiNames, delay = FALSE) {
    name2plots <- c()
    
    if (!delay) {
        for (i in 1:length(roiNames)) {
            name2plots <- c(name2plots, " ", roiNames[i])
        }
    } else {
        for (i in 1:length(roiNames)) {
            name2plots <- c(name2plots, " ", paste(roiNames[i], "(D)"))
        }
    }
    
    
    
    return(c(name2plots, " "))
}




# Extract info from file names

# description Extract rois from a list of filenames, and return a dataframe containing two columns, which contain the pairwise rois.  title Extract rois from a list of
# filenames param path The parent path of the files.  param file_fmt the format of the files, default is 'csv'.  param prefix Prefix of the file names before roi, eg.,
# 'lme_result_12_23_x_delay_1.csv', where prefix = 'lme_result_', default is empty string.  param roi_interest optional, specify a subset of rois

filename_extract_rois <- function(path, roi_interest, file_fmt = "csv", prefix = "") {
    
    
    filenames <- list.files(path, pattern = paste(prefix, ".*.", file_fmt, sep = ""))
    num_files <- length(filenames)
    
    # print(filenames)
    filenames_no_Suffix <- unlist(strsplit(filenames, paste(".", file_fmt, sep = "")))  # Ged rid of suffix
    
    filenames_main <- unlist(strsplit(filenames_no_Suffix, prefix))[seq(2, (length(filenames) * 2), 2)]  # Ged rid of prefix
    
    
    ### extract roi name from file names
    k <- length(unlist(strsplit(filenames_main[1], "_")))
    
    roi1.all <- unlist(strsplit(filenames_main, "_"))[seq(1, (num_files * k), k)]
    roi2.all <- unlist(strsplit(filenames_main, "_"))[seq(2, (num_files * k), k)]
    
    
    if (!is.character(roi_interest)) {
        roi1.int <- which(as.integer(roi1.all) %in% roi_interest)
        roi2.int <- which(as.integer(roi2.all) %in% roi_interest)
        
        selected_comp <- filenames_main[roi1.int[which(roi1.int %in% roi2.int)]]  # select those filenames with both roi1 and roi2
        roi1 <- roi1.all[roi1.int[which(roi1.int %in% roi2.int)]]
        roi2 <- roi2.all[roi2.int[which(roi2.int %in% roi1.int)]]
        
        results <- list(roi1 = roi1, roi2 = roi2, main_filenames = filenames_main, selected_comp = selected_comp, all_filenames = filenames)
    } else {
        roi1.int <- as.integer(roi1.all)
        roi2.int <- as.integer(roi2.all)
        
        roi1 <- roi1.all[roi1.int[which(roi1.int %in% roi2.int)]]
        roi2 <- roi2.all[roi2.int[which(roi2.int %in% roi1.int)]]
        results <- list(roi1 = roi1, roi2 = roi2, main_filenames = filenames_main, selected_comp = filenames_main, all_filenames = filenames)
    }
    
    
    results
    
}

# Unit function for extracting roi info from string with 'roi1_roi2_xxxxxx' format


# description Unit function for extracting roi info from string with 'roi1_roi2_xxxxxx' format title Extract rois from a list of filenames param strings A list of string
# with format like: 'xxx_roi1_roi2_xxxxxx' format

extract_roi_pair <- function(strings) {
    requireNamespace("stringr")
    roi_pair <- stringr::str_extract(strings, pattern = "[0-9]+_[0-9]+")
    roi1 <- unlist(lapply(strsplit(roi_pair, "_"), function(x) {
        x[1]
    }))
    roi2 <- unlist(lapply(strsplit(roi_pair, "_"), function(x) {
        x[2]
    }))
    return(list(roi1 = roi1, roi2 = roi2))
}


# Combine function for foreach()
comb <- function(x, ...) {
    lapply(seq_along(x), function(i) c(x[[i]], lapply(list(...), function(y) y[[i]])))
}



