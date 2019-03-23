#' @title Multivariate Linear Bootstrapping - core function
#' @param X data matrix
#' @param Boot integer, number of bootstrap samples to be generated
#' @param l numeric, the banding parameter, default is 1
#' @param eps numeric the parameters for making Gamma_kappa_matrix positive definite if necessary, default is 1;
#' @param beta numeric the parameters for making Gamma_kappa_matrix positive definite if necessary, default is 1;
#' @param l_automatic numeric the banding parameter, default is 1, data-adaptively;
#' @param l_automatic_local numeric, the banding parameter default is 0 
#' @return  a numeric matrix with n*boot rows and m columns, where n, m refer to the number of rows and number of columns in input data matrix, and boot refer to number of bootstrap samples to be generated.
#' @examples 
#' # Multivariate linear bootstrapping on a random matrix with 2 rows and 4 columns
#' MLPB3(matrix(rnorm(16),2,4), 3)
#' 
#' @export
#' @useDynLib dfConn


MLPB3 <- function(X, Boot, l = 1, eps = 1, beta = 1, l_automatic = 0, l_automatic_local = 0) {
    requireNamespace("stats")
    # Dimension of observed values
    n1 <- dim(X)[1]
    
    # sample sizes of (n1-dimensional) data observed
    n2 <- dim(X)[2]
    N <- n1 * n2
    # Boot<-100
    l_matrix <- matrix(rep(l, n1^2), n1, n1)
    
    # Defining function kappa (trapezoid [cf. McMurry & Politis (2010)])
    kappa <- function(x, l) {
        kappa <- matrix(0, dim(l)[1], dim(l)[2])
        kappa[l > 0] <- 2 - abs(x/l[l > 0])
        kappa[abs(x/l) <= 1 & l != 0] <- 1
        kappa[abs(x/l) > 2 & l != 0] <- 0
        kappa
    }
    
    # Rearranging the data columnwise in one large n1*n2-dimensional vector
    X_vec <- X
    dim(X_vec) <- NULL
    
    ######################################################################## Here starts the computation of covariances ###
    
    # Estimating all multivariate covariances and putting them in one large n1*(2n2-1) matrix C_matrix<-matrix(rep(0,(n1*n1*(2*n2-1))),n1,(2*n2-1)*n1)
    # C_matrix_kappa_trans<-matrix(rep(0,(n1*n1*(2*n2-1))),n1,(2*n2-1)*n1) R_matrix_kappa_trans<-matrix(rep(0,(n1*n1*(2*n2-1))),n1,(2*n2-1)*n1)
    C_matrix <- C_matrix_kappa_trans <- R_matrix_kappa_trans <- matrix(0, n1, (2 * n2 - 1) * n1)
    
    # Sample mean
    if (n1 == 1) 
        X_mean <- mean(X) else X_mean <- rowMeans(X)
    # if(n1!=1){ X_mean<-rep(0,n1) for(k in 1:n1){ X_mean[k]<-mean(X[k,]) } }
    
    # Computing C_matrix
    for (k in (-(n2 - 1)):(n2 - 1)) {
        blub <- 0
        if (n1 == 1) 
            C_matrix[, ((((n2 + k) - 1) * n1 + 1):(((n2 + k) - 1) * n1 + n1))] <- (1/n2) * t(X[(max(1, 1 - k) + k):(min(n2, n2 - k) + k)] - X_mean) %*% (X[(max(1, 1 - k)):(min(n2, 
                n2 - k))] - X_mean) else C_matrix[, ((((n2 + k) - 1) * n1 + 1):(((n2 + k) - 1) * n1 + n1))] <- (1/n2) * (X[, (max(1, 1 - k) + k):(min(n2, n2 - k) + k)] - X_mean) %*% t(X[, (max(1, 1 - 
            k)):(min(n2, n2 - k))] - X_mean)
        # if(n1!=1) C_matrix[,((((n2+k)-1)*n1+1):(((n2+k)-1)*n1+n1))]<-(1/n2)*(X[,(max(1,1-k)+k):(min(n2,n2-k)+k)]-X_mean)%*%t(X[,(max(1,1-k)):(min(n2,n2-k))]-X_mean)
    }
    
    if (l_automatic == 1) {
        
        ################################## Computing Autocorrelations ###
        
        R_matrix <- matrix(0, n1, (2 * n2 - 1) * n1)
        denominator_matrix <- matrix(0, n1, n1)
        diag(denominator_matrix) <- sqrt(diag(C_matrix[, (n1 * (n2 - 1) + 1):(n1 * n2)]))
        denominator_matrix <- solve(denominator_matrix)
        for (k in 1:(2 * n2 - 1)) {
            R_matrix[, ((k - 1) * n1 + 1):(k * n1)] <- denominator_matrix %*% C_matrix[, ((k - 1) * n1 + 1):(k * n1)] %*% denominator_matrix
        }
        
        M_zero <- 2
        M_zero_sqrt <- M_zero * sqrt(log10(n2)/n2)
        
        if (l_automatic_local == 0) {
            h <- 1
            marker <- 0
            while (marker == 0) {
                # Kn_vec<-seq(1,max(5,log10(n2))) if(max((abs(R_matrix[,((2*n2-1)+h*n1):((2*n2-1+n1-1)+n1*(max(5,log10(n2))+h-1))])))<M_zero*sqrt(log10(n2)/n2)) marker<-1
                if (max((abs(R_matrix[, (n1 * (n2 - 1) + 1 + h * n1):(n1 * n2 + n1 * (max(5, log10(n2)) + h - 1))]))) < M_zero_sqrt) 
                  marker <- 1
                l_adapt <- h - 1
                h <- h + 1
            }
            l <- l_adapt
            l_matrix_adapt <- matrix(rep(l, n1^2), n1, n1)
            l_matrix <- l_matrix_adapt
        }
        
        if (l_automatic_local == 1) {
            l_matrix_adapt <- matrix(0, n1, n1)
            marker <- matrix(0, n1, n1)
            for (j in 1:n1) {
                for (k in 1:n1) {
                  h <- 1
                  while (marker[j, k] == 0) {
                    if (max(abs(R_matrix[j, seq((n1 * (n2 - 1) + 1 + h * n1 + k - 1), (n1 * n2 + n1 * (max(5, log10(n2)) + h - 1) + k - 1), by = n1)])) < M_zero_sqrt) 
                      marker[j, k] <- 1
                    l_matrix_adapt[j, k] <- h - 1
                    h <- h + 1
                  }
                }
            }
            l_matrix <- l_matrix_adapt
        }
        
    }
    
    # Computing C_matrix_kappa_trans
    for (k in (-(n2 - 1)):(n2 - 1)) {
        blub <- 0
        if (n1 == 1) 
            C_matrix_kappa_trans[, ((((n2 + k) - 1) * n1 + 1):(((n2 + k) - 1) * n1 + n1))] <- (kappa(k, l_matrix)/n2) * t(X[(max(1, 1 - k) + k):(min(n2, n2 - k) + k)] - X_mean) %*% 
                (X[(max(1, 1 - k)):(min(n2, n2 - k))] - X_mean)
        if (n1 != 1 && k < 0) 
            C_matrix_kappa_trans[, ((((n2 + k) - 1) * n1 + 1):(((n2 + k) - 1) * n1 + n1))] <- t(t((kappa(k, l_matrix)/n2)) * (X[, (max(1, 1 - k) + k):(min(n2, n2 - k) + k)] - 
                X_mean) %*% t(X[, (max(1, 1 - k)):(min(n2, n2 - k))] - X_mean))
        if (n1 != 1 && k == 0) 
            C_matrix_kappa_trans[, ((((n2 + k) - 1) * n1 + 1):(((n2 + k) - 1) * n1 + n1))] <- t((1/n2) * (X[, (max(1, 1 - k) + k):(min(n2, n2 - k) + k)] - X_mean) %*% t(X[, 
                (max(1, 1 - k)):(min(n2, n2 - k))] - X_mean))
        if (n1 != 1 && k > 0) 
            C_matrix_kappa_trans[, ((((n2 + k) - 1) * n1 + 1):(((n2 + k) - 1) * n1 + n1))] <- t((kappa(k, l_matrix)/n2) * (X[, (max(1, 1 - k) + k):(min(n2, n2 - k) + k)] - 
                X_mean) %*% t(X[, (max(1, 1 - k)):(min(n2, n2 - k))] - X_mean))
    }
    ############################################################################# Computing the diagonal matrix for standardization to get autocorrelations (as above)
    denominator_matrix <- matrix(0, n1, n1)
    diag(denominator_matrix) <- sqrt(diag(C_matrix[, (n1 * (n2 - 1) + 1):(n1 * n2)]))
    denominator_matrix <- solve(denominator_matrix)
    
    # Computing R_matrix_kappa_trans (for equivariance to construct positive definite estimator)
    for (k in (-(n2 - 1)):(n2 - 1)) {
        blub <- 0
        if (n1 == 1) 
            R_matrix_kappa_trans[, ((((n2 + k) - 1) * n1 + 1):(((n2 + k) - 1) * n1 + n1))] <- (1/stats::var(as.vector(X))) * (kappa(k, l_matrix)/n2) * t(X[(max(1, 1 - k) + 
                k):(min(n2, n2 - k) + k)] - X_mean) %*% (X[(max(1, 1 - k)):(min(n2, n2 - k))] - X_mean)
        if (n1 != 1 && k < 0) 
            R_matrix_kappa_trans[, ((((n2 + k) - 1) * n1 + 1):(((n2 + k) - 1) * n1 + n1))] <- denominator_matrix %*% t(t((kappa(k, l_matrix)/n2)) * (X[, (max(1, 1 - k) + 
                k):(min(n2, n2 - k) + k)] - X_mean) %*% t(X[, (max(1, 1 - k)):(min(n2, n2 - k))] - X_mean)) %*% denominator_matrix
        if (n1 != 1 && k == 0) 
            R_matrix_kappa_trans[, ((((n2 + k) - 1) * n1 + 1):(((n2 + k) - 1) * n1 + n1))] <- denominator_matrix %*% t((1/n2) * (X[, (max(1, 1 - k) + k):(min(n2, n2 - k) + 
                k)] - X_mean) %*% t(X[, (max(1, 1 - k)):(min(n2, n2 - k))] - X_mean)) %*% denominator_matrix
        if (n1 != 1 && k > 0) 
            R_matrix_kappa_trans[, ((((n2 + k) - 1) * n1 + 1):(((n2 + k) - 1) * n1 + n1))] <- denominator_matrix %*% t((kappa(k, l_matrix)/n2) * (X[, (max(1, 1 - k) + k):(min(n2, 
                n2 - k) + k)] - X_mean) %*% t(X[, (max(1, 1 - k)):(min(n2, n2 - k))] - X_mean)) %*% denominator_matrix
    }
    
    # Arranging the autocovariance matrices in R_matrix_kappa_trans in one large matrix corresponding to X_vec
    R_kappa_matrix <- matrix(rep(0, N * N), N, N)
    for (i in 1:n2) {
        R_kappa_matrix[(((i - 1) * n1 + 1):((i - 1) * n1 + n1)), ] <- R_matrix_kappa_trans[, (n1 * (n2 - i) + 1):(n1 * (2 * n2 - i))]
    }
    
    R_kappa_matrix_hilf <- R_kappa_matrix
    
    # Computing the eigenvalues and eigenvectors of Gamma_kappa_matrix
    spectraldecomposition <- eigen(R_kappa_matrix)
    
    # orthogonal matrix
    S <- spectraldecomposition$vectors
    
    # diagonal matrix
    D <- diag(spectraldecomposition$values)
    
    D_eps_counter <- 0
    # Manipulating the diagonal matrix D to ensure positive definiteness
    if (min(diag(D)) <= eps * N^(-beta)) {
        D_eps_counter <- D_eps_counter + 1
        D_eps <- D
        D_non_definite <- D
        # for(k in 1:N){ D_eps[k,k]<-max(D[k,k],eps*N^(-beta)) }
        diag(D_eps) <- ifelse(diag(D) > eps * N^(-beta), diag(D), eps * N^(-beta))
        # Computing the banded positive definite covariance matrix estimator
        R_kappa_matrix <- S %*% D_eps %*% t(S)
    }
    
    # Get back to autocovariances!
    Gamma_kappa_matrix <- matrix(rep(0, N * N), N, N)
    nominator_matrix <- solve(denominator_matrix)
    for (i in 1:n2) {
        for (j in 1:n2) {
            Gamma_kappa_matrix[(((i - 1) * n1 + 1):((i - 1) * n1 + n1)), (((j - 1) * n1 + 1):((j - 1) * n1 + n1))] <- nominator_matrix %*% R_kappa_matrix[(((i - 1) * n1 + 
                1):((i - 1) * n1 + n1)), (((j - 1) * n1 + 1):((j - 1) * n1 + n1))] %*% nominator_matrix
        }
    }
    
    
    # Computing the (lower-left) Cholesky decomposition
    L <- t(chol(Gamma_kappa_matrix))
    L_inv <- solve(L)
    
    # Multilpication of L and centered X_vec that is X_vec_cent X_cent<-matrix(rep(0,n1*n2),n1,n2) for(k in 1:n1){ X_cent[k,]<-X[k,]-X_mean[k] }
    X_cent <- X - rowMeans(X)
    
    # Rearranging the data columnwise in one large n1*n2-dimensional vector
    X_vec_cent <- X_cent
    dim(X_vec_cent) <- NULL
    
    # 'Whitening' the centered data
    W_vec <- L_inv %*% X_vec_cent
    
    # Computing the centered and standardized residuals Z_vec
    W_vec_cent <- W_vec - mean(W_vec)
    W_vec_var <- stats::var(W_vec)
    dim(W_vec_var) <- NULL
    Z_vec <- (1/sqrt(W_vec_var)) * W_vec_cent
    
    W_matrix <- matrix(W_vec, n1, n2)
    # W_matrix_cent<-matrix(rep(0,n1*n2),n1,n2) for(k in 1:n1){ W_matrix_cent[k,]<-W_matrix[k,]-mean(W_matrix[k,]) }
    
    W_matrix_cent <- W_matrix - rowMeans(W_matrix)
    
    W_matrix_var <- (1/n2) * W_matrix_cent %*% t(W_matrix_cent)
    W_matrix_var_inv <- solve(t(chol(W_matrix_var)))
    
    # W_matrix_stud<-matrix(rep(0,n1*n2),n1,n2) for(k in 1:n2){ W_matrix_stud[,k]<-W_matrix_var_inv%*%W_matrix_cent[,k] }
    W_matrix_stud <- W_matrix_var_inv %*% W_matrix_cent
    
    Z_vec <- W_matrix_stud
    dim(Z_vec) <- NULL
    
    ########################################################################## Here starts the LPB bootstrap loop ###
    
    
    # rewrite LPB bootstrap loop using Rcpp, Xiaochun Li
    sv_draw <- replicate(Boot, sample(1:n2, n2, replace = TRUE))
    boot_func(Boot, sv_draw, matrix(Z_vec, nrow = n1), L)
}
