// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"

using namespace arma;
using namespace Rcpp;
// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
//
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
arma::uvec myseq(int start, int end, int by){
  int s = start;
  int n = 0;

  while (s <= end){
    n += 1;
    s += by;
  }

  arma::uvec res = uvec(n);
  s = start;

  for (int i = 0; i<n; i++){
    res(i) = s;
    s += by;
  }

  return res;

}

// [[Rcpp::export]]
SEXP boot_func(int Boot, arma::umat sw_mat, arma::mat Z_mat, arma::mat L){
  int mBoot = Boot;
  arma::mat mZ_mat = Z_mat;
  arma::umat mSw_mat = sw_mat-1;
  arma::mat mL = L;
  int n1 = mZ_mat.n_rows;	int n2 = mZ_mat.n_cols;
  arma::mat X_stern(n1,n2);
  arma::mat ztemp(n1,n2);
  arma::mat X_stern_all(n1*mBoot,n2);
  for (int b = 0; b < mBoot; b++){
    arma::mat temp = mZ_mat.cols(mSw_mat.col(b));
    ztemp = temp;  // temp will go out of scope at the end of the loop
    arma::vec Z_vec_stern = zeros(n1*n2);
    for (int col_i = 0; col_i < n2; col_i++) {
      int from = col_i*n1;
      int to = (col_i+1)*n1 - 1;
      Z_vec_stern.subvec(from, to)= ztemp.col(col_i);
    }
    arma::colvec X_vec_stern = mL*Z_vec_stern;
    for (int j = 0; j < n2; j++){
      int from = j*n1;
      int to = (j+1)*n1 - 1;
      X_stern.col(j) = X_vec_stern.subvec(from, to);
    }
    for (int i = 0; i < n1; i++) {
      X_stern_all.row(b*n1+i) = X_stern.row(i);
    }
  }
  return(wrap(X_stern_all));
}

// [[Rcpp::export]]
Rcpp::NumericVector mm_func(arma::mat Xjj, Function F, 
                            int mm,
                            int n_boot,
                            int window_size = 20,
                            int n_sig = 750
){

  arma::ivec temp_v = ivec(2);
  temp_v(0) = mm-1+window_size;
  temp_v(1) = n_sig;

  arma::mat X = Xjj.cols(mm-1, as_scalar(arma::min(temp_v))-1);


  Rcpp::NumericVector res = F(X, n_boot);

  return res;

}

// [[Rcpp::export]]
arma::mat boot_cjj_func(arma::mat Xjj,
                        Function MLPB3,
                  int boot_rep = 250,
                  int n_boot = 1,
                  int n_sig = 750,
                  int window_size = 20){


  arma::mat bootjj(2, n_sig);
  uvec s = myseq(1, n_sig, window_size);

  arma::rowvec x = rowvec(2*n_sig);
  int s_n = s.n_elem;

  int marker_ind = 0;
  for(int i = 1; i <= s_n; i++){

    int mm = s(i-1);

    NumericVector x_temp_num = mm_func(Xjj, MLPB3, mm, n_boot, window_size, n_sig);

    rowvec x_temp = arma::rowvec(x_temp_num);
    x(span(marker_ind, marker_ind + x_temp.n_elem - 1)) = x_temp;
    marker_ind += x_temp.n_elem;

  }

  mat xm = x;
  xm.reshape(2, n_sig);
  return xm;

}

// [[Rcpp::export]]
arma::rowvec a_func(arma::mat bootjj,
              int up_limit = 731,
              int n_sig = 750,
              int window_size = 20){

  arma::rowvec a(up_limit);

  for (int i = 1; i <= up_limit; i++){

    arma::rowvec row1 = bootjj.row(0);
    arma::rowvec row2 = bootjj.row(1);


    int k = i;
    int up;


    if (k-1+window_size >= n_sig){
      up = n_sig;
    }else{
      up = k-1+window_size;
    }


    //print(wrap(up));

    rowvec subrow1 = row1(span(k-1, up-1));
    rowvec subrow2 = row2(span(k-1, up-1));
    a(k-1) = as_scalar(cor(subrow1, subrow2));

  }

  return (a);

}

// [[Rcpp::export]]
arma::mat bootcjj(arma::mat X, Function MLPB3,
                  int up_limit = 731,
                  int window_size = 20,
                  int boot_rep = 250,
                  int n_boot = 1,
                  int n_sig = 750){

  mat boot_cjj(boot_rep, up_limit);

  for (int i = 1; i<=boot_rep; i++){

    mat bootjj = (boot_cjj_func(X, MLPB3, boot_rep, n_boot,n_sig,window_size));
    rowvec a = a_func(bootjj,up_limit, n_sig, window_size);
    boot_cjj.row(i-1) = a;
  }

  return (boot_cjj);
}


