#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
IntegerVector firstZeroIndex(NumericMatrix mat) {
  int nrow = mat.nrow(), ncol = mat.ncol();
  IntegerVector result(ncol);
  for (int j = 0; j < ncol; j++) {
    result[j] = NA_INTEGER; // Initialize with NA
    for (int i = 0; i < nrow; i++) {
      if (mat(i, j) == 0) {
        result[j] = i + 1;
        break;
      }
    }
  }
  return result;
}


// [[Rcpp::export]]
IntegerMatrix rbinom_matrix_vec(int nrows, int ncols, NumericVector probs) {
  RNGScope scope;
  
  IntegerMatrix mat(nrows, ncols);
  
  for (int j = 0; j < ncols; j++) {
    for (int i = 0; i < nrows; i++) {
      mat(i, j) = R::rbinom(1, probs[i]);
    }
  }
  
  return mat;
}

