#include <Rcpp.h>
using namespace Rcpp;
//' First zero index for each column of a matrix
//'
//' For each column, returns the 1-based index of the first zero value.
//' If a column contains no zeros, returns NA for that column.
//'
//' @param mat NumericMatrix
//' @return IntegerVector of length ncol(mat) with 1-based indices or NA.
//' @examples
//' m <- matrix(c(1, 0, 2,0, 3, 1),nrow = 3, ncol = 2)
//' firstZeroIndex(m)
//' @export
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

//' Generate a matrix of binomial random values
//'
//' Create a matrix of size nrows by ncols where each column j uses its own
//' success probability probs j.
//'
//' @param nrows Integer, number of rows.
//' @param ncols Integer, number of columns.
//' @param probs Numeric vector of length ncols, column-wise success probabilities in (0, 1).
//' @return IntegerMatrix of dimension nrows by ncols.
//' @examples
//' set.seed(1)
//' rbinom_matrix_vec(5, 3, rep(0.3, 3))
//' @export
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

