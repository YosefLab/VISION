// Utilities for fast matrix operations to be wrapped with R in FastProjectR

#include <Rcpp.h> 
using namespace Rcpp;

// [[Rcpp::export]]
void point_mult(NumericMatrix & X, NumericVector & Y) {
	// In place pointwise multiplication of matrices X and Y, both of which are N x M
	// NOTE: mutates the X matrix in memory

	unsigned int ncol = X.ncol();
	unsigned int nrow = X.nrow();
	for (unsigned int j = 0; j < ncol*nrow; j++) {
		X[j] *= Y[j];
	}
}

// [[Rcpp::export]]
NumericVector multMat(NumericMatrix X, NumericMatrix Y) {
	NumericVector res = X * Y;
	res.attr("dim") = Dimension(X.nrow(), X.ncol());
	return res;

}

// [[Rcpp::export]]
NumericMatrix calcSigScore(NumericMatrix & X, NumericVector & y, NumericMatrix & Z) {
	int ncol = X.ncol();
	int nrow = X.nrow();
	int counter = 0;
	for (int j=0; j < ncol; j++) {
		for (int i = 0; i < nrow; i++) {
				X[counter++] *= y[i];
		}
	}

	NumericMatrix res = (NumericMatrix) multMat(X, Z);
	return res;

}

//RCPP_MODULE(utilities_mod) {
//	function( "pointwise_mult", &point_mult );
//}
