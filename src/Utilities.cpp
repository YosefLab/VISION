// Utilities for fast matrix operations to be wrapped with R in FastProjectR

#include <Rcpp.h> 
#include <iostream>
#include <string>
#include <sstream>
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
	// Matrix multiplication where new matrix is created to store result 
	NumericVector res = X * Y;
	res.attr("dim") = Dimension(X.nrow(), X.ncol());
	return res;

}
