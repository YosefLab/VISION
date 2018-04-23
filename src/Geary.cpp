// Utilities for fast matrix operations to be wrapped with R in FastProjectR

#include <Rcpp.h> 
#include <iostream>
#include <string>
using namespace Rcpp;

// [[Rcpp::export]]
double geary(NumericVector X, NumericMatrix W) {
    unsigned int nrow = W.nrow();
    unsigned int ncol = W.ncol();

    double num = 0.0;
    double Wtot = 0.0;

    for(int i = 0; i < nrow; i++){
        for(int j = 0; j < ncol; j++){
            num += W(i, j) * pow(X[i] - X[j], 2);
            Wtot += W(i, j);
        }
    }

    // Now compute normalization factor 2Wtot*sum((x-mu_x)**2)/(N-1)
    
    // Compute sum((x-mu_x)**2) = xvar

    unsigned int xsize = X.size();
    double xmu = 0.0;
    for(int i = 0; i < xsize; i++){
        xmu += X[i];
    }

    xmu /= xsize;


    double xvar = 0.0;
    for(int i = 0; i < xsize; i++){
        xvar += pow(X[i] - xmu, 2);
    }


    // Full norm factor
    double norm = 2*Wtot*xvar / (xsize-1);

    return num / norm;

}

// [[Rcpp::export]]
double geary_sparse(NumericVector X, NumericMatrix ind, NumericMatrix W) {
    unsigned int nrow = ind.nrow();
    unsigned int ncol = ind.ncol();

    double num = 0.0;
    double Wtot = 0.0;

    double xi, xj;

    for(int i = 0; i < nrow; i++){
        xi = X[i];

        for(int j = 0; j < ncol; j++){
            xj = X[ind(i, j) - 1];
            num += W(i, j) * pow(xi - xj, 2);
            Wtot += W(i, j);
        }

    }

    // Now compute normalization factor 2Wtot*sum((x-mu_x)**2)/(N-1)
    
    // Compute sum((x-mu_x)**2) = xvar

    unsigned int xsize = X.size();
    double xmu = 0.0;
    for(int i = 0; i < xsize; i++){
        xmu += X[i];
    }

    xmu /= xsize;


    double xvar = 0.0;
    for(int i = 0; i < xsize; i++){
        xvar += pow(X[i] - xmu, 2);
    }


    // Full norm factor
    double norm = 2*Wtot*xvar / (xsize-1);

    return num / norm;

}

// [[Rcpp::export]]
NumericVector geary_sparse_local(NumericVector X, NumericMatrix ind, NumericMatrix W) {
    unsigned int nrow = ind.nrow();
    unsigned int ncol = ind.ncol();

    NumericVector local_vals(nrow);

    double num = 0.0;
    double Wtot = 0.0;

    double xi, xj;

    unsigned int xsize = X.size();
    double xmu = 0.0;
    for(int i = 0; i < xsize; i++){
        xmu += X[i];
    }

    xmu /= xsize;

    double xvar = 0.0;
    for(int i = 0; i < xsize; i++){
        xvar += pow(X[i] - xmu, 2);
    }

    for(int i = 0; i < nrow; i++){
        xi = X[i];
        num = 0.0;

        for(int j = 0; j < ncol; j++){
            xj = X[ind(i, j) - 1];
            num += W(i, j) * pow(xi - xj, 2);
        }

        local_vals[i] = num / xvar;
    }

    return local_vals;

}

// [[Rcpp::export]]
NumericVector geary_all(NumericMatrix X, NumericMatrix W) {
    unsigned int nrow = X.nrow();

    NumericVector geary_vals(nrow);

    for(int i = 0; i < nrow; i++){
        MatrixRow<REALSXP> currentRow = X(i, _);
        geary_vals[i] = geary(currentRow, W);
    }

    return geary_vals;

}

// [[Rcpp::export]]
NumericVector geary_sparse_all(NumericMatrix X, NumericMatrix ind, NumericMatrix W) {
    unsigned int nrow = X.nrow();

    NumericVector geary_vals(nrow);

    for(int i = 0; i < nrow; i++){
        MatrixRow<REALSXP> currentRow = X(i, _);
        geary_vals[i] = geary_sparse(currentRow, ind, W);
    }

    return geary_vals;

}
