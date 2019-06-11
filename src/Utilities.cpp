// Utilities for fast matrix operations to be wrapped with R in VISION

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

// [[Rcpp::export]]
NumericVector sigGeneInner(
        const NumericVector& sigScores,
        const NumericMatrix& exp,
        const NumericVector& geneIndices) {

    int n = geneIndices.size();
    NumericVector corr(n);

    int m = sigScores.size();

    for (int i = 0; i < n; i++) {

        double total = 0;
        int k = geneIndices[i]-1;

        for (int j = 0; j < m; j++) {
            total += exp(k, j) * sigScores[j];
        }

        corr[i] = total / (m-1);
    }

    return corr;
}

// Runs the ranksums test on two subsets of the values in 'vals'
// Ties are handled correctly, but NAs are not checked
// [[Rcpp::export]]
List wilcox_subset(
        NumericVector vals,
        NumericVector indicesA,
        NumericVector indicesB) {

    int N1 = indicesA.size();
    int N2 = indicesB.size();
    unsigned int N = N1 + N2;

    double *sortedA = new double[N1];
    for (int i = 0; i < N1; i++){
        sortedA[i] = vals[indicesA[i]-1];
    }

    double *sortedB = new double[N2];
    for (int i = 0; i < N2; i++){
        sortedB[i] = vals[indicesB[i]-1];
    }

    std::sort(sortedA, sortedA+N1);
    std::sort(sortedB, sortedB+N2);

    double *ranks = new double[N];
    double *fromA = new double[N];

    unsigned int pointerA = 0;
    unsigned int pointerB = 0;

    double last_val = 0.0;
    double current_val = 0.0;

    unsigned int tie_start = 0;
    bool in_tie = false;

    double tie_groups = 0; // sum of ti^3 - ti

    double a, b;
    double start_i, end_i, mean_i;
    bool val_is_A;
    unsigned int tie_size;

    for(int i = 0; i < N; i++){

        if (pointerA >= N1){
            b = sortedB[pointerB];
            current_val = b;
            val_is_A = false;
        } else if (pointerB >= N2){
            a = sortedA[pointerA];
            current_val = a;
            val_is_A = true;
        } else {
            a = sortedA[pointerA];
            b = sortedB[pointerB];

            if (a < b) {
                current_val = a;
                val_is_A = true;
            } else {
                current_val = b;
                val_is_A = false;
            }

        }

        // If a tie has started
        if (i > 0 && current_val == last_val && !in_tie){

            tie_start = i - 1;
            in_tie = true;

        }

        if (val_is_A) {

            fromA[i] = 1;
            pointerA = pointerA + 1;

        } else {

            fromA[i] = 0;
            pointerB = pointerB + 1;
        }

        ranks[i] = i+1;

        // Tie has ended
        if (in_tie && current_val != last_val) {

            start_i = ranks[tie_start];
            end_i = ranks[i-1];

            mean_i = (start_i + end_i) / 2;

            for (int j = tie_start; j <= i-1; j++){
                ranks[j] = mean_i;
            }

            in_tie = false;

            tie_size = i - tie_start;
            tie_groups = tie_groups + pow(tie_size, 3) - tie_size;
        }

        last_val = current_val;

    }

    // If we ended in a tie
    if (in_tie) {
        start_i = ranks[tie_start];
        end_i = ranks[N-1];

        mean_i = (start_i + end_i) / 2;

        for (int j = tie_start; j <= N-1; j++){
            ranks[j] = mean_i;
        }

        tie_size = N - tie_start;
        tie_groups = tie_groups + pow(tie_size, 3) - tie_size;

        in_tie = FALSE;
    }

    // Compute RS_A
    double RS_A = 0;
    for (int i = 0; i < N; i++){
        if (fromA[i] == 1){
            RS_A = RS_A + ranks[i];
        }
    }

    double UA = RS_A - N1 * (N1 + 1) / 2;
    double UB = N1 * N2 - UA;

    double U = std::min(UA, UB);

    double mU = N1 * N2 / 2;

    double sdU = (double)(N1 * N2 / 12) * (N + 1 - (double)tie_groups / (N * (N - 1)));
    sdU = sqrt(sdU);

    double z = (U - mU) / sdU;

    z = z + .5 / sdU; // continuity correction

    // NumericVector zV(1);
    // zV[0] = z;

    // NumericVector p = pnorm(zV, 0.0, 1.0);


    List out(2);
    out.names() = CharacterVector::create("U", "Z");
    out[0] = U;
    out[1] = z;


    delete[] sortedA;
    delete[] sortedB;
    delete[] ranks;
    delete[] fromA;

	return out;
}
