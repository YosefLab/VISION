// Utilities for fast matrix operations to be wrapped with R in FastProjectR

#include <Rcpp.h> 
#include "vptree.h"
#include <iostream>
#ifdef _OPENMP
    #include <omp.h>
    //[[Rcpp::plugins(openmp)]]
#endif
#include <string>
#include <sstream>
using namespace Rcpp;


struct PublicDistanceComparator
{
    const DataPoint& item;
    PublicDistanceComparator(const DataPoint& item) : item(item) {}
    bool operator()(const DataPoint& a, const DataPoint& b) {
        return euclidean_distance(item, a) < euclidean_distance(item, b);
    }
};


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
List ball_tree_vector_knn(NumericMatrix X, NumericVector Y, int K, int n_threads) {
	
	int N = X.nrow();
	int D = X.ncol();

	double* data;

	data = (double*) calloc(D*N, sizeof(double));
	if (data == NULL) { Rcpp::stop("memory allocation failed!\n"); }
	for (int i = 0; i < N; i++) {
		for (int j=0; j < D; j++) {
			data[i*D+j] = X(i,j);
		}
	}

	double* data2;
	
	data2 = (double*) calloc(D, sizeof(double));
	for (int i = 0; i < D; i++) {
		data2[i] = Y[i];
	}

	DataPoint obj_Y = DataPoint(D, -1, data2);


	// build ball tree on set
	VpTree<DataPoint, euclidean_distance>* tree = new VpTree<DataPoint, euclidean_distance>();
	std::vector<DataPoint> obj_X(N, DataPoint(D, -1, data));
	for (int n = 0; n < N; n++) {
		obj_X[n] = DataPoint(D, n, data + n * D);
	}
	tree -> create(obj_X, 0);


	// find nearest neighbors
	Rcpp::NumericMatrix Ind(1, K);
	Rcpp::NumericMatrix Dist(1, K);

	std::vector<int> ind_arr(K+1, 0);
	std::vector<double> dist_arr(K+1, 0.0);
			
	tree -> search(obj_Y, K+1, &ind_arr, &dist_arr);

	NumericVector ind(ind_arr.begin(), ind_arr.end());
	NumericVector dist(dist_arr.begin(), dist_arr.end());
	ind.erase(0); dist.erase(0);
	Ind(0,_) = ind;
	Dist(0,_) = dist;


	obj_X.clear();
	delete tree;


	List ret;
	ret["index"] = Ind;
	ret["dist"] = Dist;

	return ret;	
}

// [[Rcpp::export]]
List ball_tree_knn(NumericMatrix X ,int K, int n_threads) {

	int N = X.nrow();
	int D = X.ncol();

	double* data;


	data = (double*) calloc(D*N, sizeof(double));
	if (data == NULL) { Rcpp::stop("memory allocation failed!\n"); }
	for (int i = 0; i < N; i++) {
		for (int j=0; j < D; j++) {
			data[i*D+j] = X(i,j);
		}
	}


	// build ball tree on set
	VpTree<DataPoint, euclidean_distance>* tree = new VpTree<DataPoint, euclidean_distance>();
	std::vector<DataPoint> obj_X(N, DataPoint(D, -1, data));
	for (int n = 0; n < N; n++) {
		obj_X[n] = DataPoint(D, n, data + n * D);
	}
	tree -> create(obj_X, 0);


	// find nearest neighbors
	Rcpp::NumericMatrix Ind(N, K);
	Rcpp::NumericMatrix Dist(N, K);

	std::vector< std::vector<int> > ind_arr(N, std::vector<int>(K+1));
	std::vector< std::vector<double> > dist_arr(N, std::vector<double>(K+1));
	
    #ifdef _OPENMP
	    omp_set_num_threads(n_threads);
	
        #pragma omp parallel
	    {
	    #pragma omp for 
	    for (int n = 0; n < N; n++) {
		    tree -> search(obj_X[n], K+1, &ind_arr[n], &dist_arr[n]);
	    }
	    }
    #else 
	    for (int n = 0; n < N; n++) {
		    tree -> search(obj_X[n], K+1, &ind_arr[n], &dist_arr[n]);
	    }
    #endif


	for (int i = 0; i < N; i++) {
		NumericVector ind(ind_arr[i].begin(), ind_arr[i].end());
		NumericVector dist(dist_arr[i].begin(), dist_arr[i].end());
		ind.erase(0); dist.erase(0);
		Ind(i,_) = ind;
		Dist(i,_) = dist;
	}


	obj_X.clear();
	delete tree;


	List ret;
	ret["index"] = Ind;
	ret["dist"] = Dist;

	return ret;	
}

// [[Rcpp::export]]
NumericMatrix load_in_knn(NumericMatrix nn, NumericMatrix d) {
	/* An algorithm to load in the results from calling KNN into an R matrix.
	 * Takes the distances of each knn stored in D and places them in the indices specified by NN
	 *
	 * Returns: NumericMatrix (N x N) where N is the number of cells/samples of nn and d
	 */

	int K = nn.ncol();
	int N = nn.nrow();

	NumericMatrix result(N, N);

	for (int i = 0; i < N; i++) {
		for (int j = 0; j < K; j++) {
			result(i, nn(i,j)-1) = d(i, j);
		}
	}

	return result;
}

