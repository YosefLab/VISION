/*
 *  vptree.h
 *  Implementation of a vantage-point tree.
 *
 *  Created by Laurens van der Maaten.
 *  Copyright 2012, Delft University of Technology. All rights reserved.
 *
 *  FastProject version by Matthew Jones, 2017. matthew.jones@ucsf.edu
 */


#include <stdlib.h>
#include <algorithm>
#include <math.h>
#include <vector>
#include <stdio.h>
#include <queue>
#include <limits>
#include <cctype>
#include <string>
#include <Rcpp.h>

#ifndef VPTREE_H
#define VPTREE_H

#define PRINT Rprintf

class DataPoint
{
    int _D;
    int _ind;
    double* _x;

public:
    DataPoint() {
        _D = 1;
        _ind = -1;
        _x = NULL;
    }
    DataPoint(int D, int ind, double* x) {
        _D = D;
        _ind = ind;
        _x = (double*) malloc(_D * sizeof(double));
        for (int d = 0; d < _D; d++) _x[d] = x[d];
    }
    DataPoint(const DataPoint& other) {                     // this makes a deep copy -- should not free anything
        if (this != &other) {
            _D = other.dimensionality();
            _ind = other.index();
            _x = (double*) malloc(_D * sizeof(double));
            for (int d = 0; d < _D; d++) _x[d] = other.x(d);
        }
    }
    ~DataPoint() { if (_x != NULL) free(_x); }
    DataPoint& operator= (const DataPoint& other) {         // assignment should free old object
        if (this != &other) {
            if (_x != NULL) free(_x);
            _D = other.dimensionality();
            _ind = other.index();
            _x = (double*) malloc(_D * sizeof(double));
            for (int d = 0; d < _D; d++) _x[d] = other.x(d);
        }
        return *this;
    }
    int index() const { return _ind; }
    int dimensionality() const { return _D; }
    double x(int d) const { return _x[d]; }
};


double euclidean_distance(const DataPoint &t1, const DataPoint &t2) {
    double dd = .0;
    for (int d = 0; d < t1.dimensionality(); d++) {
        dd += (t1.x(d) - t2.x(d)) * (t1.x(d) - t2.x(d));
    }
    return sqrt(dd);
}

double single_euclidean_distance(int d, const DataPoint &t1, const DataPoint &t2) {
	double dd = (t1.x(d) - t2.x(d)) * (t1.x(d) - t2.x(d));
	return sqrt(dd);
}


template<typename T, double (*distance)( const T&, const T& )>
class VpTree
{
public:

    // Default constructor
    VpTree() : _root(0) {}

    // Destructor
    ~VpTree() {
        delete _root;
    }

    // Function to create a new VpTree from data
    void create(const std::vector<T>& items, int cluster) {
        delete _root;
        _items = items;
        if (cluster == 1) {
			_root = buildFromPoints_Cluster(0, items.size());
		} else {
			_root = buildFromPoints_KNN(0, items.size());
		}
    }

    // Function that uses the tree to find the k nearest neighbors of target
    void search(const T& target, unsigned int k, std::vector<int>* results, std::vector<double>* distances)
    {

        // Use a priority queue to store intermediate results on
        std::priority_queue<HeapItem> heap;

        // Variable that tracks the distance to the farthest point in our results
        double tau = DBL_MAX;

        // Perform the search
        search(_root, target, k, heap, tau);

        // Gather final results
        results->clear(); distances->clear();
        while (!heap.empty()) {
            results->push_back(_items[heap.top().index].index()+1);
            distances->push_back(heap.top().dist);
            heap.pop();
        }

        // Results are in reverse order
        std::reverse(results->begin(), results->end());
        std::reverse(distances->begin(), distances->end());
    }

    void find_partitions(std::vector< std::vector<int> >* clusters, std::vector<double>* cluster_radii, int L) {

		std::vector< std::vector<int> > results;
		find_partitions(_root, &results, cluster_radii, L);	
		
		// Gather results
		for (unsigned int i = 0; i < results.size(); i++) {
			std::vector<int> clust;
			for (unsigned int j = 0; j < results.at(i).size(); j++) {
				clust.push_back(_items[results.at(i).at(j)].index() + 1);
			}
			clusters->push_back(clust);
		}


		PRINT("%d\n", clusters->size());

	}

private:
    std::vector<T> _items;

    // Single node of a VP tree (has a point and radius; left children are closer to point than the radius)
    struct Node
    {
        int index;              // index of point in node
        double threshold;       // radius(?)
        int num_points;		// Number of points contained within this node
        Node* left;             // points closer by than threshold
        Node* right;            // points farther away than threshold

        Node() :
            index(0), threshold(0.), num_points(0), left(0), right(0) {}

        ~Node() {               // destructor
            delete left;
            delete right;
        }
    }* _root;



    // An item on the intermediate result queue
    struct HeapItem {
        HeapItem( int index, int num_points, double dist) :
            index(index), num_points(num_points), dist(dist) {}
        int index;
        double num_points;
        double dist;
        bool operator<(const HeapItem& o) const {
            return dist < o.dist;
        }
    };

    // Distance comparator for use in std::nth_element
    struct DistanceComparator
    {
        const T& item;
        DistanceComparator(const T& item) : item(item) {}
        bool operator()(const T& a, const T& b) {
            return distance(item, a) < distance(item, b);
        }
    };

    struct DistanceComparator_Single {
		const T& item;
		int d;
		DistanceComparator_Single(const T& item, int d) : item(item), d(d) {}
		bool operator()(const T& a, const T&b) {
			return single_euclidean_distance(d, item, a) < single_euclidean_distance(d, item, b);
		}
	};

    // Function that (recursively) fills the tree
    Node* buildFromPoints_KNN( int lower, int upper )
    {
        if (upper == lower) {     // indicates that we're done here!
            return NULL;
        }

        // Lower index is center of current node
        Node* node = new Node();
        node->index = lower;
        node->num_points = 1; 

        if (upper - lower > 1) {      // if we did not arrive at leaf yet

            // Choose an arbitrary point and move it to the start
            int i = (int) ((double)R::runif(0,1) * (upper - lower - 1)) + lower;
            std::swap(_items[lower], _items[i]);

            // Partition around the median distance
            int median = (upper + lower) / 2;
            std::nth_element(_items.begin() + lower + 1,
                             _items.begin() + median,
                             _items.begin() + upper,
                             DistanceComparator(_items[lower]));

            // Threshold of the new node will be the distance to the median
            node->threshold = distance(_items[lower], _items[median]);

            // Recursively build tree
            node->num_points = upper-lower;
            node->index = lower;
            node->left = buildFromPoints_KNN(lower+1, median);
            node->right = buildFromPoints_KNN(median, upper);
        }

        // Return result
        return node;
    }

    // Function that (recursively) fills the tree
    Node* buildFromPoints_Cluster( int lower, int upper )
    {

		Rcpp::Environment base("package:stats");
		Rcpp::Function kmeans_r = base["kmeans"];

        if (upper == lower) {     // indicates that we're done here!
            return NULL;
        }

        // Lower index is center of current node
        Node* node = new Node();
        node->index = lower;
        node->num_points = 1; 

        if (upper - lower > 1) {      // if we did not arrive at leaf yet

            // Choose an arbitrary point and move it to the start
            int i = (int) ((double)R::runif(0,1) * (upper - lower - 1)) + lower;
            std::swap(_items[lower], _items[i]);

            //DataPoint elem = _items[lower];
			//int max_dim = 0;
			//double max_dist = 0;
			int median = (upper + lower) / 2;


			// Find dimensionality with greatest extent to partition around	
			/*for (int d = 0; d < elem.dimensionality(); d++) {
				std::vector<DataPoint>::iterator max_elem = std::max_element(_items.begin() + lower,
																   _items.begin() + upper,
																   DistanceComparator_Single(elem, d));
				double dist = single_euclidean_distance(d, elem, *max_elem);
				if (dist > max_dist) {
					max_dist = dist;
					max_dim = d;
				}
			}*/
			
            // Partition around the median distance
            std::nth_element(_items.begin() + lower,
                             _items.begin() + median,
                             _items.begin() + upper,
                            DistanceComparator(_items[lower]));
			
			// Partition around median distance in dimension of greatest extent
			/*std::nth_element(_items.begin() + lower,
							 _items.begin() + median, 
							 _items.begin() + upper,
							 DistanceComparator_Single(elem, max_dim));
			*/
			

            // Threshold of the new node will be the distance to the median
            node->threshold = distance(_items[lower], _items[median]);

            // Recursively build tree
            node->num_points = median-lower;
            node->index = lower;
            node->left = buildFromPoints_Cluster(lower, median);
            node->right = buildFromPoints_Cluster(median, upper);
        }

        // Return result
        return node;
    }

	// Helper function that partitions the tree according to the threshold L, starting from node NODE 
	void find_partitions(Node* node, std::vector< std::vector<int> >* clusters, std::vector<double>* cluster_radii, int L) {

		if (node == NULL) return;	// indicates that we're done here
		if (node->right == NULL && node->left==NULL) {
			std::vector<int> partition;
			partition.push_back(node->index);
			clusters->push_back(partition);
		}	

		// If node meets threshold criteria, fill partition with all of points in node
		if (node->num_points <= L && node->num_points > 1) {
			cluster_radii->push_back(node->threshold);
			std::vector<int> partition;
			fill_partition(node->left, &partition);
			clusters->push_back(partition);
			find_partitions(node->right, clusters, cluster_radii, L);
		} else {	
			// If node didn't meet criteria, recursively partition tree to get rest of nodes
			find_partitions(node->left, clusters, cluster_radii, L);
			find_partitions(node->right, clusters, cluster_radii, L);
		}
	
	}

	// Helper function to fill a partition once a node passes the threshold criteria
	void fill_partition(Node* node, std::vector<int>* partition) {

		if (node == NULL) { return; }  
		if (node->right == NULL && node->left == NULL) {
			partition->push_back(node->index);
		}

		// Now recursively call function for all of the node's children
		fill_partition(node->left, partition);
		fill_partition(node->right, partition);

	}

    // Helper function that searches the tree
    void search(Node* node, const T& target, int k, std::priority_queue<HeapItem>& heap, double& tau)
    {
        if (node == NULL) return;    // indicates that we're done here

        // Compute distance between target and current node
        double dist = distance(_items[node->index], target);

        // If current node within radius tau
        if (dist < tau) {
            if (heap.size() == k) heap.pop();                // remove furthest node from result list (if we already have k results)
            heap.push(HeapItem(node->index, node->num_points, dist));           // add current node to result list
            if (heap.size() == k) tau = heap.top().dist;    // update value of tau (farthest point in result list)
        }

        // Return if we arrived at a leaf
        if (node->left == NULL && node->right == NULL) {
            return;
        }

        // If the target lies within the radius of ball
        if (dist < node->threshold) {
            if (dist - tau <= node->threshold) {        // if there can still be neighbors inside the ball, recursively search left child first
                search(node->left, target, k, heap, tau);
            }

            if (dist + tau >= node->threshold) {        // if there can still be neighbors outside the ball, recursively search right child
                search(node->right, target, k, heap, tau);
            }

            // If the target lies outsize the radius of the ball
        } else {
            if (dist + tau >= node->threshold) {        // if there can still be neighbors outside the ball, recursively search right child first
                search(node->right, target, k, heap, tau);
            }

            if (dist - tau <= node->threshold) {         // if there can still be neighbors inside the ball, recursively search left child
                search(node->left, target, k, heap, tau);
            }
        }
    }
};

#endif
