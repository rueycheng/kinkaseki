#include <iostream>
#include <fstream>
#include <sstream>

#include "magicbox.h"
#include "common.h"
#include "compat.h"

#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/format.hpp>
#include <boost/random.hpp>
#include <ext/algorithm>

//--------------------------------------------------
// instance
//-------------------------------------------------- 
struct feature_vector {
    typedef std::pair<unsigned int, float> component_type;
    typedef std::vector<component_type> vector_type;

    vector_type data;
};

//--------------------------------------------------
// Distance functions
//-------------------------------------------------- 
float cosine_distance_prenormalized(const feature_vector& lhs, const feature_vector& rhs) {
    const feature_vector::vector_type& d1 = lhs.data, d2 = rhs.data;
    feature_vector::vector_type::const_iterator first1 = d1.begin(), end1 = d1.end();
    feature_vector::vector_type::const_iterator first2 = d2.begin(), end2 = d2.end();
    float sum = 0.0;

    while (first1 != end1 && first2 != end2) {
	if (first1->first < first2->first) ++first1;
	else if (first2->first < first1->first) ++first2;
	else sum += (first1++)->second * (first2++)->second;
    }

    return sum;
}

float l1_distance(const feature_vector& lhs, const feature_vector& rhs) {
    const feature_vector::vector_type& d1 = lhs.data, d2 = rhs.data;
    feature_vector::vector_type::const_iterator first1 = d1.begin(), end1 = d1.end();
    feature_vector::vector_type::const_iterator first2 = d2.begin(), end2 = d2.end();
    float sum = 0.0;
    float diff = 0.0;

    while (first1 != end1 && first2 != end2) {
	if (first1->first < first2->first) sum += (first1++)->second;
	else if (first2->first < first1->first) sum += (first2++)->second;
	else {
	    diff = (first1++)->second - (first2++)->second;
	    sum += (diff > 0? diff: -diff);
	}
    }

    while (first1 != end1) sum += (first1++)->second;
    while (first2 != end2) sum += (first2++)->second;

    return sum;
}

float l2_distance(const feature_vector& lhs, const feature_vector& rhs) {
    const feature_vector::vector_type& d1 = lhs.data, d2 = rhs.data;
    feature_vector::vector_type::const_iterator first1 = d1.begin(), end1 = d1.end();
    feature_vector::vector_type::const_iterator first2 = d2.begin(), end2 = d2.end();
    float sum_of_square = 0.0;
    float diff = 0.0;

    while (first1 != end1 && first2 != end2) {
	if (first1->first < first2->first) diff = (first1++)->second;
	else if (first2->first < first1->first) diff = (first2++)->second;
	else diff = (first1++)->second - (first2++)->second;

	sum_of_square += diff * diff;
    }

    while (first1 != end1) {
	diff = (first1++)->second;
	sum_of_square += diff * diff;
    }

    while (first2 != end2) {
	diff = (first2++)->second;
	sum_of_square += diff * diff;
    }

    return std::sqrt(sum_of_square);
}

using namespace std;
using namespace magicbox;

//--------------------------------------------------
// Main program
//-------------------------------------------------- 
int main(int argc, char** argv) {
    // Getopt
    unsigned int num_cluster = 0;
    string distance = "cosine";

    Getopt g(argc, argv);
    g   << $(&num_cluster, "num-cluster,k", "The number of clusters (K)")
	<< $(&distance, "distance", 
		"The distance function: cosine , l1, or l2.  Defaults 'cosine'.")
	<< $$$("[options..]");

    // Go!
    vector<feature_vector> instances;
    vector<unsigned int> labels;
    unsigned int num_feature = 0;

    //--------------------------------------------------
    // Step 1: Parse the input
    //-------------------------------------------------- 
    string line;
    istringstream line_in;

    while (getline(cin, line)) {
	line_in.str(line);
	line_in.clear();

	instances.push_back(feature_vector());
	labels.push_back(0);

	feature_vector& fv = instances.back();

	unsigned int key;
	float value;
	char sep;

	while (line_in >> key >> sep >> value) {
	    if (key > num_feature) num_feature = key;
	    fv.data.push_back(make_pair(key, value));
	}
    }

    //--------------------------------------------------
    // Step 2: Initialize the centroids
    //-------------------------------------------------- 
    vector<feature_vector> centroids;

    random_sample_n(instances.begin(), instances.end(),
	    back_inserter(centroids), num_cluster);

    float (*distance_function)(const feature_vector&, const feature_vector&);
    if (distance == "l1") distance_function = l1_distance;
    else if (distance == "l2") distance_function = l2_distance;
    else distance_function = cosine_distance_prenormalized; 
    // NOTE: I assume the features weights were prenormalized when using cosine
    
    unsigned int num_instance = instances.size();

    for (unsigned int i = 1; ; ++i) {
	float sum_of_errors = 0.0;

	for (unsigned int j = 0; j < num_instance; ++j) {
	    const feature_vector& inst = instances[j];

	    float min_dist = 1.0; // The maximum in the first octant
	    unsigned int min_label = 0;

	    for (unsigned int kk = 0; kk < num_cluster; ++kk) {
		const feature_vector& cent = centroids[kk];
		float current_dist = distance_function(inst, cent);

		if (current_dist < min_dist) {
		    min_dist = current_dist;
		    min_label = kk;
		}
	    }

	    sum_of_errors += min_dist;
	    labels[j] = min_label;
	}

	cerr << "Iteration #" << i << " err=" << sum_of_errors << "\n";
    }

    return 0;
}
