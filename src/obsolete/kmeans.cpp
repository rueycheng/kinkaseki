#include <iostream>
#include <fstream>
#include <sstream>
#include <numeric>
#include <limits>

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

    return 1.0 - sum;
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
    float epsilon = 0.0001;

    Getopt g(argc, argv);
    g   << $(&num_cluster, "num-cluster,k", "The number of clusters (K).")
	<< $(&distance, "distance", 
		"The distance function: cosine , l1, or l2.  Defaults 'cosine'.")
	<< $(&epsilon, "epsilon,e", "The lowest tolerable error.  Defaults '0.0001'.")
	<< $("use-kmedoids", 
		"Employ k-medoids algorithm.\n" 
		"In other words, new centroids are chosen as the instance closest to the cluster center."
		)
	<< $$$("[options..]");

    // Go!
    if (num_cluster <= 0) die("The number of cluster (--num-cluster,-k) should be greater than 0");

    // Set up the normalizer and distance function
    float (*distance_function)(const feature_vector&, const feature_vector&);
    bool normalized_by_l2 = false;

    // NOTE: I assume the features weights were prenormalized when using cosine
    if (distance == "l1") distance_function = l1_distance;
    else if (distance == "l2") distance_function = l2_distance;
    else {
	distance_function = cosine_distance_prenormalized; 
	normalized_by_l2 = true;
    }
    
    vector<feature_vector> instances;
    vector<unsigned int> labels;
    unsigned int num_feature = 0;

    //--------------------------------------------------
    // Step 1: Parse the input
    //-------------------------------------------------- 
    string line;
    istringstream line_in;

    cerr << "kmeans: Process the input" << '\n';

    unsigned int key;
    float value;
    char sep;
    feature_vector fv;

    while (getline(cin, line)) {
	line_in.str(line);
	line_in.clear();

	fv.data.clear();
	while (line_in >> key >> sep >> value) {
	    if (key > num_feature) num_feature = key;
	    fv.data.push_back(make_pair(key, value));
	}

	if (!fv.data.empty()) {
	    instances.push_back(fv);
	    labels.push_back(0);
	}
    }

    //--------------------------------------------------
    // Step 2: Initialize the centroids
    //-------------------------------------------------- 
    vector<feature_vector> centroids;

    random_sample_n(instances.begin(), instances.end(),
	    back_inserter(centroids), num_cluster);

    cerr << "kmeans: Select arbitrary centroids" << '\n';

    unsigned int num_instance = instances.size();

    // Start with maximized error
    float sum_of_errors_prev = std::numeric_limits<float>::max();

    //--------------------------------------------------
    // Step 3: Enter the main loop
    //-------------------------------------------------- 
    cerr << "kmeans: num_cluster=" << num_cluster << '\n';
    cerr << "kmeans: num_instance=" << num_instance << '\n';
    cerr << "kmeans: num_feature=" << num_feature << '\n';

    for (unsigned int i = 1; ; ++i) {
	float sum_of_errors = 0.0;

	// Start with empty bins
	vector<vector<unsigned int> > assignments(num_cluster, vector<unsigned int>());

	for (unsigned int j = 0; j < num_instance; ++j) {
	    const feature_vector& inst = instances[j];

	    float min_dist = std::numeric_limits<float>::max(); // The maximum in the first octant
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

	    // Push the instance ID into the right bin (important!)
	    assignments[min_label].push_back(j);
	}

	cerr << "kmeans: Estimation #" << i << " err=" << sum_of_errors << "\n";

	//--------------------------------------------------
	// for (unsigned int k = 0; k < num_cluster; ++k) 
	//     cerr << "kmeans: Cluster #" << i << " size=" << assignments[k].size() << '\n';
	//-------------------------------------------------- 

	// Recalculate the centroids based on the content of bins
	for (unsigned int k = 0; k < num_cluster; ++k) {
	    const vector<unsigned int>& bin = assignments[k];

	    // HACK: Employ dense representation for fast computation
	    vector<float> rep(num_feature + 1, 0.0); // NOTE: Feature ID begins at 1
	    feature_vector new_centroid;

	    foreach (const unsigned int j, bin) {
		const feature_vector& inst = instances[j];
		foreach (const feature_vector::component_type& comp, inst.data) {
		    rep[comp.first] += comp.second;
		}
	    }

	    // Produce the sparse representation
	    unsigned int bin_size = bin.size();
	    for (unsigned int t = 0; t < num_feature; ++t) {
		if (rep[t] == 0.0) continue;
		new_centroid.data.push_back(make_pair(t, rep[t] / bin_size));
	    }

	    if (normalized_by_l2) {
		float sum_of_square = 0.0;
		foreach (const feature_vector::component_type& comp, new_centroid.data) {
		    sum_of_square += comp.second * comp.second;
		}

		sum_of_square = std::sqrt(sum_of_square);
		foreach (feature_vector::component_type& comp, new_centroid.data) {
		    comp.second /= sum_of_square;
		}
	    }

	    if (g["use-kmedoids"]) {
		unsigned int nearest_dist = std::numeric_limits<float>::max();
		unsigned int nearest_instance = bin[0];  // Dirty!

		foreach (const unsigned int j, bin) {
		    float dist = distance_function(instances[j], new_centroid);
		    if (dist < nearest_dist) {
			nearest_dist = dist;
			nearest_instance = j;
		    }
		}

		// Now use the nearest instance instead
		new_centroid = instances[nearest_instance];
	    }

	    // Write back
	    centroids[k] = new_centroid;
	}

	cerr << "kmeans: Maximization #" << i << "\n";

	if ((sum_of_errors_prev - sum_of_errors >= 0)
		&& (sum_of_errors_prev - sum_of_errors < epsilon)) break; // Congrats!
	sum_of_errors_prev = sum_of_errors;
    }

    foreach (const unsigned int l, labels) {
	cout << l << '\n';
    }

    return 0;
}
