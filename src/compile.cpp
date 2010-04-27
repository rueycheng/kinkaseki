#include <iostream>
#include <fstream>
#include <sstream>
#include <numeric>

#include "magicbox.h"
#include "common.h"
#include "compat.h"

#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/format.hpp>

//--------------------------------------------------
// record (client code)
//-------------------------------------------------- 
struct record {
    std::string key;
    unsigned int count;
};

bool operator<(const record& lhs, const record& rhs) {
    return lhs.key < rhs.key;
}

std::istream& operator>>(std::istream& in, record& r) {
    std::string line;
    
    if (std::getline(in, line)) {
	std::istringstream v(line);
	v >> r.key >> r.count;
    }

    return in;
}

std::ostream& operator<<(std::ostream& out, const record& r) {
    return out << r.key << ' ' << r.count;
}

//--------------------------------------------------
// Normalizers
//-------------------------------------------------- 
template<typename InputIterator>
void l1_normalization(InputIterator first, InputIterator last) {
    typename InputIterator::value_type norm = std::accumulate(first, last, 0.0);
    while (first != last) *first++ /= norm;
}

template<typename InputIterator>
void l2_normalization(InputIterator first, InputIterator last) {
    typename InputIterator::value_type norm = std::sqrt(std::inner_product(first, last, first, 0.0));
    while (first != last) *first++ /= norm;
}

template<typename InputIterator>
void inf_normalization(InputIterator first, InputIterator last) {
    InputIterator e = std::max_element(first, last);
    typename InputIterator::value_type norm = (e == last)? 0: *e;
    while (first != last) *first++ /= norm;
}

template<typename InputIterator>
void no_op(InputIterator first, InputIterator last) { }

using namespace std;
using namespace magicbox;

//--------------------------------------------------
// Main program
//-------------------------------------------------- 
int main(int argc, char** argv) {
    // Getopt
    string eol_token = "#eol";
    string df_file;
    string normalization;

    Getopt g(argc, argv);
    g   << $(&eol_token, "eol-token", "The special token indicating the end of line")
	<< $(&df_file, "df-file", "The document-frequency map used to compute idf")
	<< $(&normalization, "normalization", "The normalization method: l1, l2, or inf.  Defaults 'no-op'")
	<< $$$("[options..]");

    namespace fs = boost::filesystem;

    // Go!

    //--------------------------------------------------
    // Step 1: Load the document frequencies back in
    //-------------------------------------------------- 
    typedef unordered_map<string, unsigned int> freq_map;

    freq_map term_no;
    vector<string> vocab_term;
    vector<unsigned int> vocab_df;

    {
	fs::ifstream dfin(df_file);
	record r;

	while (dfin >> r) {
	    vocab_term.push_back(r.key);
	    vocab_df.push_back(r.count);
	    term_no[r.key] = vocab_term.size() - 1;
	}
    }

    unsigned int N = vocab_df.at(term_no[eol_token]);
    float logN = log(N);

    //--------------------------------------------------
    // Step 2: Process input line-by-line
    //-------------------------------------------------- 
    typedef unordered_set<unsigned int> term_set;

    string line, word;
    istringstream line_in;
    term_set terms;
    vector<unsigned int> term_id;
    vector<float> term_idf;

    // Set up normalizer
    typedef void (*normalizer)(vector<float>::iterator, vector<float>::iterator);
    normalizer norm;

    if (normalization == "l1") norm = l1_normalization;
    else if (normalization == "l2") norm = l2_normalization;
    else if (normalization == "inf") norm = inf_normalization;
    else norm = no_op;

    while (getline(cin, line)) {
	line_in.str(line);
	line_in.clear();

	terms.clear();
	while (line_in >> word) 
	    if (term_no.find(word) != term_no.end()) terms.insert(term_no[word]);

	term_id.clear();
	term_idf.clear();
	copy(terms.begin(), terms.end(), back_inserter(term_id));
	sort(term_id.begin(), term_id.end());
	foreach (const unsigned int id, term_id) {
	    term_idf.push_back(logN - log(vocab_df[id]));
	}

	// Apply the normalizer
	norm(term_idf.begin(), term_idf.end());

	unsigned int num_term = term_id.size();
	for (unsigned int i = 0; i < num_term; ++i)
	    cout << term_id[i] << ':' << term_idf[i] << ' ';

	cout << '\n';
    }

    return 0;
}
