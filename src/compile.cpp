#include <iostream>
#include <fstream>
#include <sstream>

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

using namespace std;
using namespace magicbox;

//--------------------------------------------------
// Main program
//-------------------------------------------------- 
int main(int argc, char** argv) {
    // Getopt
    string eol_token = "#eol";
    string df_file = "";

    Getopt g(argc, argv);
    g   << $(&eol_token, "eol-token", "The special token indicating the end of line")
	<< $(&df_file, "df-file", "The document-frequency map used to compute idf")
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
    double logN = log(N);

    //--------------------------------------------------
    // Step 2: Process input line-by-line
    //-------------------------------------------------- 
    typedef unordered_set<unsigned int> term_set;

    string line, word;
    istringstream line_in;
    term_set terms;

    while (getline(cin, line)) {
	line_in.str(line);
	line_in.clear();

	terms.clear();
	while (line_in >> word) 
	    if (term_no.find(word) != term_no.end()) terms.insert(term_no[word]);

	vector<unsigned int> terms_sorted;
	copy(terms.begin(), terms.end(), back_inserter(terms_sorted));
	sort(terms_sorted.begin(), terms_sorted.end());

	double idf = 0;
	foreach (const unsigned int id, terms_sorted) {
	    idf = logN - log(vocab_df[id]);
	    cout << id << ':' << idf << ' ';
	}
	cout << '\n';
    }

    return 0;
}
