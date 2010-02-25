#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <cmath>
#include <numeric>

#include "magicbox.h"
#include "common.h"
#include "compat.h"

using namespace std;
using namespace magicbox;

typedef pair<unsigned int,unsigned short> posting;
typedef vector<posting> posting_list;
typedef unordered_map<string, posting_list> inverted_table;
typedef unordered_map<string, unsigned short> freq_map;

namespace std {
    ostream& operator<<(ostream& out, const posting_list& pl) {
	foreach (const posting& p, pl) { out << p.first << ' ' << p.second << ' '; }
	return out;
    }
}

//--------------------------------------------------
// Main program
//-------------------------------------------------- 
int main(int argc, char** argv) {

    // Getopt
    vector<string> input;
    string model = "model.unnamed";

    Getopt g(argc, argv);
    g   << $(&model, "model,m", "Specify the model directory")
	<< $(&input, "input", "", -1)
	<< $$$("files..");

    if (input.empty()) input.push_back("-");

    // Create model directory
    namespace fs = boost::filesystem;

    fs::create_directory(model);
    if (!fs::exists(model)) die("Cannot create directory " + model);
    fs::path basedir(model);

    // Bookkeeping..
    vector<string> docnos;
    inverted_table inv;

    foreach (const string& filename, input) {
	AutoIn in(filename);
	string line;

	while (getline(in(), line)) {
	    istringstream iss(line);
	    string docno;
	    iss >> docno;

	    if (!boost::ends_with(docno, ":facet")) {
		// Case 1: Regular texts

		docnos.push_back(docno);
		unsigned int docid = docnos.size();

		{
		    // Get counts
		    freq_map counts;
		    string word;
		    while (iss >> word) ++counts[word];

		    foreach (const freq_map::value_type& p, counts) {
			inv[p.first].push_back(posting(docid, p.second));
		    }
		}
	    }
	    else {
		// Case 2: Facets
	    }
	}
    }

    // Sort the keys and output
    vector<string> keys;
    foreach (const inverted_table::value_type& entry, inv) { keys.push_back(entry.first); }
    stable_sort(keys.begin(), keys.end());
    foreach (const string& key, keys) { cout << key << ' ' << inv[key] << '\n'; }

    return 0;
}
