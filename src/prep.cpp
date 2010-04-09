#include <iostream>
#include <sstream>

#include "magicbox.h"
#include "common.h"
#include "compat.h"

using namespace std;
using namespace magicbox;

int main(int argc, char** argv) {

    // Getopt
    unsigned int max_entry = 1024 * 1024;
    vector<string> input;

    Getopt g(argc, argv);
    g   << $(&max_entry, "max-entry,M", "Maximum number of vocabulary stored in memory")
	<< $(&input, "input", "", -1)
	<< $$$("[options..]");

    // Main program
    typedef unordered_map<string,unsigned int> freq_map;
    freq_map df;

    string line, word;
    while (getline(cin, line)) {
	istringstream in(line);

	in >> word; // The first token is never used here

	if (df.size() > max_entry) {
	//--------------------------------------------------
	//     vector<string> keys;
	//     foreach (freq_map::value_type& p, df) { keys.push_back(p.first); }
	//     sort(keys.begin(), keys.end());
	//     foreach (string key, keys) { cout << key << ' ' << df[key] << '\n'; }
	//-------------------------------------------------- 
	    cout << "Flush" << '\n';

	    df.clear();
	}
	    
	while (in >> word) ++df[word];
    }

    if (df.size()) {
	//--------------------------------------------------
	// vector<string> keys;
	// foreach (freq_map::value_type& p, df) { keys.push_back(p.first); }
	// sort(keys.begin(), keys.end());
	// foreach (string key, keys) { cout << key << ' ' << df[key] << '\n'; }
	//-------------------------------------------------- 
	cout << "Flush" << '\n';

	df.clear();
    }

    return 0;
}
