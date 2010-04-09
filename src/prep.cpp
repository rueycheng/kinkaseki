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
    typedef unordered_map<string, unsigned int> freq_map;
    typedef unordered_set<string> term_set;

    freq_map df;
    term_set terms;

    string line, word;
    istringstream line_in;

    while (getline(cin, line)) {
	// Now is time to collect unique terms in a text unit
	terms.clear();

	line_in.str(line);
	line_in.clear();
	line_in >> word; // Get rid of the first token, i.e., ID
	while (line_in >> word) terms.insert(word);

	// Put those counts back in
	foreach (const string& term, terms) { ++df[term]; }
    }

    return 0;
}
