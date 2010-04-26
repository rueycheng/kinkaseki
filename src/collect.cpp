#include <iostream>
#include <sstream>

#include "magicbox.h"
#include "common.h"
#include "compat.h"

#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/format.hpp>
#include <unistd.h>

using namespace std;
using namespace magicbox;

int main(int argc, char** argv) {
    // Getopt
    bool drop_first_token = false;
    unsigned int max_entry = 1024 * 1024;

    Getopt g(argc, argv);
    g   << $(&max_entry, "max-entry,M", "Maximum number of vocabulary stored in memory")
	<< $(&drop_first_token, "drop-first-token", 
		"Drop the first token in each line.  Useful when the first token is used as docno.")
	<< $$$("[options..]");

    //--------------------------------------------------
    // Main program
    //-------------------------------------------------- 
    typedef unordered_map<string, unsigned int> freq_map;
    typedef unordered_set<string> term_set;

    namespace fs = boost::filesystem;

    freq_map df;
    term_set terms;

    //--------------------------------------------------
    // Step 1: Produce multiple disk-runs
    //-------------------------------------------------- 
    unsigned int runno = 0;
    unsigned int pid = getpid();
    string line, word;
    istringstream line_in;

    while (getline(cin, line)) {
	line_in.str(line);
	line_in.clear();

	// Get rid of the first token as necessary
	if (drop_first_token) line_in >> word; 

	// Now is time to collect unique terms in a text unit
	terms.clear();
	while (line_in >> word) terms.insert(word);

	// Put those counts back in
	foreach (const string& term, terms) { ++df[term]; }

	// Flush as necesary
	if (df.size() > max_entry) {
	    string runfile = boost::str(boost::format("collect-%d.run.%05d") % pid % ++runno);
	    fs::ofstream runout(runfile);

	    vector<string> terms;
	    foreach (const freq_map::value_type &p, df) { terms.push_back(p.first); }
	    stable_sort(terms.begin(), terms.end());
	    foreach (const string &term, terms) { runout << term << ' ' << df[term] << '\n'; }
	    df.clear(); // Done
	}
    }

    // Flush as necesary (shameless copy)
    if (!df.empty()) {
	string runfile = boost::str(boost::format("collect-%d.run.%05d") % pid % ++runno);
	fs::ofstream runout(runfile);

	vector<string> terms;
	foreach (const freq_map::value_type &p, df) { terms.push_back(p.first); }
	stable_sort(terms.begin(), terms.end());
	foreach (const string &term, terms) { runout << term << ' ' << df[term] << '\n'; }
	df.clear(); // Done
    }

    //--------------------------------------------------
    // Step 2: Merge disk-runs
    //-------------------------------------------------- 

    return 0;
}
