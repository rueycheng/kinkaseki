#include <iostream>
#include <fstream>
#include <sstream>

#include "magicbox.h"
#include "common.h"
#include "compat.h"

#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/format.hpp>
#include <boost/ptr_container/ptr_vector.hpp>
#include <unistd.h>

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
// ifstream_run
//-------------------------------------------------- 
typedef std::istream_iterator<record> ifstream_run;
bool operator<(const ifstream_run& lhs, const ifstream_run& rhs) { 
    return (*lhs) < (*rhs); 
}

bool operator>(const ifstream_run& lhs, const ifstream_run& rhs) { 
    return (*rhs) < (*lhs); 
}

//--------------------------------------------------
// ifstream_pool
//-------------------------------------------------- 
typedef boost::ptr_vector<std::ifstream> ifstream_pool;

//--------------------------------------------------
// multiway_merge (algorithm)
//-------------------------------------------------- 
template<typename RandomAccessIterator, typename OutputIterator, typename StrictWeakOrdering>
void multiway_merge(RandomAccessIterator _first, RandomAccessIterator _last, 
	OutputIterator _result, StrictWeakOrdering _comp) {
    typedef typename RandomAccessIterator::value_type Run;

    Run _eor;
    std::make_heap(_first, _last, _comp);
    while (_first != _last) {
	std::pop_heap(_first, _last, _comp);
	Run& _curr = *(_last - 1);
	*_result++ = *_curr++;
	if (_curr != _eor) std::push_heap(_first, _last, _comp);
	else --_last;
    }
}

using namespace std;
using namespace magicbox;

//--------------------------------------------------
// Main program
//-------------------------------------------------- 
int main(int argc, char** argv) {
    // Getopt
    unsigned int max_entry = 1024 * 1024;
    string eol_token = "#eol";
    bool drop_first_token = false;

    Getopt g(argc, argv);
    g   << $(&max_entry, "max-entry,M", "The maximum number of vocabulary stored in memory")
	<< $(&eol_token, "eol-token", "The special token indicating the end of line")
	<< $(&drop_first_token, "drop-first-token", 
		"Drop the first token in each line.  Useful when the first token is used as docno.")
	<< $$$("[options..]");

    typedef unordered_map<string, unsigned int> freq_map;
    typedef unordered_set<string> term_set;

    namespace fs = boost::filesystem;

    // Go!
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
	terms.insert(eol_token);

	// Put those counts back in
	foreach (const string& term, terms) { ++df[term]; }

	// Flush as necesary
	if (df.size() > max_entry) {
	    string runfile = boost::str(boost::format("collect-%d.run.%05d") % pid % ++runno);
	    fs::ofstream runout(runfile);

	    cerr << "collect: flush to " << runfile << '\n';

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
	cerr << "collect: flush to " << runfile << '\n';

	vector<string> terms;
	foreach (const freq_map::value_type &p, df) { terms.push_back(p.first); }
	stable_sort(terms.begin(), terms.end());
	foreach (const string &term, terms) { runout << term << ' ' << df[term] << '\n'; }
	df.clear(); // Done
    }

    //--------------------------------------------------
    // Step 2: Merge disk-runs
    //-------------------------------------------------- 
    
    string mergefile = boost::str(boost::format("collect-%d.merge") % pid);

    // Prepare the merge output
    //
    // NOTE: This step can be integrated into the multiway framework as a `reduce' function.
    //       But I don't want to deal with that right now.

    {
	fs::ofstream mergeout(mergefile);

	// Put every runfile back into the pool, as istreams
	ifstream_pool p;
	for (unsigned int i = 1; i <= runno; ++i) 
	    p.push_back(new ifstream(boost::str(boost::format("collect-%d.run.%05d") % pid % i).c_str()));

	cerr << "collect: produce the merge run" << '\n';

	ifstream_run end;
	vector<ifstream_run> runs;
	foreach (ifstream& in, p) { runs.push_back(in); }
	multiway_merge(runs.begin(), runs.end(), 
		ostream_iterator<record>(mergeout, "\n"), greater<ifstream_run>());
    }

    // Produce the result
    
    {
	fs::ifstream mergein(mergefile);
	record r_prev, r;

	cerr << "collect: output" << '\n';

	if (mergein >> r_prev) {
	    while (mergein >> r) {
		if (r.key == r_prev.key) r_prev.count += r.count;
		else {
		    cout << r_prev << '\n';
		    r_prev = r;
		}
	    }
	    cout << r_prev << '\n';
	}
    }

    cerr << "collect: done" << '\n';
    return 0;
}
