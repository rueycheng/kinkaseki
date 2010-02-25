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

typedef vector< vector<unsigned int> > forward_table;

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

    // Reguar bookkeeping stuff
    vector<string> docnos;
    vector<unsigned int> dlen;
    inverted_table inv;
    freq_map tf;
    freq_map df;
    unsigned int N = 0;

    // For facets
    vector<string> facets;
    forward_table fwd;
    freq_map ftab;

    foreach (const string& filename, input) {
	AutoIn in(filename);
	string line;

	while (getline(in(), line)) {
	    istringstream iss(line);
	    string docno;
	    iss >> docno;
	    cerr << docno << '\n';

	    if (!boost::ends_with(docno, ":facet")) {
		// Case 1: Regular texts
		docnos.push_back(docno);
		dlen.push_back(0);
		unsigned int docid = docnos.size();

		{
		    // Get counts
		    freq_map counts;
		    string word;
		    while (iss >> word) ++counts[word];

		    foreach (const freq_map::value_type& p, counts) {
			++df[p.first];
			tf[p.first] += p.second;
			inv[p.first].push_back(posting(docid, p.second));
			N += p.second;
			dlen.back() += p.second;
		    }
		}
	    }
	    else {
		// Case 2: Facets
		fwd.push_back(vector<unsigned int>());

		{
		    string word;
		    while (iss >> word) {
			if (ftab.find(word) == ftab.end()) {
			    facets.push_back(word);
			    ftab[word] = facets.size();
			}

			fwd.back().push_back(ftab[word]);
		    }
		}
	    }
	}
    }

    cerr << "Model: " << model << '\n';

    // Output inv and vocab
    {
	// Sort the keys and output
	vector<string> vocab;
	foreach (const inverted_table::value_type& entry, inv) { vocab.push_back(entry.first); }
	stable_sort(vocab.begin(), vocab.end());

	fs::ofstream vocab_out(basedir / "vocab");
	fs::ofstream inv_out(basedir / "inv");
	cerr << "Save vocab and inv" << '\n';

	vocab_out << vocab.size() << '\n';
	inv_out << N << '\n';

	foreach (const string& v, vocab) { 
	    vocab_out << v << '\n';
	    inv_out << df[v] << ' ' << tf[v] << ' ' << inv[v] << '\n'; 
	}
    }

    // Output docno
    {
	fs::ofstream docno_out(basedir / "docno");
	fs::ofstream dlen_out(basedir / "dlen");
	cerr << "Save docno and dlen" << '\n';

	docno_out << docnos.size() << '\n';
	dlen_out << dlen.size() << '\n';
	foreach (const string& d, docnos) { docno_out << d << '\n'; }
	foreach (const unsigned int len, dlen) { dlen_out << len << '\n'; }
    }

    // Output facet
    {
	fs::ofstream facet_out(basedir / "facet");
	fs::ofstream ffwd_out(basedir / "ffwd");
	cerr << "Save facet and ffwd" << '\n';

	facet_out << facets.size() << '\n';
	ffwd_out << fwd.size() << '\n';
	foreach (const string& f, facets) { facet_out << f << '\n'; }
	foreach (const forward_table::value_type& flist, fwd) {
	    ffwd_out << flist.size() << ' ';
	    foreach (const unsigned int id, flist) { ffwd_out << id << ' '; }
	    ffwd_out << '\n';
	}
    }

    return 0;
}
