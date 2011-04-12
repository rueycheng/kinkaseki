#include <fstream>
#include <iostream>
#include <iterator>
#include <string>

#include "kinkaseki/CLI.h"

int main(int argc, char** argv) {
    using namespace std;

    kinkaseki::CLI cli(argc, argv);

    cli
	.setUsage("TRUTH-FILE TEST-FILE")
	.setSynopsis("Evaluate the precision/recall of the segmented text\n")
	.setTexts(
	    "  TRUTH-FILE\tThe gold standard segmented text\n"
	    "  TEST-FILE\tThe segmented text\n"
	);

    vector<string> args = cli.parse();

    if (args.size() != 2) {
	cli.showHelp();
	return 0;
    }

    ifstream truth(args[0].c_str());
    ifstream test(args[1].c_str());

    {
	int word_matched = 0, word_relevant = 0, word_retrieved = 0;
	int boundary_matched = 0, boundary_relevant = 0, boundary_retrieved = 0;

	string l0, l1;
	unsigned int lineno = 0;
	while (getline(truth, l0) && getline(test, l1)) {
	    ++lineno;

	    vector<int> b0, b1;
	    unsigned int pos = 0;

	    string::iterator iter0 = l0.begin();
	    string::iterator iter1 = l1.begin();

	    // Implicit boundary
	    b0.push_back(pos);
	    b1.push_back(pos);

	    while (iter0 != l0.end() && iter1 != l1.end()) {
		if (*iter0 == *iter1) {
		    if (isspace(*iter0)) {
			b0.push_back(pos);
			b1.push_back(pos);
		    }
		    else ++pos;

		    ++iter0;
		    ++iter1;
		}
		else if (isspace(*iter0)) {
		    b0.push_back(pos);
		    ++iter0;
		}
		else if (isspace(*iter1)) {
		    b1.push_back(pos);
		    ++iter1;
		}
		else 
		    break;
	    }

	    if (iter0 != l0.end() || iter1 != l1.end()) {
		cerr << "Content mismatched at line " << lineno << "\n";
		continue;
	    }

	    // Calculate word precision/recall
	    bool matched = false;
	    int w = 0, b = 0;
	    vector<int>::iterator i0 = b0.begin();
	    vector<int>::iterator i1 = b1.begin();

	    while (i0 != b0.end() && i1 != b1.end()) {
		if (*i0 == *i1) {
		    if (matched) ++w;
		    else matched = true;

		    ++b;
		    ++i0;
		    ++i1;
		}
		else {
		    matched = false;

		    if (*i0 < *i1) ++i0;
		    else ++i1;
		}
	    }

	    word_matched += w;
	    word_relevant += (b0.size() - 1);
	    word_retrieved += (b1.size() - 1);

	    // NOTE: trivial boundaries are discounted
	    boundary_matched += (b - 2);
	    boundary_relevant += (b0.size() - 2);
	    boundary_retrieved += (b1.size() - 2);
	}

	if (getline(truth, l0) || getline(test, l1)) {
	    cerr << "One of the files contains more lines than the other\n";
	    return 1;
	}

	// Report precision/recall
	float word_recall = float(word_matched) / word_relevant;
	float word_precision = float(word_matched) / word_retrieved;
	float word_f = 2 * word_recall * word_precision / (word_recall + word_precision);
	float boundary_recall = float(boundary_matched - 2) / boundary_relevant; 
	float boundary_precision = float(boundary_matched - 2) / boundary_retrieved;
	float boundary_f = 2 * boundary_recall * boundary_precision / (boundary_recall + boundary_precision);

	cout 
	    << "P " << word_precision 
	    << " R " << word_recall 
	    << " F " << word_f << "\n";
	cout 
	    << "BP " << boundary_precision 
	    << " BR " << boundary_recall 
	    << " BF " << boundary_f << "\n";
    }

    return 0;
}
