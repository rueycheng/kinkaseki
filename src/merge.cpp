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

//--------------------------------------------------
// Main program
//-------------------------------------------------- 
int main(int argc, char** argv) {

    // Getopt
    vector<string> input;

    Getopt g(argc, argv);
    g   << $("bypass", "Do not tokenize or perform lemmatization")
	<< $(&input, "input", "", -1)
	<< $$$("files..");

    if (input.empty()) input.push_back("-");

    namespace fs = boost::filesystem;

    {
	string current_doc;

	foreach (const string& filename, input) {
	    AutoIn in(filename);
	    string line;

	    while (getline(in(), line)) {
		istringstream iss(line);
		string word;

		iss >> word;
		strip_after_first(word, ':');

		if (current_doc != word) {
		    if (!current_doc.empty()) cout << '\n';
		    current_doc = word;
		    cout << current_doc << ' ';
		}

		while (iss >> word) cout << word << ' ';
	    }
	    cout << '\n';
	}
    }

    return 0;
}
