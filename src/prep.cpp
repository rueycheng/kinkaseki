#include <iostream>
#include <sstream>

#include "magicbox.h"
#include "common.h"
#include "compat.h"

using namespace std;
using namespace magicbox;

int main(int argc, char** argv) {

    // Getopt
    vector<string> input;

    Getopt g(argc, argv);
    g   << $("df", "Collect document-frequency")
	<< $("tf", "Collect term-frequency")
	<< $(&input, "input", "", -1)
	<< $$$("files..");

    // Main program
    unordered_map<string,unsigned int> df;

    string line, word;
    while (getline(cin, line)) {
	istringstream in(line);
	while (in >> word) ++df[word];
    }

    cout << df.size() << '\n';

    return 0;
}
