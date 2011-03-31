#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <cmath>
#include <numeric>

#include <boost/algorithm/string.hpp>
#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH

#include "compat.h"
#include "kinkaseki/CLI.h"
#include "kinkaseki/PorterStemmer.h"
#include "kinkaseki/CJKVSeparator.h"

using namespace std;

//--------------------------------------------------
// Main program
//-------------------------------------------------- 
int main(int argc, char** argv) {

    kinkaseki::CLI cli(argc, argv);

    bool do_stem = false;
    bool do_utf8 = false;

    cli
	.bind("no-act,n", "Do nothing")
	.bind(do_stem, "stem", "Explicitly lemmatize all tokens")
	.bind(do_utf8, "utf8", "Use CJKVSeparator instead")
	.setSynopsis("Tokenize the input");

    cli.parse();

    // Bypass as necessary
    if (cli["no-act"]) {
	cout << cin.rdbuf();
	return 0;
    }

    // Stopword list
    unordered_set<string> stopword;

    {
	const char* sw[] = {
	    "a", "able", "about", "across", "after", "all", "almost", "also", "am", "among", "an",
	    "and", "any", "are", "as", "at", "be", "because", "been", "but", "by", "can", "cannot",
	    "could", "dear", "did", "do", "does", "either", "else", "ever", "every", "for", "from",
	    "get", "got", "had", "has", "have", "he", "her", "hers", "him", "his", "how", "however",
	    "i", "if", "in", "into", "is", "it", "its", "just", "least", "let", "like", "likely",
	    "may", "me", "might", "most", "must", "my", "neither", "no", "nor", "not", "of", "off",
	    "often", "on", "only", "or", "other", "our", "own", "rather", "said", "say", "says",
	    "she", "should", "since", "so", "some", "than", "that", "the", "their", "them", "then",
	    "there", "these", "they", "this", "tis", "to", "too", "twas", "us", "wants", "was",
	    "we", "were", "what", "when", "where", "which", "while", "who", "whom", "why", "will",
	    "with", "would", "yet", "you", "your", 0 };

	const char** w = sw;
	while (*w != 0) stopword.insert(*w++);
    }

    {
	kinkaseki::PorterStemmer ps;
	string line;

	while (getline(cin, line)) {
	    boost::to_lower(line);

	    if (do_utf8) {
		// CJKVTokenizer tok(line);
		boost::tokenizer<kinkaseki::CJKVSeparator> tok(line);
		foreach (const string& t, tok) { 
 		    if (stopword.find(t) != stopword.end() || t.size() < 3 || t.size() > 25) continue;
		    if (!do_stem) cout << t << ' ';
		    else {
			if ((t[0] & 0x80) == 0x00) cout << ps.stem(t) << ' ';
			else cout << t << ' ';
		    }
		}
	    }
	    else {
		boost::char_separator<char> sep;
		boost::tokenizer<boost::char_separator<char> > tok(line, sep);
		foreach (const string& t, tok) { 
		    if (stopword.find(t) != stopword.end() || t.size() < 3 || t.size() > 25) continue;
		    if (!do_stem) cout << t << ' ';
		    else {
			if ((t[0] & 0x80) == 0x00) cout << ps.stem(t) << ' ';
			else cout << t << ' ';
		    }
		}
	    }
	    cout << "\n";
	}
    }

    return 0;
}
