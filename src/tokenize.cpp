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
// PorterStemmer
//-------------------------------------------------- 
extern "C" {
#include "porter.h"
}

class PorterStemmer {
    struct stemmer* _stem;

public:
    PorterStemmer() { _stem = create_stemmer(); }
    ~PorterStemmer() { free_stemmer(_stem); }
    std::string stem(const char* p, int k) {
	if (k >= 1024) return std::string(); // Skip this term

	char buf[1024];
	buf[k + 1] = 0;
	for (int i = 0; i <= k; i++) 
	    buf[i] = (p[i] >= 'A' && p[i] <= 'Z')? p[i] - 'A' + 'a': p[i];

	int kk = ::stem(_stem, buf, k); 
	return std::string(buf, kk + 1);
    }

    std::string stem(const std::string& s) {
	return stem(s.c_str(), s.size() - 1);
    }
};

//--------------------------------------------------
// Main program
//-------------------------------------------------- 
int main(int argc, char** argv) {

    // Getopt
    Getopt g(argc, argv);
    g   << $("no-act,n", "Do nothing")
	<< $("stem", "Explicitly lemmatize all the tokens")
	<< $("utf8", "Use cjkv_separator instead")
	<< $$$("[options..]");

    // Bypass as necessary
    if (g["no-act"]) {
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
	PorterStemmer ps;
	string line;

	bool do_stem = g["stem"];
	bool do_utf8 = g["utf8"];

	while (getline(cin, line)) {
	    if (do_utf8) {
		boost::to_lower(line);
		CJKVTokenizer tok(line);
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
		boost::to_lower(line);
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
