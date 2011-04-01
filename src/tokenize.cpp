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
#include "kinkaseki/StopList.h"

using namespace std;

template<typename Separator, 
    typename StopList, 
    typename Stemmer>
class Pipeline {
    Separator separator;
    StopList stoplist;
    Stemmer stemmer;

public:
    Pipeline(Separator separator = Separator(), 
	     StopList stoplist = StopList(),
	     Stemmer stemmer = Stemmer()):
	separator(separator), stoplist(stoplist), stemmer(stemmer) 
    {}

    void run(istream& in, ostream& out) {
	string line;

	while (getline(in, line)) {
	    boost::to_lower(line);
	    boost::tokenizer<Separator> tok(line, separator);
	    foreach (const string& t, tok) { 
		if (stoplist.has(t)) continue;
		out << stemmer.stem(t) << ' ';
	    }
	    out << "\n";
	}
    }
};

template<typename Separator, typename StopList, typename Stemmer>
Pipeline<Separator, StopList, Stemmer> 
make_pipeline(Separator& separator, StopList& stoplist, Stemmer& stemmer) {
    return Pipeline<Separator, StopList, Stemmer>(separator, stoplist, stemmer);
}

//--------------------------------------------------
// Main program
//-------------------------------------------------- 
int main(int argc, char** argv) {

    kinkaseki::CLI cli(argc, argv);

    cli
	.bind("bypass", "Bypass")
	.bind("no-stem", "Do not employ stemmer")
	.bind("no-utf8", "Use default separator.  Do not use CJKV-enabled separator")
	.setSynopsis("Tokenize the input");

    cli.parse();

    // Bypass as necessary
    if (cli["bypass"]) {
	cout << cin.rdbuf();
	return 0;
    }

    // No stem
    if (cli["no-stem"]) {
	cerr << "Option --no-stem is currently unsupported." << endl;
	return 0;
    }

    // Stopword list
    const char* stopwords[] = {
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
	"with", "would", "yet", "you", "your" };

    kinkaseki::StopList stoplist(stopwords);
    kinkaseki::PorterStemmer porter;

    if (cli["no-utf8"]) {
	boost::char_separator<char> sep;
	make_pipeline(sep, stoplist, porter).run(cin, cout);
    }
    else {
	kinkaseki::CJKVSeparator sep;
	make_pipeline(sep, stoplist, porter).run(cin, cout);
    }

    return 0;
}
