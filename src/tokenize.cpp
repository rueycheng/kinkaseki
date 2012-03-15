#include <iostream>
#include <string>
#include <boost/algorithm/string.hpp>
#include <boost/foreach.hpp>

#include "kinkaseki/CLI.hpp"
#include "kinkaseki/PorterStemmer.hpp"
#include "kinkaseki/CJKVSeparator.hpp"
#include "kinkaseki/StopList.hpp"

template<typename Separator, typename StopList, typename Stemmer>
class Pipeline {
    Separator separator;
    StopList stoplist;
    Stemmer stemmer;

public:
    Pipeline(Separator separator = Separator(), 
	     StopList stoplist = StopList(),
	     Stemmer stemmer = Stemmer())
	:separator(separator), stoplist(stoplist), stemmer(stemmer) {}

    void run(std::istream& in, std::ostream& out) {
	std::string line;

	while (getline(in, line)) {
	    boost::to_lower(line);
	    boost::tokenizer<Separator> tok(line, separator);
	    BOOST_FOREACH (const std::string& t, tok) { 
		if (stoplist.match(t)) continue;
		out << stemmer.stem(t) << ' ';
	    }
	    out << "\n";
	}
    }
};

template<typename Separator, typename StopList, typename Stemmer>
Pipeline<Separator, StopList, Stemmer> 
pipeline(Separator separator, StopList stoplist, Stemmer stemmer) {
    return Pipeline<Separator, StopList, Stemmer>(separator, stoplist, stemmer);
}

//--------------------------------------------------
// Main program
//-------------------------------------------------- 
int main(int argc, char** argv) {
    using namespace std;

    kinkaseki::CLI cli(argc, argv);

    cli
	.bind("bypass", "Bypass")
	.bind("no-stoplist", "Do not filter stopwords")
	.bind("no-stemmer", "Do not employ stemmer")
	.bind("no-cjkv", "Use the default whitespace separator.")
	.setSynopsis("Tokenize the input");

    cli.parse();

    // Bypass as necessary
    if (cli["bypass"]) {
	cout << cin.rdbuf();
	return 0;
    }

    // Stopword list
    const char* stopwords[] = { 
	"a", "able", "about", "across", "after", "all", "almost", "also", "am",
	"among", "an", "and", "any", "are", "as", "at", "be", "because",
	"been", "but", "by", "can", "cannot", "could", "dear", "did", "do",
	"does", "either", "else", "ever", "every", "for", "from", "get", "got",
	"had", "has", "have", "he", "her", "hers", "him", "his", "how",
	"however", "i", "if", "in", "into", "is", "it", "its", "just", "least",
	"let", "like", "likely", "may", "me", "might", "most", "must", "my",
	"neither", "no", "nor", "not", "of", "off", "often", "on", "only",
	"or", "other", "our", "own", "rather", "said", "say", "says", "she",
	"should", "since", "so", "some", "than", "that", "the", "their",
	"them", "then", "there", "these", "they", "this", "tis", "to", "too",
	"twas", "us", "wants", "was", "we", "were", "what", "when", "where",
	"which", "while", "who", "whom", "why", "will", "with", "would", "yet",
	"you", "your" 
    };

    kinkaseki::StopList stoplist(stopwords, 
				 cli["no-stoplist"]? 0: 
				 sizeof(stopwords) / sizeof(const char*));

    kinkaseki::PorterStemmer stemmer(cli["no-stemmer"]);

    // FIXME: If-else looks really ugly
    if (cli["no-utf8"]) {
	boost::char_separator<char> sep;
	pipeline(sep, stoplist, stemmer).run(cin, cout);
    }
    else {
	kinkaseki::CJKVSeparator sep;
	pipeline(sep, stoplist, stemmer).run(cin, cout);
    }

    return 0;
}
