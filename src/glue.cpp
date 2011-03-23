#include <iostream>
#include <boost/foreach.hpp>
#include "compat.h"
#include "kinkaseki/CLI.h"

#define foreach BOOST_FOREACH

#include <iterator>
#include <boost/functional/hash.hpp>
#include <unordered_map>
#include <algorithm>

template<typename UnigramMap, typename BigramMap>
struct MoreEntropyGain {
    UnigramMap& umap;
    BigramMap& bmap;

    MoreEntropyGain(UnigramMap& um, BigramMap& bm): umap(um), bmap(bm) {}
    bool operator()(typename BigramMap::key_type x, typename BigramMap::key_type y) {
	unsigned int fx = bmap[x].size();
	unsigned int fy = bmap[y].size();

	return fx > fy;
    }
};

int main(int argc, char** argv) {
    using namespace std;
    using kinkaseki::CLI;

    CLI cli(argc, argv);

    cli
	.bind("verbose", "Show verbose output")
	.setSynopsis("Segment input texts using the glue algorithm\n")
	.setTexts(
	    "  No concrete example so far.\n"
	);

    vector<string> args = cli.parse();

    // Main program starts here
    //
    // Step 1: House keeping
    vector<unsigned int> text;
    unordered_map<string, unsigned int> lexicon;
    unsigned int nextID = 0;

    lexicon["<>"] = 0; // means 'empty'
    lexicon["<eol>"] = 1; // mean 'end-of-line'

    text.push_back(1); // The first token is an <eol>

    string line;
    while (getline(cin, line)) {
	istringstream sin(line);
	string token;
	unsigned int tokenID;

	while (sin >> token) {
	    // FIXME: Optimize this line
	    if (lexicon.find(token) == lexicon.end()) 
		lexicon.insert(make_pair(token, tokenID = nextID++));
	    else
		tokenID = lexicon[token];

	    text.push_back(tokenID);
	}

	text.push_back(1); // Of course there is an <eol>
    }

    // Step 2: Create posting lists
    //
    // We assume the size of the text stream fits into a 4-byte integer
    // From now on, we'll call each token as a 'Unigram'
    typedef unsigned int Unigram;
    typedef pair<unsigned int, unsigned int> Bigram;
    typedef vector<unsigned int> PostingList;

    unordered_map<Unigram, PostingList, boost::hash<Unigram> > unigram;
    unordered_map<Bigram, PostingList, boost::hash<Bigram> > bigram;

    {
	vector<unsigned int>::iterator iter = text.begin(), first = text.begin();
	vector<unsigned int>::iterator last = text.end();

	++iter; // Go one step ahead

	while (iter != last) {
	    Unigram u = *iter;
	    Bigram b = Bigram(*(iter - 1), *iter);

	    unigram[u].push_back(distance(first, iter));
	    bigram[b].push_back(distance(first, iter - 1));
	    ++iter;
	}
    }

    // Step 3: Create a heap
    //
    // More detail later
    vector<Bigram> heap;
    heap.reserve(bigram.size());

    MoreEntropyGain<
	unordered_map<Unigram, PostingList, boost::hash<Unigram> >,
	unordered_map<Bigram, PostingList, boost::hash<Bigram> > 
	>
	compare(unigram, bigram);

    {
	unordered_map<Bigram, PostingList, boost::hash<Bigram> >
	    ::const_iterator iter = bigram.begin(), last = bigram.end();

	while (iter != last) heap.push_back(iter->first);

	make_heap(heap.begin(), heap.end(), compare);
    }

    while (!heap.empty()) {
	Bigram& top = heap.front();
	cout << bigram[top].size() << ' ' << top.first << ' ' << top.second << ' ' << "\n";
	pop_heap(heap.begin(), heap.end(), compare);
    }

    return 0;
}
