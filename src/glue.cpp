#include <iostream>
#include <boost/foreach.hpp>
#include "kinkaseki/CLI.h"

#define foreach BOOST_FOREACH

#include "compat.h"
#include <iterator>
#include <boost/functional/hash.hpp>
#include <algorithm>
#include <ext/functional>

#include "kinkaseki/Lexicon.h"
#include "kinkaseki/LineReader.h"

//--------------------------------------------------
// Helper functions
//-------------------------------------------------- 
template<typename Pair, typename Predicate>
struct Choose1st {
    Predicate pred;

    bool operator()(const Pair& x, const Pair& y) {
	return pred(x.first, y.first);
    }
};

template<typename Pair, typename Predicate>
Choose1st<Pair, Predicate> choose1st(Predicate pred) {
    return Choose1st<Pair, Predicate>();
}

template<typename Pair, typename Predicate>
struct Choose2nd {
    Predicate pred;

    bool operator()(const Pair& x, const Pair& y) {
	return pred(x.second, y.second);
    }
};

template<typename Pair, typename Predicate>
Choose2nd<Pair, Predicate> choose2nd(Predicate pred) {
    return Choose2nd<Pair, Predicate>();
}

//--------------------------------------------------
// Main program goes here
//-------------------------------------------------- 
int main(int argc, char** argv) {
    using namespace std;

    kinkaseki::CLI cli(argc, argv);

    int numIteration = 100000;
    int topK = 5;
    int minSupport = 3;

    cli
	.bind("verbose", "Show verbose output")
	.bind(numIteration, "iteration,i", "Specify the number of iteration")
	.bind(topK, "top,t", "Process the top K bigrams per iteration")
	.bind(minSupport, "support,s", "Specify the minimum support")
	.setSynopsis("Segment input texts using the glue algorithm\n")
	.setTexts(
	    "  No concrete example so far.\n"
	);

    vector<string> args = cli.parse();

    // Step 1: Read input and store them as a token stream
    kinkaseki::Lexicon lexicon;
    const int UNK = lexicon.encode("");
    const int EOL = lexicon.encode("\n");

    vector<int> text;
    text.push_back(EOL); // The first token is an <eol>

    kinkaseki::LineReader reader;
    while (istream& line = reader.getline(cin)) {
	string token;
	while (line >> token) 
	    text.push_back(lexicon.encode(token));
	text.push_back(EOL);
    }

    // Step 2: Create posting lists
    //
    // We assume the size of the text stream fits into a 4-byte integer.
    // From now on, we'll call each token as a 'Unigram'.
    typedef int Unigram;
    typedef pair<int, int> Bigram;
    typedef vector<unsigned int> PostingList;

    // NOTE: We keep track of the positions for each unigram (and thus we know
    // the frequencies).  For bigram, we only save the counts.
    typedef vector<PostingList> UnigramIndex;
    typedef unordered_map<Bigram, int, boost::hash<Bigram> > BigramIndex;

    UnigramIndex unigram;
    unigram.reserve(lexicon.size());

    BigramIndex bigram(lexicon.size());

    {
	typedef vector<Unigram>::iterator iterator;

	// Note the range: [begin, end - 1)
	iterator first = text.begin(), last = text.end() - 1;
	iterator iter = first + 1;

	Unigram prev = *first;
	if (prev != EOL) unigram[prev].push_back(0);

	while (iter != last) {
	    Unigram curr = *iter;
	    if (curr != EOL) {
		unigram[curr].push_back(distance(first, iter));
		if (prev != EOL) bigram[Bigram(prev, curr)]++;
	    }

	    prev = curr;
	    ++iter;
	}
    }

    // FIXME: Calculate the entropy and the rest

    // Now, enter the loop
    int iteration = 0;
    while (++iteration <= numIteration) {

	if (iteration % 1000 == 0)
	    cerr << "Iteration " << iteration << "\n";

	// Step 3: Populate the top-k bigrams
	//
	// More detail later
	typedef pair<Bigram, float> BigramScore;

	vector<BigramScore> topBigram;

	{
	    vector<BigramScore> score;
	    score.reserve(bigram.size());

	    BigramIndex::iterator iter = bigram.begin(), last = bigram.end();
	    while (iter != last) { 
		Unigram x = iter->first.first, y = iter->first.second;

		int f_xy = iter->second;
		int f_x = unigram[x].size();
		int f_y = unigram[y].size();

		if (x != y && f_xy >= minSupport) {
		    using std::log;
		    // float g = std::log(f_x) + std::log(f_y) - std::log(f_xy);
		    float g = f_x * log(f_x) + f_y * log(f_y) - f_xy * log(f_xy);
		    if (f_x > f_xy) g -= (f_x - f_xy) * log(f_x - f_xy);
		    if (f_y > f_xy) g -= (f_y - f_xy) * log(f_y - f_xy);

		    score.push_back(BigramScore(iter->first, g));
		}

		++iter;
	    }

	    partial_sort(
		score.begin(), 
		score.begin() + topK, 
		score.end(), 
		choose2nd<BigramScore>(std::less<float>())
	    );

	    copy(
		score.begin(), 
		score.begin() + topK, 
		back_inserter(topBigram)
	    );
	}

	// Step 4: Scan-Rewrite procedure
	//
	// More detail later
	unsigned int maxPos = text.size();

	foreach (const BigramScore& bs, topBigram) {
	    int x = bs.first.first;
	    int y = bs.first.second;
	    float score = bs.second;

	    if (cli["verbose"]) 
		cerr << lexicon.decode(x) << lexicon.decode(y) << ' ' 
		//--------------------------------------------------
		// << bigram[Bigram(x, y)] << ' ' 
		// << unigram[x].size() << ' '
		// << unigram[y].size() << ' '
		//-------------------------------------------------- 
		<< score << "\n";

	    // (1) Prepare the posting lists for x, y, and xy
	    PostingList::iterator 
		xiter = unigram[x].begin(), xlast = unigram[x].end();
	    PostingList::iterator 
		yiter = unigram[y].begin(), ylast = unigram[y].end();
	    PostingList pl_x, pl_y, pl_xy;

	    while (xiter != xlast && yiter != ylast) {
		if (*xiter < *yiter) {
		    unsigned int nextPos = *xiter + 1;
		    while (nextPos < maxPos && text[nextPos] == UNK) ++nextPos;

		    if (nextPos == *yiter) {
			pl_xy.push_back(*xiter);
			++xiter;
			++yiter;
		    }
		    else
			pl_x.push_back(*xiter++);
		}
		else
		    pl_y.push_back(*yiter++);
	    }

	    if (xiter != xlast)
		pl_x.insert(pl_x.end(), xiter, xlast); // oops

	    if (yiter != ylast)
		pl_y.insert(pl_y.end(), yiter, ylast); // oops

	    // (2) Prepare the update, the decrement, and the increment lists
	    unordered_set<Bigram, boost::hash<Bigram> > update;
	    unordered_set<unsigned int> decrement;
	    unordered_set<unsigned int> increment;

	    // (3) Mark immediate neighbors of ``xy'' and correct the counts
	    // i.e., the decrement list
	    foreach (unsigned int pos, pl_xy) {
		unsigned int prevPos = pos - 1, nextPos = pos + 1;
		while (prevPos >= 0 && text[prevPos] == UNK) --prevPos;
		while (nextPos < maxPos && text[nextPos] == UNK) ++nextPos;

		unsigned int nextPos2 = nextPos + 1;
		while (nextPos2 < maxPos && text[nextPos2] == UNK) ++nextPos2;

		if (prevPos >= 0 && text[prevPos] != 1) 
		    decrement.insert(prevPos); // do ``ax'' if a exists

		decrement.insert(pos); // always do ``xy''

		if (nextPos2 < maxPos && text[nextPos2] != 1) 
		    decrement.insert(nextPos); // do ``yb'' if b exists
	    }

	    foreach (unsigned int pos, decrement) {
		unsigned int nextPos = pos;
		while (text[++nextPos] == UNK) ; // just being slack

		Bigram b(text[pos], text[nextPos]);
		update.insert(b);
		bigram[b]--;
	    }

	    // (4) Rewrite ``xy'' as ``z0'' (0 as the padding symbol)
	    
	    // FIXME: Check if this is consistent to the lexicon
	    int z = lexicon.encode(lexicon.decode(x) + lexicon.decode(y));

	    foreach (unsigned int pos, pl_xy) {
		unsigned int nextPos = pos;
		while (text[++nextPos] == UNK) ; // being slack again

		text[pos] = z;
		text[nextPos] = UNK;

		unsigned int prevPos = pos - 1, nextPos2 = nextPos + 1;
		while (prevPos >= 0 && text[prevPos] == UNK) --prevPos;
		while (nextPos2 < maxPos && text[nextPos2] == UNK) ++nextPos2;

		if (prevPos >= 0 && text[prevPos] != 1) 
		    increment.insert(prevPos); // do ``az'' if a exists

		if (nextPos2 < maxPos && text[nextPos2] != 1)
		    increment.insert(pos); // do ``zb'' if b exists
	    }

	    foreach (unsigned int pos, increment) {
		unsigned int nextPos = pos;
		while (text[++nextPos] == UNK) ; // just being slack here

		Bigram b(text[pos], text[nextPos]);
		update.insert(b);
		bigram[b]++;
	    }

	    // (5) Renew the posting lists
	    unigram.push_back(PostingList());
	    unigram[x] = pl_x;
	    unigram[y] = pl_y;
	    unigram[z] = pl_xy;
	}
    }

    // Step 6: Output
    foreach (Unigram u, text) {
	cout << lexicon.decode(u) << ' ';
    }

    return 0;
}
