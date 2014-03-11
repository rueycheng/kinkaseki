#include <algorithm>
#include <iostream>
#include <fstream>
#include <iterator>
#include <boost/format.hpp>
#include <boost/foreach.hpp>
#include <boost/functional/hash.hpp>

#include "compat.h"
#include "kinkaseki/functors.hpp"
#include "kinkaseki/CLI.hpp"
#include "kinkaseki/Lexicon.hpp"
#include "kinkaseki/LineReader.hpp"

typedef int Unigram;
typedef std::pair<int, int> Bigram;

void 
output_segmented_text(std::ostream& out, std::vector<Unigram>& text, 
		      kinkaseki::Lexicon& lexicon, Unigram unk, Unigram eol)
{
    bool newline = true;
    BOOST_FOREACH (Unigram u, text) {
	if (u == unk) continue;

	if (u == eol) {
	    out << "\n";
	    newline = true;
	}
	else {
	    if (newline)
		newline = false;
	    else
		out << ' ';
	    out << lexicon.decode(u);
	}
    }
}

//--------------------------------------------------
// Main program goes here
//-------------------------------------------------- 
int main(int argc, char** argv) {
    using namespace std;

    kinkaseki::CLI cli(argc, argv);

    float alpha = 1.0;
    float beta = 1.0;
    float dir = 0.0;
    int numIteration = 10000;
    double ratio = 0.0;
    int topK = 1;
    // int minSupport = 3;
    // int charLimit = 10000;
    // int subcharLimit = 2;
    // int objectiveType = 1;
    // string punct = "";
    int outputEvery = 100;
    string outputPath = "";

    cli
	.bind("verbose,v", "Show verbose output")
	.bind(alpha, "alpha,a", "Specify parameter 'alpha'")
	.bind(beta, "beta,b", "Specify parameter 'beta'")
	.bind(dir, "dir,d", "Specify parameter 'dir'")
	.bind(numIteration, "iteration,i", "Specify the number of iteration")
	.bind(ratio, "ratio,r", "Specify the expected word/token ratio as the terminal condition")
	.bind(topK, "top,t", "Process the top K bigrams per iteration")
	.bind(outputEvery, "outputEvery", "Output temporary result every N iterations")
	.bind(outputPath, "outputPath", "Output result to this directory")
	// .bind(minSupport, "support,s", "Specify the minimum support")
	// .bind(charLimit, "limit", "Specify the maximum number of characters in a word")
	// .bind(subcharLimit, "sublimit", "Specify the maximum number of characters in a subword")
	// .bind(objectiveType, "type", "Specify the objective function (1=ratio, 2=freq, 3=logfreq)")
	// .bind(punct, "punct", "Specify the punctuation marks, separated by whitespace")
	.setSynopsis("Segment input texts using the glue algorithm\n")
	.setTexts(
	    "  No concrete example so far.\n"
	);

    vector<string> args = cli.parse();

    // Step 1: Read input and store them as a token stream
    //
    // We assume the size of the text stream fits into a 4-byte integer.
    // From now on, we'll call each token as a 'Unigram'.

    kinkaseki::Lexicon lexicon;
    const Unigram UNK = lexicon.encode("");
    const Unigram EOL = lexicon.encode("\n");

    lexicon.encode("\xef\xbc\x81");
    lexicon.encode("\xef\xbc\x8c");
    lexicon.encode("\xef\xbc\x88");
    lexicon.encode("\xef\xbc\x89");
    lexicon.encode("\xef\xbc\x9a");
    lexicon.encode("\xef\xbc\x9f");
    lexicon.encode("\xef\xb8\xb0");
    lexicon.encode("\xe3\x80\x81");
    lexicon.encode("\xe3\x80\x82");
    lexicon.encode("\xe3\x80\x8c");
    const Unigram SPECIAL_END = 
	lexicon.encode("\xe3\x80\x8d"); // Mark the last special token

    vector<Unigram> text;
    kinkaseki::LineReader reader;

    while (istream& line = reader.getline(cin)) {
	string token;
	while (line >> token) 
	    text.push_back(lexicon.encode(token));
	text.push_back(EOL);
    }

    // Step 2: Create posting lists
    //
    // NOTE: We keep track of the positions for each unigram (and thus we know
    // the frequencies).  For bigram, we only save the counts.
    typedef vector<unsigned int> PostingList;
    typedef vector<PostingList> UnigramIndex;
    typedef vector<int> UnigramSizeIndex;
    typedef unordered_map<Bigram, int, boost::hash<Bigram> > BigramIndex;

    // Set up 
    UnigramIndex unigram(lexicon.size(), PostingList());
    UnigramSizeIndex unigramSize(lexicon.size(), 1); // initially 1
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
		if (prev != EOL)
		    bigram[Bigram(prev, curr)]++;
	    }

	    prev = curr;
	    ++iter;
	}
    }

    int numTokens = 0;
    BOOST_FOREACH (const PostingList& pl, unigram) {
	numTokens += pl.size();
    }

    // Now, enter the loop
    int iteration = 0;
    while (++iteration) {
	// N denotes the total number of `regular' tokens
	// H denotes H(W), J denotes J(W)
	//
	// NOTE: W stands for current support (more details later)
	unsigned int N = 0;
	// float H = 0.0, J = 0.0;

	{
	    typedef vector<PostingList>::iterator iterator;

	    iterator iter = unigram.begin(), last = unigram.end();
	    while (iter != last) {
		int size = iter->size();
		if (size > 0) {
		    // J += size * std::log(size);
		    N += size;
		}

		++iter;
	    }

	    // H = std::log(N) - J / N;

	    // Show the current entropy
	    // cerr << "N " << N << " J " << J << " H " << H << "\n";
	}

	double currentRatio = static_cast<double>(N) / numTokens;

	if (currentRatio < ratio || (ratio <= 0.0 && iteration > numIteration))
	    break;

	// if (cli["verbose"])
	    // cerr << iteration << " " << currentRatio << "\n";

	if (iteration % 100 == 0)
	    cerr << "Iteration " << iteration << " (ratio: " << currentRatio << ")\n";

	// Step 3: Populate the top-k bigrams
	//
	// More detail later
	typedef pair<Bigram, float> BigramScore;

	vector<BigramScore> topBigram;

	{
	    using std::log;

	    vector<BigramScore> score;
	    score.reserve(bigram.size());

	    float log_N = std::log(N);
	    // float average_J = J / N;

	    BigramIndex::iterator iter = bigram.begin(), last = bigram.end();
	    for (; iter != last; ++iter) { 
		Unigram x = iter->first.first;
		Unigram y = iter->first.second;

		int f_xy = iter->second + dir;
		int f_x = unigram[x].size() + dir;
		int f_y = unigram[y].size() + dir;

		if (f_xy < 1) continue;

		if (x <= SPECIAL_END) continue;
		if (y <= SPECIAL_END) continue;
		//--------------------------------------------------
		// if (unigramSize[x] + unigramSize[y] > charLimit) continue;
		//
		// int smaller = unigramSize[x] > unigramSize[y]? 
		//     unigramSize[y]: unigramSize[x];
		// if (smaller > subcharLimit) continue;
		//
		// if (f_xy < minSupport) continue;
		//-------------------------------------------------- 
		
		// float delta_H = log((f_x - f_xy) * (f_y - f_xy) / f_xy) - log_N;
		// float objective = objectiveFunction(beta, f_xy, N, numTokens, delta_H, unigram.size());
		int unigramSize_xy = unigramSize[x] + unigramSize[y];

		float objective = 
		    f_xy * (1 + log(f_x * f_y / f_xy) - log_N - alpha) + 
		    beta * f_xy * (unigramSize_xy * log(unigramSize_xy) - 
				   unigramSize[x] * log(unigramSize[x]) - unigramSize[y] * log(unigramSize[y]));
		    
		score.push_back(BigramScore(iter->first, objective));
	    }

	    if (score.size() > static_cast<unsigned int>(topK)) {
		partial_sort(score.begin(), score.begin() + topK, 
			     score.end(), kinkaseki::choose2nd<BigramScore>(std::less<float>()));

		copy(score.begin(), score.begin() + topK, back_inserter(topBigram));
	    }
	    else {
		copy(score.begin(), score.end(), back_inserter(topBigram));
	    }
	}

	if (topBigram.empty()) break;

	// Step 4: Scan-Rewrite procedure
	//
	// More detail later
	unsigned int maxPos = text.size();

	BOOST_FOREACH (const BigramScore& bs, topBigram) {
	    int x = bs.first.first;
	    int y = bs.first.second;
	    float score = bs.second;

	    if (cli["verbose"]) 
		cerr 
		    << iteration << ' ' 
		    << currentRatio << ' '
		    << lexicon.decode(x) << lexicon.decode(y) << ' ' 
		    << unigram[x].size() << ' '
		    << unigram[y].size() << ' '
		    << bigram[bs.first] << ' '
		    << score << "\n";

	    // (1) Prepare and renew the posting lists for x, y, and xy (i.e., z)
	    int z = lexicon.encode(lexicon.decode(x) + lexicon.decode(y));
	    unigram.push_back(PostingList()); // Now unigram[z] exists
	    unigramSize.push_back(unigramSize[x] + unigramSize[y]);

	    if (x != y) {
		PostingList pl_x, pl_y, pl_xy;

		PostingList::iterator 
		    xiter = unigram[x].begin(), xlast = unigram[x].end();
		PostingList::iterator 
		    yiter = unigram[y].begin(), ylast = unigram[y].end();

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

		unigram[x] = pl_x;
		unigram[y] = pl_y;
		unigram[z] = pl_xy;
	    }
	    else {
		PostingList pl_x, pl_xx;

		// FIXME: Assert that unigram[x] contains at least one elements
		PostingList::iterator 
		    xiter = unigram[x].begin(), xlast = unigram[x].end();

		while (xiter != xlast - 1) { // NOTE: [begin, last - 1)
		    unsigned int nextPos = *xiter + 1;
		    while (nextPos < maxPos && text[nextPos] == UNK) ++nextPos;

		    if (*(xiter + 1) == nextPos) {
			pl_xx.push_back(*xiter);
			xiter += 2;

			if (xiter == xlast) break;
		    }
		    else
			pl_x.push_back(*xiter++);
		}

		pl_x.insert(pl_x.end(), xiter, xlast); // Always keep the last ones

		unigram[x] = pl_x;
		unigram[z] = pl_xx;
	    }

	    // (2) Prepare the update, the decrement, and the increment lists
	    unordered_set<Bigram, boost::hash<Bigram> > update;
	    unordered_set<unsigned int> decrement;
	    unordered_set<unsigned int> increment;

	    // (3) Mark immediate neighbors of ``xy'' and correct the counts
	    // i.e., the decrement list
	    BOOST_FOREACH (unsigned int pos, unigram[z]) {
		int prevPos = pos - 1; 
		unsigned int nextPos = pos + 1;
		while (prevPos >= 0 && text[prevPos] == UNK) --prevPos;
		while (nextPos < maxPos && text[nextPos] == UNK) ++nextPos;

		unsigned int nextPos2 = nextPos + 1;
		while (nextPos2 < maxPos && text[nextPos2] == UNK) ++nextPos2;

		if (prevPos >= 0 && text[prevPos] != EOL) 
		    decrement.insert(prevPos); // do ``ax'' if a exists

		decrement.insert(pos); // always do ``xy''

		if (nextPos2 < maxPos && text[nextPos2] != EOL) 
		    decrement.insert(nextPos); // do ``yb'' if b exists
	    }

	    BOOST_FOREACH (unsigned int pos, decrement) {
		unsigned int nextPos = pos;
		while (text[++nextPos] == UNK) ; // just being slack

		Bigram b(text[pos], text[nextPos]);
		update.insert(b);
		bigram[b]--;
	    }

	    // (4) Rewrite ``xy'' as ``z0'' (0 as the padding symbol)
	    
	    BOOST_FOREACH (unsigned int pos, unigram[z]) {
		unsigned int nextPos = pos;
		while (text[++nextPos] == UNK) ; // being slack again

		text[pos] = z;
		text[nextPos] = UNK;

		int prevPos = pos - 1; 
		unsigned int nextPos2 = nextPos + 1;
		while (prevPos >= 0 && text[prevPos] == UNK) --prevPos;
		while (nextPos2 < maxPos && text[nextPos2] == UNK) ++nextPos2;

		if (prevPos >= 0 && text[prevPos] != 1) 
		    increment.insert(prevPos); // do ``az'' if a exists

		if (nextPos2 < maxPos && text[nextPos2] != 1)
		    increment.insert(pos); // do ``zb'' if b exists
	    }

	    BOOST_FOREACH (unsigned int pos, increment) {
		unsigned int nextPos = pos;
		while (text[++nextPos] == UNK) ; // just being slack here

		Bigram b(text[pos], text[nextPos]);
		update.insert(b);
		bigram[b]++;
	    }
	}

	if (iteration % outputEvery == 0) {
	    string filename = boost::str(boost::format("output.%05d") % iteration);

	    cerr << boost::format("Save output to %s") % filename << '\n';
	    ofstream fout(filename);
	    output_segmented_text(fout, text, lexicon, UNK, EOL);
	}
    }

    // Step 6: Output
    output_segmented_text(cout, text, lexicon, UNK, EOL);
    //--------------------------------------------------
    // BOOST_FOREACH (Unigram u, text) {
    // 	if (u == UNK) continue;
    // 
    // 	if (u == EOL) 
    // 	    cout << "\n";
    // 	else
    // 	    cout << lexicon.decode(u) << ' ';
    // }
    //-------------------------------------------------- 

    return 0;
}
