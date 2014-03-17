#include <algorithm>
#include <iostream>
#include <fstream>
#include <iterator>
#include <boost/filesystem.hpp>
#include <boost/format.hpp>
#include <boost/foreach.hpp>
#include <boost/functional/hash.hpp>
#include <boost/regex.hpp>

#include "compat.h"
#include "kinkaseki/functors.hpp"
#include "kinkaseki/CLI.hpp"
#include "kinkaseki/Lexicon.hpp"
#include "kinkaseki/LineReader.hpp"

typedef int Unigram;
typedef std::pair<int, int> Bigram;

void 
output_segmented_text(std::ostream& out, 
		      std::vector<Unigram>::const_iterator first,
		      std::vector<Unigram>::const_iterator last,
		      kinkaseki::Lexicon& lexicon, Unigram unk, Unigram eol)
{
    bool newline = true;

    for (; first != last; ++first) {
	if (*first == unk) continue;

	if (*first == eol) {
	    out << "\n";
	    newline = true;
	} else {
	    if (newline) newline = false;
	    else out << ' ';
	    out << lexicon.decode(*first);
	}
    }
}

void die(const std::string& message) {
    std::cerr << message << '\n';
    exit(1);
}

void die(const boost::format& fmt) {
    die(boost::str(fmt));
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
    double threshold = 0.0;
    double ratio = 0.0;
    int topK = 1;
    // int minSupport = 3;
    // int charLimit = 10000;
    // int subcharLimit = 2;
    // int objectiveType = 1;
    // string punct = "";
    int outputEvery = 100;
    string prefix = "output";
    string outputPath = ".";
    bool useX2Penalty = false;

    cli
	.bind("verbose,v", "Show verbose output")
	.bind("linguisticRule,l", "Prevent merging non-wordlike pairs")
	.bind(alpha, "alpha,a", "Specify parameter 'alpha'")
	.bind(beta, "beta,b", "Specify parameter 'beta'")
	.bind(dir, "dir,d", "Specify parameter 'dir'")
	.bind(numIteration, "iteration,I", "Specify the number of iteration")
	.bind(threshold, "threshold,T", "Stop when the gain is below this threshold")
	.bind(ratio, "ratio,r", "Specify the expected word/token ratio as the terminal condition")
	.bind(topK, "top,t", "Process the top K bigrams per iteration")
	.bind(outputEvery, "outputEvery", "Output temporary result every N iterations")
	.bind(prefix, "prefix", "Prefix for the output files")
	.bind(outputPath, "outputPath,D", "Output result to this directory")
	.bind(useX2Penalty, "useX2Penalty", "Use X^2 penalty function")
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

    if (!boost::filesystem::is_directory(outputPath) && !boost::filesystem::create_directory(outputPath))
	die(boost::format("Cannot create output path: '%s'") % outputPath);

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

    text.push_back(EOL);
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
	// iterator first = text.begin(), last = text.end() - 1;
	iterator first = text.begin(), last = text.end();
	iterator iter = first + 1;

	Unigram prev = *first;
	unigram[prev].push_back(0);

	while (iter != last) {
	    Unigram curr = *iter;
	    unigram[curr].push_back(distance(first, iter));
	    bigram[Bigram(prev, curr)]++;

	    prev = curr;
	    ++iter;
	}
    }

    int numTokens = 0;
    BOOST_FOREACH (const PostingList& pl, unigram) {
	numTokens += pl.size();
    }

    // Precompute xlogx values up to 100
    vector<double> xlogx, x2;
    xlogx.push_back(0); // 0 log 0 = 0
    x2.push_back(0); // 0^2 = 0
    for (int i = 1; i <= 100; ++i) {
	xlogx.push_back(i * log(i));
	x2.push_back(i * i);
    }

    // Set up regular expressions
    std::string vowel = "[#%&()*3679AEIOQRUaeilou~]";
    std::string consonant = "[DGLMNSTWZbcdfghklmnprstvwyz]";
    std::string syllable = consonant + "*" + vowel + "+" + consonant + "*";

    const boost::regex word_pattern("^(" + syllable + "){1,4}$");

    // Now, enter the loop
    bool quit = false;
    int iteration = 0;

    while (!quit && ++iteration) {
	unsigned int N = 0;

	{
	    typedef vector<PostingList>::iterator iterator;

	    iterator iter = unigram.begin(), last = unigram.end();
	    while (iter != last) 
		N += iter++->size();
	}

	double currentRatio = static_cast<double>(N) / numTokens;

	if (currentRatio < ratio || (ratio <= 0.0 && iteration > numIteration))
	    break;

	if (iteration % 100 == 0)
	    cerr << "Iteration " << iteration << " (ratio: " << currentRatio << ")\n";

	// Step 3: Populate the top-k bigrams
	//
	// More detail later
	typedef pair<Bigram, float> BigramScore;

	vector<BigramScore> topBigram;

	{
	    // using std::log2;
	    using std::log;

	    vector<BigramScore> score;
	    score.reserve(bigram.size());

	    float log_N = log(N);
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

		if (cli["linguisticRule"]) { 
		    string word_xy = lexicon.decode(x) + lexicon.decode(y);
		    if (!boost::regex_match(word_xy, word_pattern)) continue;
		}

		int W = unigram.size();

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

		if (static_cast<unsigned int>(unigramSize_xy) >= xlogx.size()) 
		    die("Something's wrong.  The word is way too big");

		float penalty = useX2Penalty? 
		    (x2[unigramSize_xy] - x2[unigramSize[x]] - x2[unigramSize[y]]):
		    (xlogx[unigramSize_xy] - xlogx[unigramSize[x]] - xlogx[unigramSize[y]]);

		float objective = 
		    f_xy * (1 + log(f_x * f_y / f_xy) - log_N - alpha) + 
		    0.5 * ((W + 1) * log(N - f_xy) - W * log_N) +
		    beta * f_xy * penalty;
		    // beta * f_xy * (unigramSize_xy * log(unigramSize_xy) - 
				   // unigramSize[x] * log(unigramSize[x]) - unigramSize[y] * log(unigramSize[y]));
		    
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

	    if (cli["verbose"]) {
		cerr << 
		    boost::format("%|6| %|-16| %|10|%|10.2f| %|10|") 
		    % iteration
		    % (lexicon.decode(x) + " " + lexicon.decode(y))
		    % N
		    % score
		    % bigram[bs.first] 
		    << '\n';
		    // << iteration << ' ' 
		    // << currentRatio << ' '
		    // << lexicon.decode(x) << lexicon.decode(y) << ' ' 
		    // << unigram[x].size() << ' '
		    // << unigram[y].size() << ' '
		    // << bigram[bs.first] << ' '
		    // << score << "\n";
	    }

	    // quit execution later if the threshould has been reached
	    if (score > threshold) quit = true;

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
	    string filename = boost::str(boost::format("%s.%05d") % prefix % iteration);
	    boost::filesystem::path filepath(outputPath);
	    filepath /= filename;

	    cerr << boost::format("Save output to %s") % filepath << '\n';
	    ofstream fout(filepath.c_str());
	    output_segmented_text(fout, text.begin() + 1, text.end(), lexicon, UNK, EOL);
	}
    }

    // Step 6: Output
    {
	string filename = boost::str(boost::format("%s.final") % prefix);
	boost::filesystem::path filepath(outputPath);
	filepath /= filename;

	cerr << boost::format("Save output to %s") % filepath << '\n';
	ofstream fout(filepath.c_str());
	output_segmented_text(fout, text.begin() + 1, text.end(), lexicon, UNK, EOL);
    }

    return 0;
}
