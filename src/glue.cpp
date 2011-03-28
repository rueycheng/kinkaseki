#include <iostream>
#include <boost/foreach.hpp>
#include "compat.h"
#include "kinkaseki/CLI.h"

#define foreach BOOST_FOREACH

#include "compat.h"
#include <iterator>
#include <boost/functional/hash.hpp>
#include <algorithm>
#include <ext/functional>

//--------------------------------------------------
// Home-made priority queue
//-------------------------------------------------- 
using std::vector;
using std::less;

template<typename _Tp, typename _Sequence = vector<_Tp>, 
    typename _Map = unordered_map<_Tp, typename _Sequence::distance_type>,
    typename _Compare = less<typename _Sequence::value_type> >
class PriorityQueue {
public:
    typedef typename _Sequence::value_type value_type;
    typedef typename _Sequence::reference reference;
    typedef typename _Sequence::const_reference const_reference;
    typedef typename _Sequence::size_type size_type;
    typedef _Sequence container_type;

protected:
    _Sequence c;
    _Map m;
    _Compare comp;

    typedef typename _Sequence::iterator _Iterator;
    typedef typename _Sequence::difference_type _Distance;
    typedef typename _Sequence::value_type _Value;

    void __push_heap(_Iterator __first, _Distance __holeIndex,
		     _Distance __topIndex, _Tp __value)
    {
	_Distance __parent = (__holeIndex - 1) / 2;
	while (__holeIndex > __topIndex
	       && __comp(*(__first + __parent), __value)) {
	    *(__first + __holeIndex) = *(__first + __parent);
	    __holeIndex = __parent;
	    __parent = (__holeIndex - 1) / 2;
	}

	*(__first + __holeIndex) = __value;
    }

    void push_heap(_Iterator __first, _Iterator __last)
    {
	_Value __value = *(__last - 1);
	__push_heap(__first, _Distance((__last - __first) - 1),
		    _Distance(0), __value);
    }

    void __pop_heap(_Iterator __first, _Iterator __last, _Iterator __result)
    {
	_Value __value = *__result;
	*__result = *__first;
	__adjust_heap(__first, _Distance(0), _Distance(__last - __first), __value);
    }

    void pop_heap(_Iterator __first, _Iterator __last)
    {
	--__last;
	__pop_heap(__first, __last, __last);
    }


    void __adjust_heap(_Iterator __first, _Distance __holeIndex,
		     _Distance __len, _Tp __value)
    {
	const _Distance __topIndex = __holeIndex;
	_Distance __secondChild = __holeIndex;
	while (__secondChild < (__len - 1) / 2) {
	    __secondChild = 2 * (__secondChild + 1);
	    if (comp(*(__first + __secondChild),
		     *(__first + (__secondChild - 1))))
		__secondChild--;
	    *(__first + __holeIndex) = *(__first + __secondChild);
	    __holeIndex = __secondChild;
	}

	if ((__len & 1) == 0 && __secondChild == (__len - 2) / 2) {
	    __secondChild = 2 * (__secondChild + 1);
	    *(__first + __holeIndex) = *(__first + __secondChild - 1);
	    __holeIndex = __secondChild - 1;
	}

	__push_heap(__first, __holeIndex, __topIndex, __value);
    }

    void make_heap(_Iterator __first, _Iterator __last)
    {
	if (__last - __first < 2) return;

	const _Distance __len = __last - __first;
	_Distance __parent = (__len - 2) / 2;
	while (true) {
	    _Value __value = *(__first + __parent);
	    __adjust_heap(__first, __parent, __len, __value);
	    if (__parent == 0) return;
	    __parent--;
	}
    }

public:
    explicit PriorityQueue(const _Compare& __comp = _Compare(),
			   const _Sequence& __s = _Sequence(),
			   const _Map& __m = _Map())
    : c(__s), m(__m), comp(__comp) 
    { 
	make_heap(c.begin(), c.end(), comp); 
    }

    template<typename _InputIterator> 
    PriorityQueue(_InputIterator __first, _InputIterator __last,
		  const _Compare& __comp = _Compare(), 
		  const _Sequence& __s = _Sequence(),
		  const _Map& __m = _Map()) 
    : c(__s), m(__m), comp(__comp) 
    {
	c.insert(c.end(), __first, __last);
	make_heap(c.begin(), c.end(), comp);
    }

    bool empty() const { return c.empty(); }
    size_type size() const { return c.size(); }
    const_reference top() const { return c.front(); }

    void push(const value_type& __x) {
	c.push_back(__x);
	push_heap(c.begin(), c.end(), comp);
    }

    void pop() {
	pop_heap(c.begin(), c.end(), comp);
	c.pop_back();
    }
};

//--------------------------------------------------
// Helper functions
//-------------------------------------------------- 
template<typename UnigramMap, typename BigramMap>
struct MoreEntropyGain {
    UnigramMap& umap;
    BigramMap& bmap;

    MoreEntropyGain(UnigramMap& um, BigramMap& bm): umap(um), bmap(bm) {}
    bool operator()(typename BigramMap::key_type x, typename BigramMap::key_type y) {
	unsigned int fx = bmap[x].size();
	unsigned int fy = bmap[y].size();

	return fx < fy;
    }
};

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

    Choose2nd() {}

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
    // Step 1: Read input and store them as a token stream
    //         Note the ``special tokens''
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
    // We assume the size of the text stream fits into a 4-byte integer.
    // From now on, we'll call each token as a 'Unigram'.
    typedef unsigned int Unigram;
    typedef pair<unsigned int, unsigned int> Bigram;
    typedef vector<unsigned int> PostingList;

    // NOTE: We keep track of the positions for each unigram (and thus we know
    // the frequencies).  For bigram, we only save the counts.
    typedef vector<PostingList> UnigramIndex;
    typedef unordered_map<Bigram, unsigned int, boost::hash<Bigram> > BigramIndex;

    UnigramIndex unigram(lexicon.size());
    BigramIndex bigram(lexicon.size());

    {
	vector<unsigned int>::iterator first = text.begin(), last = text.end();
	vector<unsigned int>::iterator iter = first;

	Unigram prev = 1, curr;

	while (iter != last) {
	    curr = *iter;

	    if (curr != 1) {
		unigram[curr].push_back(distance(first, iter));
		if (prev != 1) bigram[Bigram(prev, curr)]++;
	    }

	    prev = curr;
	    ++iter;
	}
    }

    // Step 3: Populate the top-k bigrams
    //
    // More detail later
    typedef pair<Bigram, float> BigramScore;

    vector<BigramScore> score(bigram.size());

    {
	BigramIndex::iterator iter = bigram.begin(), last = bigram.end();
	while (iter != last) { 
	    int f_x = unigram[iter->first.first].size();
	    int f_y = unigram[iter->first.second].size();
	    int f_xy = iter->second;
	    float g = std::log(f_x - f_xy) + std::log(f_y - f_xy) - std::log(f_xy);

	    score.push_back(BigramScore(iter->first, g));
	}

	partial_sort(
	    score.begin(), 
	    score.begin() + 10, 
	    score.end(), 
	    choose2nd<BigramScore>(std::less<float>())
	);
    }

//--------------------------------------------------
//     vector<Bigram> heap;
//     heap.reserve(bigram.size());
// 
//     MoreEntropyGain<
// 	unordered_map<Unigram, PostingList, boost::hash<Unigram> >,
// 	unordered_map<Bigram, PostingList, boost::hash<Bigram> > 
// 	>
// 	compare(unigram, bigram);
// 
//     {
// 	unordered_map<Bigram, PostingList, boost::hash<Bigram> >
// 	    ::const_iterator iter = bigram.begin(), last = bigram.end();
// 
// 	while (iter != last)
// 	    heap.push_back(iter++->first);
// 
// 	make_heap(heap.begin(), heap.end(), compare);
// 
// 	while (!heap.empty()) {
// 	    Bigram& top = heap.front();
// 	    cout << bigram[top].size() << ' ' << heap.size() << ' ' << top.first << ' ' << top.second << ' ' << "\n";
// 
// 	    pop_heap(heap.begin(), heap.end(), compare);
// 	    heap.pop_back();
// 	}
//     }
//-------------------------------------------------- 

    return 0;
}
