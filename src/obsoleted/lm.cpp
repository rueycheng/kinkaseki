#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <cmath>
#include <numeric>

#include "expatpp.h"
#include "expatpp_callback.h"
#include "magicbox.h"
#include "lm.h"

#include <boost/algorithm/string.hpp>
#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>

//--------------------------------------------------
// autoconf stuff
//-------------------------------------------------- 
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#ifdef HAVE_UNORDERED_MAP
#include <unordered_map>
#else
#ifdef HAVE_TR1_UNORDERED_MAP
#include <tr1/unordered_map>
using namespace std::tr1;
#endif
#endif

using namespace std;
using namespace magicbox;

//--------------------------------------------------
// TRECDocument
//-------------------------------------------------- 
struct TRECDocument {
    string content;
    string docno;
    vector<string> facets;

    void clear() {
	content.clear();
	docno.clear();
	facets.clear();
    }

    void normalize() {
	boost::trim(content);
	boost::trim(docno);
	foreach (string& facet, facets) {
	    boost::trim(facet);
	}
    }
};

//--------------------------------------------------
// TRECReader
//-------------------------------------------------- 
class TRECReader: public ContextualSAXHandler {
    TRECDocument doc;
    vector<TRECDocument> collection;

    istream& in;
    ostream& err;

    bool sentStart;
    bool sentEnd;

public:
    TRECReader(istream& in, ostream& err = std::cerr): in(in), err(err), sentStart(false), sentEnd(false) {}

    void enterContext(string context, string name, attribute_map attr) {
	boost::to_lower(context);
	if (context == "/root/doc/facet") doc.facets.push_back("");
    }

    void leaveContext(string context, string name) {
	boost::to_lower(context);
	if (context == "/root/doc") {
	    doc.normalize();
	    collection.push_back(doc);
	    doc.clear();
	}
    }

    void text(string context, string text) {
	boost::to_lower(context);
	if (context == "/root/doc/docno") doc.docno += text;
	else if (context == "/root/doc/facet") doc.facets.back() += text;
	else if (context == "/root/doc/lang" || context == "/root/doc/date") ;
	else if (context.find("/root/doc") == 0) doc.content += text;
    }

    const vector<TRECDocument>& getCollection() { return collection; }
    vector<TRECDocument>::size_type size() { return collection.size(); }
    void clear() { collection.clear(); }

    bool foundNewDocument() {
	char buf[BUFSIZE];

	string start = "<?xml version=\"1.0\" encoding=\"UTF-8\"?><root>";
	string end = "</root>";

	if (sentEnd) return false;

	if (!sentStart) {
	    parser.parse(start.c_str(), start.size(), false);
	    sentStart = true;
	}

	vector<TRECDocument>::size_type alreadyFound = collection.size();
	while (alreadyFound == collection.size() && in.read(buf, BUFSIZE))
	    parser.parse(buf, in.gcount(), false) || carp(err, buf);

	if (alreadyFound - collection.size() > 0) return true;

	parser.parse(buf, in.gcount(), false);

	if (!sentEnd) {
	    parser.parse(end.c_str(), end.size(), true);
	    sentEnd = true;
	}

	return alreadyFound - collection.size() > 0;
    }
};

//--------------------------------------------------
// Vocabulary
//-------------------------------------------------- 
template<typename __word = unsigned int>
class Vocabulary {
    unordered_map<string, __word> vocab;
    vector<string> rvocab;

public:
    Vocabulary() { rvocab.push_back("<unk>"); }

    Vocabulary(istream& in) {
	Vocabulary();

	string word;
	while (getline(in, word)) {
	    rvocab.push_back(word);
	    vocab[word] = rvocab.size() - 1;
	}
    }

    void save(ostream& o) { 
	foreach (const string& word, rvocab) { o << word << endl; } 
    }

    __word encode(const string& word) {
	if (vocab.find(word) != vocab.end()) return vocab[word];
	rvocab.push_back(word);
	return vocab[word] = rvocab.size() - 1;
    }

    string decode(unsigned int code) { return rvocab.at(code); }

    __word size() { return rvocab.size(); }
};

//--------------------------------------------------
// OrderedPostingList
//-------------------------------------------------- 
template<typename __doc = unsigned int, 
    typename __count = unsigned short, 
    class __tokenizer = CJKVTokenizer>
class OrderedPostingList {
public:
    typedef __doc doc_type;
    typedef __count count_type;
    typedef __tokenizer tokenizer_type;

private:
    typedef pair<doc_type, count_type> item;
    typedef vector<item> posting;
    typedef unordered_map<string, posting> posting_list;
    typedef unordered_map<string, count_type> freq_map;

    posting_list pl;

public:
    unsigned int addDocument(doc_type docid, const string& text) {
	tokenizer_type tok(text);
	freq_map tf;

	unsigned int size = 0;
	foreach (const CJKVTokenizer::value_type& token, tok) { 
	    ++tf[boost::to_lower_copy(token)]; 
	    ++size;
	}

	foreach (const typename freq_map::value_type& p, tf) {
	    pl[p.first].push_back(item(docid, p.second));
	}

	return size;
    }

    void save(ostream& out) {
	vector<typename posting_list::key_type> key;
	foreach (const typename posting_list::value_type& p, pl) { key.push_back(p.first); }
	stable_sort(key.begin(), key.end());

	foreach (const typename posting_list::key_type& k, key) {
	    out << k << ' ' << pl[k].size() << ' ';
	    foreach (const item& pp, pl[k]) { out << pp.first << ' ' << pp.second << ' '; }
	    out << "\n";
	}
    }
};

//--------------------------------------------------
// ifstream_pool
//-------------------------------------------------- 
class ifstream_pool {
    unordered_map<string, ifstream*> m;

public:
    ifstream* operator[](string filename) {
	if (m.find(filename) == m.end()) 
	    m[filename] = new ifstream(filename.c_str());
	return m[filename];
    }

    void close(string filename) {
	if (m.find(filename) == m.end()) return;
	delete m[filename];
	m.erase(filename);
    }

    ~ifstream_pool() {
	typedef pair<string, ifstream*> string_ifstreamp;
	foreach (const string_ifstreamp& p, m) { delete p.second; }
	m.clear();
    }
};

//--------------------------------------------------
// ifstream_node
//-------------------------------------------------- 
template<typename T> struct ifstream_node {
    ifstream* in;
    T data;

    ifstream_node(ifstream* in, const T& data): in(in), data(data) {}
    bool operator<(const ifstream_node<T>& b) const { return data < b.data; }
};

template<typename T> struct ifstream_node_reverse_cmp {
    bool operator()(const ifstream_node<T>& lhs, const ifstream_node<T>& rhs) const { return rhs < lhs; }
};

//--------------------------------------------------
// record (client code)
//-------------------------------------------------- 
struct record {
    string key;
    vector<unsigned int> values;
    vector<short> counts;

    bool operator<(const record& b) const {
	return key < b.key || 
	    (key == b.key && values.size() > 0 && b.values.size() > 0 && values.front() < b.values.front());
    }

    record& operator+=(const record& b) {
	if (key != b.key) throw 0;
	copy(b.values.begin(), b.values.end(), back_inserter(values));
	copy(b.counts.begin(), b.counts.end(), back_inserter(counts));
	return *this;
    }

    void load(istream& in) {
	key.clear();
	values.clear();
	counts.clear();

	if (in) {
	    unsigned int size = 0;
	    in >> key >> size;

	    for (unsigned int i = 0; i < size; ++i) {
		unsigned int d;
		short c;

		in >> d >> c;
		values.push_back(d);
		counts.push_back(c);
	    }
	}
    }

    void save(ostream& out) const {
	out << key << ' ' << values.size() << ' ';
	vector<unsigned int>::const_iterator viter = values.begin();
	vector<short>::const_iterator citer = counts.begin();

	while (viter != values.end()) out << *viter++ << ' ' << *citer++ << ' ';
    }
};

namespace std {
    istream& operator>>(istream& in, record& r) {
	r.load(in);
	return in;
    }

    ostream& operator<<(ostream& out, const record& r) {
	r.save(out);
	return out;
    }
}

//--------------------------------------------------
// Utility
//-------------------------------------------------- 
template<typename T, typename S> 
struct cmp_first {
    bool operator()(const pair<T,S>& lhs, const pair<T,S>& rhs) {
	return lhs.first < rhs.first;
    }
};

template<typename T, typename S> 
struct cmp_second {
    bool operator()(const pair<T,S>& lhs, const pair<T,S>& rhs) {
	return lhs.second < rhs.second;
    }
};

template<typename T, typename S> 
struct rcmp_first {
    bool operator()(const pair<T,S>& lhs, const pair<T,S>& rhs) {
	return rhs.first < lhs.first;
    }
};

template<typename T, typename S> 
struct rcmp_second {
    bool operator()(const pair<T,S>& lhs, const pair<T,S>& rhs) {
	return rhs.second < lhs.second;
    }
};

//--------------------------------------------------
// Main program
//-------------------------------------------------- 
int main(int argc, char** argv) {

    // Getopt
    string model;
    vector<string> input;
    float mu = 2500;

    Getopt g(argc, argv);
    g   << $("train", "Train the model and save to the index")
	<< $("test", "Test the saved model with the input")
	<< $(&mu, "mu", "Specify the smoothing parameter `mu'.")
	<< $(&model, "model,m", "Path to the model")
	<< $(&input, "input", "", -1)
	<< $$$("{ --train [options..] files.. | --test [options..] }");

    if (g["train"] && g["test"]) die("Specify either --train or --test.  Not both");
    if (!g["train"] && !g["test"]) { cout << g; exit(0); }

    // 0) Stuff
    if (model.empty()) die("Missing model directory");
    if (input.empty()) input.push_back("-");

    namespace fs = boost::filesystem;

    // 1) Dotch?
    if (g["train"]) {
	//--------------------------------------------------
	// Training
	//-------------------------------------------------- 
	// 2.0) Create model directory
	fs::create_directory(model);
	if (!fs::exists(model)) die("Cannot create directory " + model);
	fs::path basedir(model);

	// 2.1) Produce a number of inverted runs
	{
	    vector<string> docno;
	    vector<unsigned int> docsize;
	    unsigned int runno = 0;

	    foreach (const string& filename, input) {
		AutoIn in(filename);
		TRECReader reader(in());

		while (reader.foundNewDocument()) { 
		    if (reader.size() < 100000) continue;

		    fs::ofstream pl_out(basedir / str(boost::format("runfile.%04d") % ++runno));
		    cerr << boost::format("Save run #%d") % runno << "\n";
		    OrderedPostingList<> pl;
		    foreach (const TRECDocument& doc, reader.getCollection()) { 
			docno.push_back(doc.docno);
			docsize.push_back(pl.addDocument(docno.size(), doc.content)); 
		    }

		    pl.save(pl_out);
		    reader.clear();
		}

		if (reader.size()) { 
		    fs::ofstream pl_out(basedir / str(boost::format("runfile.%04d") % ++runno));
		    cerr << boost::format("Save run #%d") % runno << "\n";
		    OrderedPostingList<> pl;
		    foreach (const TRECDocument& doc, reader.getCollection()) { 
			docno.push_back(doc.docno);
			docsize.push_back(pl.addDocument(docno.size(), doc.content)); 
		    }

		    pl.save(pl_out);
		    reader.clear();
		}
	    }

	    // Output docid map
	    cerr << "Save docno" << "\n";
	    fs::ofstream docno_out(basedir / "docno");

	    docno_out << "<collection>" << "\n";
	    foreach (const string& filename, docno) { docno_out << filename << "\n"; }
	    docno.clear();

	    // Output docsize
	    cerr << "Save docsize" << "\n";
	    fs::ofstream docsize_out(basedir / "docsize");

	    unsigned int colsize = accumulate(docsize.begin(), docsize.end(), 0);
	    docsize_out << pack(colsize);
	    foreach (const unsigned int& size, docsize) { docsize_out << pack(size); }
	    docsize.clear();
	}

	// 2.2) Merge all the runs
	vector<fs::path> runfiles;

	// NOTE: boost::foreach does not seem to work with directory_iterator (or any other ways to do it?)
	//       Consider std::remove_copy_if()
	for (fs::directory_iterator item(basedir), end; item != end; ++item)
	    if (boost::starts_with(item->leaf(), "runfile.")) runfiles.push_back(item->path());
	stable_sort(runfiles.begin(), runfiles.end());

	if (runfiles.size() == 1) {
	    cerr << "Copy and save run #1" << "\n";
	    fs::ofstream run_out(basedir / "all.runfile");
	    fs::ifstream run_in(runfiles.front());
	    run_out << run_in.rdbuf();
	}
	else {
	    cerr << "Merge all runs" << "\n";
	    typedef ifstream_node<record> record_node;

	    fs::ofstream run_out(basedir / "all.runfile");
	    ifstream_pool pool;
	    vector<record_node> pq;

	    foreach (const fs::path& p, runfiles) { 
		ifstream* in = pool[p.string()];
		record r;
		if (in && (*in >> r)) pq.push_back(record_node(in, r));
	    }

	    // Output
	    string last_key;  // Starts empty
	    vector<record> buf;

	    ifstream_node_reverse_cmp<record> smaller_first;
	    make_heap(pq.begin(), pq.end(), smaller_first);
	    while (pq.size()) {
		pop_heap(pq.begin(), pq.end(), smaller_first);

		if (pq.back().data.key != last_key) {
		    if (buf.size()) {
			vector<record>::const_iterator iter = buf.begin();
			record bb = *iter++;
			while (iter != buf.end()) bb += *iter++;
			run_out << bb << "\n";
			buf.clear();
		    }
		    last_key = pq.back().data.key;
		}

		buf.push_back(pq.back().data);

		if (*(pq.back().in) >> pq.back().data) push_heap(pq.begin(), pq.end(), smaller_first);
		else pq.pop_back();
	    }

	    if (buf.size()) {
		vector<record>::const_iterator iter = buf.begin();
		record bb = *iter++;
		while (iter != buf.end()) bb += *iter++;
		run_out << bb << "\n";
		buf.clear();
	    }
	}

	// 2.3) Produce the inverted index
	{
	    cerr << "Save inverted_index" << "\n";
	    fs::ofstream inv_out(basedir / "inverted_index");
	    fs::ofstream vocab_out(basedir / "vocabulary");
	    fs::ifstream run_in(basedir / "all.runfile");

	    unsigned int offset = 0;

	    record r;
	    while (run_in >> r) {
		vocab_out << r.key << ' ' << offset << "\n";

		unsigned int df = r.values.size();
		inv_out << pack(df);

		unsigned int tf = accumulate(r.counts.begin(), r.counts.end(), 0);
		inv_out << pack(tf);

		vector<unsigned int>::const_iterator viter = r.values.begin();
		vector<short>::const_iterator citer = r.counts.begin();
		while (viter != r.values.end()) inv_out << pack(*viter++) << pack(*citer++);

		offset += (2 + df)  * sizeof(unsigned int) + df * sizeof(short);
	    }
	}
    }
    else if (g["test"]) {
	//--------------------------------------------------
	// Test
	//-------------------------------------------------- 
	// 3.0) Load model directory
	if (!fs::exists(model)) die("The model directory " + model + " does not exist");
	fs::path basedir(model);

	// 3.1) Load vocabulary and open the inverted index
	unordered_map<string, unsigned int> term_offset;

	{
	    fs::ifstream vocab_in(basedir / "vocabulary");
	    string term;
	    unsigned int offset;
	    while (vocab_in >> term >> offset) term_offset[term] = offset;
	}

	fs::ifstream inv_in(basedir / "inverted_index");

	vector<unsigned int> docsize;

	{
	    fs::ifstream docsize_in(basedir / "docsize");
	    unsigned int size;
	    while (docsize_in >> unpack(size)) docsize.push_back(size);
	}

	vector<string> docno;

	{
	    fs::ifstream docno_in(basedir / "docno");
	    copy(istream_iterator<string>(docno_in), istream_iterator<string>(), back_inserter(docno));
	}

	// 3.2) Process the query
	ostringstream o;
	o << cin.rdbuf();
	string line = o.str();
	CJKVTokenizer tok(line);

	vector<unsigned int> doc1, doc2;
	vector<float> score1, score2;

	bool turn = true;

	foreach (const CJKVTokenizer::value_type& token, tok) { 
	    string term = token;
	    boost::to_lower(term);

	    bool oov = term_offset.find(term) == term_offset.end();

	    if (oov) {
		//--------------------------------------------------
		// cout << term << " (OOV)\n";
		//-------------------------------------------------- 
		;
	    }
	    else {
		vector<unsigned int>& result_doc = turn? doc1: doc2;
		vector<unsigned int>& result_doc_prev = turn? doc2: doc1;
		vector<float>& result_score = turn? score1: score2;
		vector<float>& result_score_prev = turn? score2: score1;
		turn ^= true;

		// Go.
		//--------------------------------------------------
		// cout << term << ": " << term_offset[term] << "\n";
		//-------------------------------------------------- 
		inv_in.seekg(term_offset[term], ios_base::beg);

		vector<unsigned int> docs;
		vector<float> scores;
		unsigned int df, tf;
		inv_in >> unpack(df) >> unpack(tf);

		//--------------------------------------------------
		// cout << "  " << df << ' ' <<  tf << "\n";
		//-------------------------------------------------- 
		float pref = float(tf) / docsize[0];
		float log_pref = log(pref);

		for (unsigned int i = 0; i < df; i++) {
		    unsigned int doc;
		    short count;
		    inv_in >> unpack(doc) >> unpack(count);

		    float p = log((count + mu * pref) / (docsize[doc] + mu));
		    float pdiff = p - log_pref;

		    docs.push_back(doc);
		    scores.push_back(pdiff);
		}

		// Merge
		result_doc.clear();
		result_score.clear();

		vector<unsigned int>::const_iterator d1 = result_doc_prev.begin();
		vector<unsigned int>::const_iterator d1_end = result_doc_prev.end();
		vector<unsigned int>::const_iterator d2 = docs.begin();
		vector<unsigned int>::const_iterator d2_end = docs.end();
		vector<float>::const_iterator s1 = result_score_prev.begin();
		vector<float>::const_iterator s2 = scores.begin();

		while (d1 != d1_end && d2 != d2_end) {
		    if (*d1 < *d2) {
			result_doc.push_back(*d1++);
			result_score.push_back(*s1++);
		    }
		    else if (*d2 < *d1) {
			result_doc.push_back(*d2++);
			result_score.push_back(*s2++);
		    }
		    else {
			result_doc.push_back(*d1++);
			result_score.push_back(*s1++ + *s2++);
			++d2;
		    }
		}

		if (d1 == d1_end) {
		    while (d2 != d2_end) {
			result_doc.push_back(*d2++);
			result_score.push_back(*s2++);
		    }
		}
		else {
		    while (d1 != d1_end) {
			result_doc.push_back(*d1++);
			result_score.push_back(*s1++);
		    }
		}
	    }
	}


	vector<unsigned int>& doc_prev = turn? doc2: doc1;
	vector<float>& score_prev = turn? score2: score1;

	vector< pair<unsigned int, float> > ranked;

	vector<unsigned int>::const_iterator dp = doc_prev.begin();
	vector<float>::const_iterator sp = score_prev.begin();
	while (dp != doc_prev.end()) 
	    ranked.push_back(make_pair(*dp++, *sp++));

	stable_sort(ranked.begin(), ranked.end(), rcmp_second<unsigned int, float>());

	typedef pair<unsigned int, float> item ;
	foreach (const item& ri, ranked) {
	    cout << docno[ri.first] << ' ' << ri.second << "\n";
	}
    }

    return 0;
}
