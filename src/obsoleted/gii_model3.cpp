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
// Vocabulary
//-------------------------------------------------- 
template<typename __word = unsigned int>
class vocabulary {
    unordered_map<string, __word> vocab;
    vector<string> rvocab;
    vector<unsigned int> counts;
    vector<unsigned int> doc_counts;
    unsigned int total;

    __word grow(const string& w, unsigned int c = 0, unsigned int dc = 0) {
	rvocab.push_back(w);
	counts.push_back(c);
	doc_counts.push_back(dc);
	return vocab[w] = rvocab.size() - 1;
    }

public:
    vocabulary(): total(0) { grow("<unk>"); }
    vocabulary(istream& in, unsigned int size, unsigned int total) { load(in, size, total); }

    void update_counts(__word id, unsigned short num) { 
	++doc_counts.at(id); 
	counts.at(id) += num;
	total += num;
    }

    unsigned int df(__word id) { return doc_counts.at(id); }
    unsigned int tf(__word id) { return counts.at(id); }
    unsigned int total_word_count() const { return total; }
    unsigned int total_doc_count() const { return doc_counts.size(); }
    __word size() { return rvocab.size(); }

    void save(ostream& o) { 
	vector<string>::const_iterator w_iter = rvocab.begin();
	vector<unsigned int>::const_iterator c_iter = counts.begin();
	vector<unsigned int>::const_iterator dc_iter = doc_counts.begin();
	while (w_iter != rvocab.end())
	    o << *w_iter++ << ' ' << *c_iter++ << ' ' << *dc_iter++ << "\n";
    }

    void load(istream& in, unsigned int size, unsigned int total) {
	vocab.clear();
	rvocab.clear();
	counts.clear();
	doc_counts.clear();

	string line;
	while (size-- > 0 && getline(in, line)) {
	    istringstream iss(line);
	    string word;
	    unsigned int count = 0;
	    unsigned int doc_count = 0;
	    iss >> word >> count >> doc_count;
	    grow(word, count, doc_count);
	}

	this->total = total;
    }

    __word encode(const string& word) {
	if (vocab.find(word) != vocab.end()) return vocab[word];
	return grow(word);
    }

    __word safe_encode(const string& word) {
	if (vocab.find(word) != vocab.end()) return vocab[word];
	return 0;
    }

    string decode(unsigned int code) { return rvocab.at(code); }
};

//--------------------------------------------------
// posting
//-------------------------------------------------- 
struct posting {
    typedef pair<unsigned int, unsigned short> count;
    typedef pair<unsigned int, float> weight;

    string name;
    unsigned short total;
    vector<count> counts;
};

struct document_posting: public posting {};
struct term_posting: public posting {};

namespace std {
    ostream& operator<<(ostream& o, const posting& p) {
	o << p.name << ' ' << p.counts.size() << ' ' << p.total << ' ';
	foreach (const posting::count& c, p.counts) { o << c.first << ' ' << c.second << ' '; }
	return o;
    }

    istream& operator>>(istream& in, posting& p) {
	unsigned int size;
	in >> p.name >> size >> p.total;
	
	unsigned int word_id, count;
	for (unsigned int i = 0; i < size && in >> word_id >> count; ++i)
	    p.counts.push_back(posting::count(word_id, count));
	return in;
    }
}

//--------------------------------------------------
// Main program
//-------------------------------------------------- 
int main(int argc, char** argv) {

    // Getopt
    vector<string> input;
    float mu = 2500;
    string method = "tfidf";
    string model = "model.unnamed";
    unsigned int top_n = 100;

    Getopt g(argc, argv);
    g   << $("index", "Create searchable model")
        << $("query", "Query the saved model")
	<< $(&model, "model,m", "Path to the model directory")
	<< $(&method, "method", "Ranking method: tfidf or lm")
	<< $(&mu, "mu", "Specify the parameter 'mu' (for LM)")
	<< $(&top_n, "top-n,n", "Show only top N results")
	<< $(&input, "input", "@@ IT DOES NOT MATTER WHAT I PUT IN HERE @@", -1)
	<< $$$("files..");

    if (input.empty()) input.push_back("-");
    namespace fs = boost::filesystem;

    if (!g["index"] && !g["query"]) bye("Specify either --index or --query");
    if (g["index"] && g["query"]) bye("Cannot run in both mode.  Specify either one");

    //--------------------------------------------------
    // Indexing
    //-------------------------------------------------- 

    //--------------------------------------------------
    // NOTE: Employing in-memory inversion algorithm since this is quite in a hurry
    //-------------------------------------------------- 
    if (g["index"]) {
	// 0) Make directory
	fs::create_directory(model);
	if (!fs::exists(model)) die("Cannot create directory " + model);
	fs::path basedir(model);

	// Collection-wide
	vocabulary<> vocab;
	vector<document_posting> docsize;
	vector<document_posting> facets;
	vector<term_posting> terms;

	foreach (const string& filename, input) {
	    AutoIn in(filename);
	    string line;

	    while (getline(in(), line)) {
		istringstream iss(line);
		string docno;

		iss >> docno;

		// Process docno
		cerr << docno << "\n";
		if (docno.empty()) continue;

		if (!boost::ends_with(docno, ":facet")) {
		    // Case 1: Regular texts
		    unordered_map<unsigned int, unsigned short> counts;
		    unsigned short total = 0; // Document size

		    // Accumulate word counts
		    {
			string word;
			while (iss >> word) {
			    ++counts[vocab.encode(word)];
			    ++total; 
			}
		    }

		    // Update tf, df, and total_word_count all at once
		    foreach (const posting::count& c, counts) { vocab.update_counts(c.first, c.second); }

		    // Save docsize
		    docsize.push_back(document_posting());
		    document_posting& dest = docsize.back();
		    dest.name = docno;
		    dest.total = total;

		//--------------------------------------------------
		//     copy(counts.begin(), counts.end(), back_inserter(dest.counts));
		//     stable_sort(dest.counts.begin(), dest.counts.end(), first_cmp());
		//-------------------------------------------------- 
		}
		else {
		    // Case 2: Facets
		    // 
		}
	    }
	}

	// Save vocabulary
	{
	    cerr << "Save vocabulary" << "\n";
	    fs::ofstream vocab_out(basedir / "vocabulary");

	    vocab_out << "@vocab " << vocab.size() << ' ' << vocab.total_word_count() << "\n";
	    vocab.save(vocab_out);
	}

	// Save docsize
	{
	    cerr << "Save docsize" << "\n";
	    fs::ofstream docsize_out(basedir / "docsize");

	    docsize_out << "@docsize " << docsize.size() << "\n";
	    foreach (const posting& dp, docsize) { docsize_out << dp << "\n"; }
	}

	return 0;
    }

//--------------------------------------------------
//     //--------------------------------------------------
//     // Querying the model
//     //-------------------------------------------------- 
//     if (g["test"]) {
// 	// 0) Load directory
// 	if (!fs::exists(model)) die("Cannot open directory " + model);
// 	fs::path basedir(model);
// 
// 	// Load vocabulary
// 	vocabulary<> vocab;
// 	unsigned int total_word_count;
// 
// 	{
// 	    cerr << "Load vocabulary" << "\n";
// 	    fs::ifstream vocab_in(basedir / "vocabulary");
// 
// 	    string line, head;
// 	    if (!getline(vocab_in, line)) throw InvalidFormat();
// 
// 	    unsigned int vocab_size;
// 	    istringstream iss(line);
// 	    if (!(iss >> head >> vocab_size >> total_word_count)) throw InvalidFormat();
// 	    vocab.load(vocab_in, vocab_size, total_word_count);
// 	}
// 
// 	// Translate query
// 	vector<document_posting> queries;
// 
// 	{
// 	    cerr << "Make query" << "\n";
// 
// 	    string line;
// 	    while (getline(cin, line)) {
// 		istringstream iss(line);
// 		string docno;
// 
// 		iss >> docno;
// 		if (boost::ends_with(docno, ":facet")) continue;
// 
// 		unordered_map<unsigned int, unsigned short> counts;
// 		unsigned short total = 0; // Document size
// 
// 		// Accumulate word counts
// 		{
// 		    string word;
// 		    while (iss >> word) {
// 			++counts[vocab.safe_encode(word)]; // NOTE: Use safe_encode to prevent adding new words
// 			++total; 
// 		    }
// 		}
// 
// 		// Save it back into the posting
// 		queries.push_back(document_posting());
// 		document_posting& current = queries.back();
// 
// 		current.name = docno;
// 		current.total = total;
// 		copy(counts.begin(), counts.end(), back_inserter(current.counts));
// 		stable_sort(current.counts.begin(), current.counts.end(), first_cmp());
// 	    }
// 	}
// 
// 	// Load collection
// 	vector<document_posting> documents;
// 	unsigned int total_document_count;
// 
// 	{
// 	    cerr << "Load collection" << "\n";
// 	    fs::ifstream collection_in(basedir / "collection");
// 
// 	    string line, head;
// 	    if (!getline(collection_in, line)) throw InvalidFormat();
// 	    
// 	    istringstream iss(line);
// 	    if (!(iss >> head >> total_document_count)) throw InvalidFormat();
// 
// 	    for (unsigned int i = 0; i < total_document_count; ++i) {
// 		getline(collection_in, line);
// 		istringstream iss2(line);
// 
// 		documents.push_back(document_posting());
// 		iss2 >> documents.back();
// 	    }
// 	}
// 
// 	// Now we are ready
// 	if (method == "tfidf") {
// 	    float logN = log(total_document_count);
// 
// 	    // Process queries
// 	    vector<string> queryno;
// 	    vector<vector<posting::weight> > queryv;
// 	    vector<float> queryn;
// 
// 	    foreach (const document_posting& q, queries) {
// 		queryno.push_back(q.name);
// 		queryv.push_back(vector<posting::weight>());
// 		vector<posting::weight>& qv = queryv.back();
// 		queryn.push_back(0);
// 		float& nq = queryn.back();
// 
// 		foreach (const posting::count& c, q.counts) {
// 		    float part = log(1 + c.second); // lnc
// 		    qv.push_back(posting::weight(c.first, part));
// 		    nq += part; // Get sum of squares
// 		}
// 		nq = sqrt(nq); // .. and take squares
// 	    }
// 
// 	    // Process documents
// 	    vector<string> docno;
// 	    vector<vector<posting::weight> > docv;
// 	    vector<float> docn;
// 
// 	    foreach (const document_posting& d, documents) {
// 		docno.push_back(d.name);
// 		docv.push_back(vector<posting::weight>());
// 		vector<posting::weight>& dv = docv.back();
// 		docn.push_back(0);
// 		float& nd = docn.back();
// 
// 		foreach (const posting::count& c, d.counts) {
// 		    float part = log(1 + c.second) * (logN - log(vocab.df(c.first))); // ltc
// 		    dv.push_back(posting::weight(c.first, part));
// 		    nd += part; // Get sum of squares
// 		}
// 		nd = sqrt(nd); // .. and take squares
// 	    }
// 
// 	    // Fantastic.  Lose memory
// 	    queries.clear(); 
// 	    documents.clear();
// 
// 	    for (unsigned int i = 0; i < queryno.size(); ++i) {
// 		cerr << queryno.at(i) << '\n';
// 
// 		typedef pair<string, float> item;
// 		vector<item> ranklist;
// 
// 		for (unsigned int j = 0; j < docno.size(); ++j) {
// 		    // Now do the `dot'
// 		    vector<posting::weight>::const_iterator qq = queryv.at(i).begin();
// 		    vector<posting::weight>::const_iterator qq_end = queryv.at(i).end();
// 		    vector<posting::weight>::const_iterator dd = docv.at(j).begin();
// 		    vector<posting::weight>::const_iterator dd_end = docv.at(j).end();
// 
// 		    float dotsum = 0;
// 		    while (qq != qq_end && dd != dd_end) {
// 			if ((*qq).first < (*dd).first) ++qq;
// 			else if ((*qq).first > (*dd).first) ++dd;
// 			else dotsum += (*qq++).second * (*dd++).second;
// 		    }
// 
// 		    ranklist.push_back(item(docno.at(j), dotsum / queryn.at(i) / docn.at(j)));
// 		}
// 
// 		if (ranklist.size() > top_n) {
// 		    nth_element(ranklist.begin(), ranklist.begin() + top_n, ranklist.end(), second_rcmp());
// 		    ranklist.erase(ranklist.begin() + top_n, ranklist.end());
// 		}
// 
// 		stable_sort(ranklist.begin(), ranklist.end(), second_rcmp());
// 		foreach (const item& r, ranklist) {
// 		    cout << queryno.at(i) << ' ' << r.first << ' ' << r.second << '\n';
// 		}
// 	    }
// 	}
// 	else if (method == "lm") {
// 	    // Process queries
// 	    vector<string> queryno;
// 	    vector<vector<posting::weight> > queryv;
// 	    vector<float> queryn;
// 
// 	    foreach (const document_posting& q, queries) {
// 		queryno.push_back(q.name);
// 		queryv.push_back(vector<posting::weight>());
// 		vector<posting::weight>& qv = queryv.back();
// 		queryn.push_back(q.total); // |q|
// 
// 		foreach (const posting::count& c, q.counts) {
// 		    qv.push_back(posting::weight(c.first, c.second)); // f_{iq}
// 		}
// 	    }
// 
// 	    // Process documents
// 	    vector<string> docno;
// 	    vector<vector<posting::weight> > docv;
// 	    vector<float> docn;
// 
// 	    foreach (const document_posting& d, documents) {
// 		docno.push_back(d.name);
// 		docv.push_back(vector<posting::weight>());
// 		vector<posting::weight>& dv = docv.back();
// 		docn.push_back(log(d.total + mu)); // \log(|d_j| + \mu)
// 
// 		foreach (const posting::count& c, d.counts) {
// 		    float part = log(c.second + mu * vocab.tf(c.first) / vocab.total_word_count() ); // \log(f_{ij} + \mu * f_i / f)
// 		    dv.push_back(posting::weight(c.first, part));
// 		}
// 	    }
// 
// 	    // Fantastic.  Lose memory
// 	    queries.clear(); 
// 	    documents.clear();
// 
// 	    // Go.
// 	    for (unsigned int i = 0; i < queryno.size(); ++i) {
// 		cerr << queryno.at(i) << '\n';
// 
// 		typedef pair<string, float> item;
// 		vector<item> ranklist;
// 
// 		for (unsigned int j = 0; j < docno.size(); ++j) {
// 		    // Now do the `dot'
// 		    vector<posting::weight>::const_iterator qq = queryv.at(i).begin();
// 		    vector<posting::weight>::const_iterator qq_end = queryv.at(i).end();
// 		    vector<posting::weight>::const_iterator dd = docv.at(j).begin();
// 		    vector<posting::weight>::const_iterator dd_end = docv.at(j).end();
// 
// 		    float dotsum = 0;
// 		    while (qq != qq_end && dd != dd_end) {
// 			if ((*qq).first < (*dd).first) {
// 			    unsigned int word_id = (*qq).first;
// 			    dotsum += (*qq++).second * log(mu * (1 + vocab.tf(word_id)) / vocab.total_word_count());
// 			}
// 			else if ((*qq).first > (*dd).first) ++dd;
// 			else dotsum += (*qq++).second * (*dd++).second;
// 		    }
// 
// 		    ranklist.push_back(item(docno.at(j), dotsum - queryn.at(i) * docn.at(j)));
// 		}
// 
// 		if (ranklist.size() > top_n) {
// 		    nth_element(ranklist.begin(), ranklist.begin() + top_n, ranklist.end(), second_rcmp());
// 		    ranklist.erase(ranklist.begin() + top_n, ranklist.end());
// 		}
// 
// 		stable_sort(ranklist.begin(), ranklist.end(), second_rcmp());
// 		foreach (const item& r, ranklist) {
// 		    cout << queryno.at(i) << ' ' << r.first << ' ' << r.second << '\n';
// 		}
// 	    }
// 	}
//     }
//-------------------------------------------------- 

    return 0;
}
