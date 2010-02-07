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
// Vocabulary
//-------------------------------------------------- 
template<typename __word = unsigned int>
class vocabulary {
    unordered_map<string, __word> vocab;
    vector<string> rvocab;
    vector<unsigned int> counts;

public:
    vocabulary() { 
	rvocab.push_back("<unk>"); 
	counts.push_back(0);
    }

    vocabulary(istream& in, unsigned int vocab_size) {
	vocabulary();
	string line;

	while (vocab_size-- > 0 && getline(in, line)) {
	    istringstream iss(line);
	    string word;
	    unsigned int count = 0;
	    iss >> word >> count;

	    rvocab.push_back(word);
	    counts.push_back(count);
	    vocab[word] = rvocab.size() - 1;
	}
    }

    void count(__word id) { ++counts.at(id); }
    unsigned int get_count(__word id) { return counts.at(id); }

    void save(ostream& o) { 
	vector<string>::const_iterator w_iter = rvocab.begin();
	vector<unsigned int>::const_iterator c_iter = counts.begin();

	while (w_iter != rvocab.end())
	    o << *w_iter++ << ' ' << *c_iter++ << "\n";
    }

    __word encode(const string& word) {
	if (vocab.find(word) != vocab.end()) return vocab[word];
	rvocab.push_back(word);
	counts.push_back(0);
	return vocab[word] = rvocab.size() - 1;
    }

    string decode(unsigned int code) { return rvocab.at(code); }
    __word size() { return rvocab.size(); }
    unsigned int total() { return accumulate(counts.begin(), counts.end(), 0); }
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
    vector<string> facets;
};

namespace std {
    ostream& operator<<(ostream& o, const posting& p) {
	o << p.name << ' ' << p.counts.size() << ' ' << p.total << ' ';
	foreach (const posting::count& c, p.counts) { o << c.first << ' ' << c.second << ' '; }
	foreach (const string& f, p.facets) { o << f << ' '; }
	return o;
    }
}

typedef posting term_posting;
typedef posting document_posting;

//--------------------------------------------------
// Main program
//-------------------------------------------------- 
int main(int argc, char** argv) {

    // Getopt
    vector<string> input;
    float mu = 5000;
    string method = "lm";

    Getopt g(argc, argv);
    g   << $("test", "Test the model")
	<< $("index-sentences", "Index sentences rather than documents")
	<< $("retrieve-facets", "Retrieve facets")
	<< $(&method, "method,m", "Ranking method: tfidf or lm")
	<< $(&mu, "mu", "Specify the parameter 'mu' (for LM)")
	<< $(&input, "input", "@@ IT DOES NOT MATTER WHAT I PUT IN HERE @@", -1)
	<< $$$("files..");

    namespace fs = boost::filesystem;

    if (!g["test"]) {
	if (input.empty()) input.push_back("-");

	// Collection-wide
	vocabulary<> vocab;
	vector<document_posting> documents;

	// Current documents
	string docno;
	unsigned int total = 0;
	unordered_map<unsigned int, unsigned short> counts;
	vector<string> facets;

	// Tentatively ignore '--index-sentences' and facets
	if (g["index-sentences"]) throw NotSupportedFunction();

	foreach (const string& filename, input) {
	    AutoIn in(filename);
	    string line;

	    while (getline(in(), line)) {
		if (line.empty()) {
		    documents.push_back(document_posting());
		    document_posting& dest = documents.back();

		    dest.name = docno;
		    dest.total = total;
		    copy(counts.begin(), counts.end(), back_inserter(dest.counts));
		    stable_sort(dest.counts.begin(), dest.counts.end(), first_cmp());
		    copy(facets.begin(), facets.end(), back_inserter(dest.facets));

		    docno.clear();
		    total = 0;
		    counts.clear();
		    facets.clear();
		    continue;
		}

		string id;
		vector<string> words;
		istringstream iss(line);

		iss >> id;
		if (!boost::ends_with(id, ":facet")) {
		    strip_after_first(id, '#'); // Lose sentence no.
		    if (docno.empty()) docno = id;
		    if (docno != id) throw InvalidFormat();

		    copy(istream_iterator<string>(iss), 
			    istream_iterator<string>(), back_inserter(words));

		    foreach (const string& word, words) {
			unsigned int word_id = vocab.encode(word);
			++counts[word_id];
			++total;
			vocab.count(word_id);
		    }
		}
		else {
		    strip_after_first(id, ':'); // Lose the facet mark
		    if (docno.empty()) docno = id;
		    if (docno != id) throw InvalidFormat();

		    copy(istream_iterator<string>(iss), 
			    istream_iterator<string>(), back_inserter(facets));

		    stable_sort(facets.begin(), facets.end());
		}
	    }

	    if (!docno.empty()) throw InvalidFormat();
	}

	cout << "@vocab " << vocab.size() << ' ' << vocab.total() << "\n";
	vocab.save(cout);

	cout << "@document " << documents.size() << "\n";
	foreach (const posting& dp, documents) { cout << dp << "\n"; }
    }
    else {
	if (input.empty()) throw MissingArgument();

	string line;

	// Get query
	vector<string> query;
	
	{
	    ifstream in(input.front().c_str());
	    if (!getline(in, line)) throw InvalidFormat();
	    istringstream iss(line);

	    if (!boost::starts_with(line, "@query ")) throw InvalidFormat();
	    strip_before_first(line, ' ');

	    PorterStemmer ps;
	    CJKVTokenizer tok(line);
	    foreach (string t, tok) { 
		if ((t[0] & 0x80) == 0x00) {
		    boost::to_lower(t);
		    query.push_back(ps.stem(t.c_str(), t.size())); 
		} else query.push_back(t);
	    }
	}

	// Get vocabulary
	unsigned int vocab_size, vocab_total;

	{
	    if (!getline(cin, line)) throw InvalidFormat();
	    istringstream iss(line);

	    string head;
	    if (!(iss >> head >> vocab_size >> vocab_total) || head != "@vocab") 
		throw InvalidFormat();
	}

	// Load vocabulary (words and counts)
	vocabulary<> vocab(cin, vocab_size);

	// Get the number of documents
	unsigned int document_size;

	{
	    if (!getline(cin, line)) throw InvalidFormat();
	    istringstream iss(line);

	    string head;
	    if (!(iss >> head >> document_size) || head != "@document") 
		throw InvalidFormat();
	}

	// Now translate the query
	float query_normalization;
	vector<posting::weight> query_vector;
	float logN = log(document_size);

	foreach (const string& query_word, query) { 
	    unsigned int query_word_id = vocab.encode(query_word);
	    unsigned int query_word_df = vocab.get_count(query_word_id);
	    float query_weight = log(2) * (logN - log(query_word_df));
	    query_vector.push_back(posting::weight(query_word_id, query_weight));
	}

	stable_sort(query_vector.begin(), query_vector.end(), first_cmp());

	{
	    float score_part = 0;
	    foreach (const posting::weight& c, query_vector) { score_part += c.second * c.second; }
	    query_normalization = sqrt(score_part);
	}

	// NOTE: A bloody mess.
	typedef unordered_map<string, unsigned short> bag;

	vector<bag> query_word_facets;
	bag query_facets;

	for (unsigned int i = 0; i < query_vector.size(); ++i)
	   query_word_facets.push_back(bag());

	// Get documents
	for (unsigned int i = 0; i < document_size; ++i) {
	    if (!getline(cin, line)) throw InvalidFormat();
	    istringstream iss(line);

	    string docno;
	    unsigned int docsize, total;
	    if (!(iss >> docno >> docsize >> total)) throw InvalidFormat();

	    float document_normalization;
	    vector<posting::weight> document_vector;

	    for (unsigned int j = 0; j < docsize; ++j) {
		unsigned int word_id;
		short count;
		if (!(iss >> word_id >> count)) throw InvalidFormat();

		document_vector.push_back(posting::weight(word_id, log(1 + count)));
	    }

	    {
		float score_part = 0;
		foreach (const posting::weight& c, document_vector) { score_part += c.second * c.second; }
		document_normalization = sqrt(score_part);
	    }

	    // Load facets
	    vector<string> facets;
	    copy(istream_iterator<string>(iss), 
		    istream_iterator<string>(), back_inserter(facets));

	    // Dot
	    vector<posting::weight>::const_iterator q_iter = query_vector.begin();
	    vector<posting::weight>::const_iterator q_end = query_vector.end();
	    vector<posting::weight>::const_iterator d_iter = document_vector.begin();
	    vector<posting::weight>::const_iterator d_end = document_vector.end();

	    vector<bag>::iterator qf_iter = query_word_facets.begin(); // The root of evil!

	    unsigned int sum = 0;
	    while (q_iter != q_end && d_iter != d_end) {
		bool larger = first_cmp()(*d_iter, *q_iter);
		bool smaller = first_cmp()(*q_iter, *d_iter);

		if (larger) ++d_iter;
		else if (smaller) {
		    ++q_iter; ++qf_iter;
		}
		else {
		    sum += (*d_iter).second * (*q_iter).second;
		    foreach (const string& facet, facets) { 
			++(*qf_iter)[facet]; 
			++query_facets[facet];
		    }

		    ++d_iter;
		    ++q_iter;
		    ++qf_iter;
		}
	    }

	    float score = sum / document_normalization;

	    cout << docno << ' ' << boost::format("%0.6f") % score << "\n";
	}

	cout << "__FACETS__" << "\n";
	foreach (const bag& b, query_word_facets) {
	    foreach (const bag::value_type& v, b) { cout << v.first << ' ' << v.second << ' '; }
	    cout << "\n";
	}

	foreach (const bag::value_type& v, query_facets) { cout << v.first << ' ' << v.second << ' '; }
	cout << "\n";
    }

    return 0;
}
