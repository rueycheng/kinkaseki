#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <cmath>
#include <numeric>
#include <algorithm>
#include "ext/functional"

#include "magicbox.h"
#include "common.h"
#include "compat.h"

using namespace std;
using namespace magicbox;

//--------------------------------------------------
// basic_posting_list
//-------------------------------------------------- 

// General-purposed posting lists.
template<typename Key, typename Count = void> class basic_posting_list {
public:
    typedef Key key_type;
    typedef Count count_type;
    typedef pair<Key, Count> posting_type;

    basic_posting_list(): total_count(0) {}
    void add(const posting_type& p) { total_count += p.second; postings.push_back(p); }
    void add(const key_type& k, const count_type& c) { add(posting_type(k, c)); }

protected:
    count_type total_count;
    vector<posting_type> postings;
};

// Specialization for fixed-count (e.g., 1) postings
template<typename Key> class basic_posting_list<Key, void> {
public:
    typedef Key key_type;
    typedef void count_type;
    typedef Key posting_type;

    void add(const posting_type& p) { postings.push_back(p); }

protected:
    vector<posting_type> postings;
};

//--------------------------------------------------
// posting list
//--------------------------------------------------
typedef basic_posting_list<unsigned int, unsigned int> posting_list;
typedef basic_posting_list<unsigned int> fixed_count_posting_list;

// Specialized posting list for storing docnos and in-document counts
// t_i -> (D_j,f_ij)+
struct term_posting_list: public posting_list {
    unsigned int df() const { return postings.size(); }
    unsigned int tf() const { return total_count; }

    ostream& save(ostream& out) const {
	out << pack(df()) << pack(tf());
	foreach (const posting_type& p, postings) { out << pack(p.first) << pack(p.second); }
	return out;
    }
};

// Specialized posting list for storing facnos
// d_j -> (F_k)+
struct document_posting_list: public fixed_count_posting_list {
    unsigned int no_facets() const { return postings.size(); }

    const string& get_docno() const { return docno; }
    void set_docno(string s) { docno = s; }
    void set_doclen(unsigned int u) { doclen = u; }
    void set_docnorm(float f) { docnorm = f; }

    ostream& save(ostream& out) const {
	out << pack(doclen) << pack(docnorm) << pack(no_facets());
	foreach (const posting_type& p, postings) { out << pack(p); }
	return out;
    }

private:
    string docno; 
    unsigned int doclen; // for LM, okapi
    float docnorm; // for tfidf (lnc.ltc)
};

namespace std {
    ostream& operator<<(ostream& out, const term_posting_list& tpl) { return tpl.save(out); }
    ostream& operator<<(ostream& out, const document_posting_list& dpl) { return dpl.save(out); }
};

//--------------------------------------------------
// 
//-------------------------------------------------- 
template<typename T> struct writer {
    const T& subject;
    writer(const T& t): subject(t) {}
};

namespace std {
    template<typename S> ostream& operator<<(ostream& out, const writer<S>& w) { return w(out); }
};

//--------------------------------------------------
// Main program
//-------------------------------------------------- 
int main(int argc, char** argv) {

    // Getopt
    string model = "model.unnamed";

    Getopt g(argc, argv);
    g   << $(&model, "model,m", "Specify the model directory")
	<< $$$("");

    // Create model directory
    namespace fs = boost::filesystem;

    fs::create_directory(model);
    if (!fs::exists(model)) die("Cannot create directory " + model);
    fs::path basedir(model);

    // Structures
    vector<document_posting_list> documents; // A forward index + aux. info (e.g., doclen, docnorm, ..)
    unordered_map<string, term_posting_list> terms; // An inverted index

    unordered_map<string, unsigned int> facets; // An id-lookup table for facets
    vector<string> facets_fwd;

    // Now work with cin (I hate over-indented codes)
    string docno;
    string line;

    while (getline(cin, line)) {
	istringstream iss(line);

	iss >> docno;
	cerr << docno << '\n';

	if (!boost::ends_with(docno, ":facet")) {
	    // Spawn a new posting list and make it current
	    documents.push_back(document_posting_list());
	    unsigned int docid = documents.size();
	    document_posting_list& current = documents.back();

	    string t;
	    unsigned int doclen = 0;
	    float docnorm = 0.0;

	    typedef unordered_map<string, unsigned short> count_table;
	    count_table freq;
	    while (iss >> t) ++freq[t];
	    foreach (const count_table::value_type& p, freq) { terms[p.first].add(docid, p.second); }

	    current.set_docno(docno);
	    current.set_doclen(doclen);
	    current.set_docnorm(docnorm);
	}
	else {
	    // The latest is the current
	    document_posting_list& current = documents.back();

	    string f;
	    while (iss >> f) {
		if (facets.find(f) == facets.end()) {
		    facets_fwd.push_back(f);
		    facets[f] = facets_fwd.size();
		}
		current.add(facets[f]);
	    }
	}
    }

    cerr << "Done" << '\n';

    // Sort the vocabulary
    using __gnu_cxx::select1st;
    vector<string> vocab;
    transform(terms.begin(), terms.end(), back_inserter(vocab), select1st< pair<string, term_posting_list> >());
    stable_sort(vocab.begin(), vocab.end());

    // Output vocab
    {
	fs::ofstream vocab_out(basedir / "vocab");
	cerr << "Save " << basedir / "vocab" << '\n';

	vocab_out << vocab.size() << '\n';
	foreach (const string& term, vocab) { vocab_out << term << '\n'; }
    }

    // Output docno
    {
	fs::ofstream docno_out(basedir / "docno");
	cerr << "Save " << basedir / "docno" << '\n';

	docno_out << documents.size() << '\n';
	foreach (const document_posting_list& dpl, documents) { docno_out << dpl.get_docno() << '\n'; }
    }

    // Output facet
    {
	fs::ofstream facet_out(basedir / "facet");
	cerr << "Save " << basedir / "facet" << '\n';

	facet_out << facets_fwd.size() << '\n';
	foreach (const string& facet, facets_fwd) { facet_out << facet << '\n'; }
    }

    // Output t2d (binary)
    {
	fs::ofstream t2d_out(basedir / "t2d");
	cerr << "Save " << basedir / "t2d" << '\n';

	for (unsigned int i = 0; i < vocab.size(); ++i) t2d_out << pack((unsigned int)0);

	vector<unsigned int> offset;
	foreach (const string& term, vocab) {
	    const term_posting_list& tpl = terms[term];
	    t2d_out << tpl;
	    offset.push_back(t2d_out.tellp());
	}

	t2d_out.seekp(0);
	foreach (const unsigned int off, offset) { t2d_out << pack(off); }
    }

    // Output d2f (binary)
    {
	fs::ofstream d2f_out(basedir / "d2f");
	cerr << "Save " << basedir / "d2f" << '\n';

	for (unsigned int i = 0; i < documents.size(); ++i) d2f_out << pack((unsigned int)0);

	vector<unsigned int> offset;
	foreach (const document_posting_list& dpl, documents) {
	    d2f_out << dpl;
	    offset.push_back(d2f_out.tellp());
	}

	d2f_out.seekp(0);
	foreach (const unsigned int off, offset) { d2f_out << pack(off); }
    }

    return 0;
}
