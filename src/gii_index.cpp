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
// Functor
//-------------------------------------------------- 
template<typename M> struct map_lookup {
    map_lookup(const M& m, typename M::mapped_type z): m(m), z(z) {}
    typename M::mapped_type operator()(typename M::key_type key) const { 
	typename M::const_iterator p = m.find(key);
       	return (p == m.end())? z: p->second; // HACK: unforunately operator[] is never `const'
    }

private:
    const M& m;
    typename M::mapped_type z;
};

template<typename M> map_lookup<M> lookup(const M& m, typename M::mapped_type z) {
    return map_lookup<M>(m, z);
}

template<typename M> struct cmp_value_less {
    cmp_value_less(const M& m): m(m) {}
    bool operator()(typename M::size_type x, typename M::size_type y) const {
	return m[x] < m[y];
    }

private:
    const M& m;
};

template<typename M> 
cmp_value_less<M> value_less(const M& m) {
    return cmp_value_less<M>(m);
}

template<typename M> struct cmp_value_greater {
    cmp_value_greater(const M& m): m(m) {}
    bool operator()(typename M::size_type x, typename M::size_type y) const {
	return m[x] > m[y];
    }

private:
    const M& m;
};

template<typename M> 
cmp_value_greater<M> value_greater(const M& m) {
    return cmp_value_greater<M>(m);
}

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
// Subroutines
//-------------------------------------------------- 
namespace fs = boost::filesystem;

void create_model(fs::path&);
void query_model(fs::path&, bool, bool, unsigned int, unsigned int);

//--------------------------------------------------
// Main program
//-------------------------------------------------- 
int main(int argc, char** argv) {
    // Getopt
    unsigned int top_n = 1000, top_m = 100;
    string model = "model.unnamed";

    Getopt g(argc, argv);
    g   << $("query", "Query the model")
	<< $("no-result", "Do not return results (with --query)")
	<< $("no-facet", "Do not return facets (with --query)")
	<< $("silent", "Turn off error reporting")
	<< $(&top_n, "top-document,n", "Retrieve only top N results (with --query)")
	<< $(&top_m, "top-facet,m", "Retrieve only top M facet results (with --query)")
	<< $(&model, "model,m", "Specify the model directory")
	<< $$$("");

    // Create model directory
    fs::create_directory(model);
    if (!fs::exists(model)) die("Cannot open directory " + model);
    fs::path basedir(model);

    if (!g["query"]) create_model(basedir);
    else query_model(basedir, g["no-result"], g["no-facet"], top_n, top_m);

    return 0;
}

//--------------------------------------------------
// Case 1:  Create a model
//-------------------------------------------------- 
void create_model(fs::path& basedir) {
    // Structures
    vector<document_posting_list> documents; // A forward index + aux. info (e.g., doclen, docnorm, ..)
    unordered_map<string, term_posting_list> terms; // An inverted index

    unordered_map<string, unsigned int> facets; // An id-lookup table for facets
    vector<string> facets_fwd;
    unsigned int collectionlen = 0;

    // Now work with cin (I hate over-indented codes)
    string docno;
    string line;

    while (getline(cin, line)) {
	istringstream iss(line);
	iss >> docno;
	if (!boost::ends_with(docno, ":facet")) {
	    // Spawn a new posting list and make it current
	    documents.push_back(document_posting_list());
	    unsigned int docid = documents.size();
	    document_posting_list& current = documents.back();

	    cerr << boost::format("%d %s") % docid % docno << '\n';

	    string t;
	    unsigned int doclen = 0;
	    float docnorm = 0.0;

	    typedef unordered_map<string, unsigned short> count_table;
	    count_table freq;
	    while (iss >> t) ++freq[t];

	    foreach (const count_table::value_type& p, freq) { 
		terms[p.first].add(docid, p.second); 
		doclen += p.second;

		float log_count = log(p.second);
		docnorm += (1 + log_count) * (1 + log_count);
	    }

	    current.set_docno(docno);
	    current.set_doclen(doclen);
	    current.set_docnorm(docnorm);

	    collectionlen += doclen;
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

    cerr << '\n';

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
	    offset.push_back(t2d_out.tellp());
	    t2d_out << tpl;
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
	    offset.push_back(d2f_out.tellp());
	    d2f_out << dpl;
	}

	d2f_out.seekp(0);
	foreach (const unsigned int off, offset) { d2f_out << pack(off); }
    }

    // Output collection
    {
	fs::ofstream collection_out(basedir / "collection");
	cerr << "Save " << basedir / "collection" << '\n';
	collection_out << collectionlen << '\n';
    }
}

//--------------------------------------------------
// Case 2:  Query the model
//-------------------------------------------------- 
void query_model(fs::path& basedir, bool no_result, bool no_facet, unsigned int top_n, unsigned int top_m) {
    // Load vocab
    unordered_map<string, unsigned int> vocab;
    unsigned int T;

    {
	fs::ifstream vocab_in(basedir / "vocab");
	cerr << "Load " << basedir / "vocab" << '\n';

	string t;
	unsigned int id = 0;

	vocab_in >> T;
	while (vocab_in >> t) vocab[t] = ++id;
    }

    // Load docno
    unsigned int N;
    vector<string> docno;

    {
	fs::ifstream docno_in(basedir / "docno");
	cerr << "Load " << basedir / "docno" << '\n';

	string t;

	docno_in >> N;
	docno.push_back("<unk>");
	while (docno_in >> t) docno.push_back(t);
    }

    // Load facet
    unsigned int F;
    vector<string> facet;

    {
	fs::ifstream facet_in(basedir / "facet");
	cerr << "Load " << basedir / "facet" << '\n';

	string t;

	facet_in >> F;
	facet.push_back("<unk>");
	while (facet_in >> t) facet.push_back(t);
    }

    // Score accumulator (`brute-force')
    vector<float> score;
    score.resize(N + 1);

    vector<unsigned short> facet_count;
    facet_count.resize(F + 1);

    vector<float> facet_rank;
    vector<float> facet_norm;
    facet_rank.resize(F + 1);
    facet_norm.resize(F + 1);

    // Load t2d offset
    fs::ifstream t2d_in(basedir / "t2d");
    vector<unsigned int> t2d_offset;
    t2d_offset.resize(T + 1);
    cerr << "Load " << basedir / "t2d" << " (offsets)" << '\n';
    for (unsigned int i = 1; i <= T; ++i) t2d_in >> unpack(t2d_offset[i]);

    // Load d2f offset
    fs::ifstream d2f_in(basedir / "d2f");
    vector<unsigned int> d2f_offset;
    d2f_offset.resize(N + 1);
    cerr << "Load " << basedir / "d2f" << " (offsets)" << '\n';
    for (unsigned int i = 1; i <= N; ++i) d2f_in >> unpack(d2f_offset[i]);

    // Load dlen
    vector<unsigned int> dlen;
    dlen.resize(N + 1);
    cerr << "Load " << basedir / "d2f" << " (dlen)" << '\n';
    for (unsigned int i = 1; i <= N; ++i) {
	d2f_in.seekg(d2f_offset[i]);
	d2f_in >> unpack(dlen[i]);
    }

    // Load collection
    unsigned int L;
    
    {
	fs::ifstream collection_in(basedir / "collection");
	collection_in >> L;
    }

    // HACK: They said programmers are lazy people
    const float mu = 2500.0;

    // Now, we are ready to work with queries
    unsigned int topic_id = 0;
    string topicno;
    string line;

    while (getline(cin, line)) {
	istringstream iss(line);
	iss >> topicno;
	++topic_id;

	// Make sure we had a good topic number
	vector<unsigned int> query;

	if (boost::ends_with(topicno, ":topic")) strip_after_first(topicno, ':');
	else {
	    query.push_back(lookup(vocab, 0)(topicno)); // Otherwise, throw it back into the query
	    topicno = str(boost::format("%03d") % topic_id);
	}

	// Read the rest of the query
	transform(istream_iterator<string>(iss), istream_iterator<string>(), 
		back_inserter(query), lookup(vocab, 0));

	// Reset the accumulators
	fill(score.begin(), score.end(), 0.0);

	//--------------------------------------------------
	// Tthe full-blown scoring function
	//
	// \log \Pr(Q|d) = \sum_{q \in Q \cap d} \log(1 + f_{q,d} |C| / \mu f_q)    ...(1)
	//               + |Q| \log \mu - |Q| \log |C| + \sum_{q \in Q} \log f_q    ...(2)
	//               - |Q| \log(|d| + \mu)                                      ...(3)
	//
	// NOTE: (1) is done in posting-list traversal,
	//       (2) is a globally-determined constant in the session,
	//       and (3) is a document-dependent normalizing factor.
	//--------------------------------------------------
	
	// Placeholder for the result of (2)
	float sum_of_logtf = 0.0;

	// Iterate through each term
	foreach (unsigned int term_id, query) {
	    if (term_id == 0) continue;

	    // Make a jump
	    t2d_in.seekg(t2d_offset[term_id]);

	    unsigned int df, tf;
	    t2d_in >> unpack(df) >> unpack(tf);

	    // Accumulate counts for computation of (2)
	    sum_of_logtf += log(tf);

	    for (unsigned int i = 0; i < df; ++i) {
		unsigned int doc_id, count;
		t2d_in >> unpack(doc_id) >> unpack(count);

		// Compute (1)
		score[doc_id] += log(1 + float(count * L) / (mu * tf));
	    }
	}

	// Actually, (2) is not necessary here (even in the later facet ranking)
	unsigned int qlen = query.size();

	// Collect non-zero id's
	vector<unsigned int> rank;

	for (unsigned int i = 1; i <= N; ++i) {
	    if (!score[i]) continue;
	    rank.push_back(i);


	    // Add (3) back in
	    score[i] -= qlen * log(dlen[i] + mu);
	}

	// Make a copy for faster facet lookup
	vector<unsigned int> copied;

	if (top_n == 0 || rank.size() <= top_n) {
	    // Now copy as it is
	    copy(rank.begin(), rank.end(), back_inserter(copied));
	}
	else {
	    // HACK: Now we need to copy the top N stuff back into candidates, in ascending order
	    nth_element(rank.begin(), rank.begin() + top_n, rank.end(), value_greater(score));
	    rank.erase(rank.begin() + top_n, rank.end());
	    copy(rank.begin(), rank.end(), back_inserter(copied));
	    stable_sort(copied.begin(), copied.end());
	}

	// NOTE: I made it a little bit weird here.  Should subject to change.
	if (!no_result) {
	    stable_sort(rank.begin(), rank.end(), value_greater(score));
	    foreach (unsigned int doc_id, rank) {
		cout << topicno << ' ' << docno[doc_id] << ' ' << score[doc_id] << '\n';
	    }
	}

	// Now, produce facets if you will
	if (no_facet) return;
	
	// Reset facet_counts
	fill(facet_count.begin(), facet_count.end(), 0);
	fill(facet_rank.begin(), facet_rank.end(), 0.0);
	fill(facet_norm.begin(), facet_norm.end(), 0.0);

	foreach (unsigned int doc_id, copied) {
	    unsigned int nop1, ff; // Does `facet frequency' sound weird to you?
	    float nop2;
	    d2f_in.seekg(d2f_offset[doc_id]);
	    d2f_in >> unpack(nop1) >> unpack(nop2) >> unpack(ff);
	    for (unsigned int i = 0; i < ff; ++i) {
		unsigned int facet_id;
		d2f_in >> unpack(facet_id);

		//--------------------------------------------------
		// Phase 1: counts
		//-------------------------------------------------- 
		++facet_count[facet_id];

		 //--------------------------------------------------
		// Phase 2: rank
		//-------------------------------------------------- 
		facet_rank[facet_id] += exp(score[doc_id]) * (1.0 / ff);
		facet_norm[facet_id] += (1.0 / ff);
	    }
	}

	//--------------------------------------------------
	// Output results for counts
	//-------------------------------------------------- 
	vector<unsigned int> facet_candidate;

	for (unsigned int i = 1; i <= F; ++i) {
	    if (!facet_count[i]) continue;
	    facet_candidate.push_back(i);
	}

	if (top_m != 0 && facet_candidate.size() > top_m) {
	    nth_element(facet_candidate.begin(), facet_candidate.begin() + top_m,
		    facet_candidate.end(), value_greater(facet_count));
	    facet_candidate.erase(facet_candidate.begin() + top_m, facet_candidate.end());
	}

	stable_sort(facet_candidate.begin(), facet_candidate.end(), value_greater(facet_count));
	foreach (unsigned int facet_id, facet_candidate) {
	    cout << topicno << ":facet-count" << ' ' << facet[facet_id] << ' ' << facet_count[facet_id] << '\n';
	}
	
	//--------------------------------------------------
	// Do it all over again for ranks
	//-------------------------------------------------- 
	facet_candidate.clear();

	for (unsigned int i = 1; i <= F; ++i) {
	    if (!facet_norm[i]) continue;
	    if (facet_count[i] < 2) continue;
	    facet_candidate.push_back(i);
	    facet_rank[i] = log(facet_rank[i]) - log(facet_norm[i]);
	}

	if (top_m != 0 && facet_candidate.size() > top_m) {
	    nth_element(facet_candidate.begin(), facet_candidate.begin() + top_m,
		    facet_candidate.end(), value_greater(facet_rank));
	    facet_candidate.erase(facet_candidate.begin() + top_m, facet_candidate.end());
	}

	stable_sort(facet_candidate.begin(), facet_candidate.end(), value_greater(facet_rank));
	foreach (unsigned int facet_id, facet_candidate) {
	    cout << topicno << ":facet-rank" << ' ' << facet[facet_id] << ' ' << facet_rank[facet_id] << ' ' << facet_count[facet_id] << '\n';
	}
    }
}
