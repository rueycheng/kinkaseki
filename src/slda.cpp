#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <algorithm>
#include <iterator>
#include <ctime>
#include <boost/random.hpp>

#include "magicbox.h"
#include "common.h"
#include "compat.h"

using namespace std;
using namespace magicbox;

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
// RNG
//-------------------------------------------------- 
template<typename __float = float, 
    class __Generator = boost::mt19937, 
    class __Distribution = boost::uniform_01<__Generator, __float> >
class RNG {
    __Distribution d;

public:
    RNG(uint32_t seed = std::time(0)): d(__Generator(seed)) {}
    template<typename __number>
	__number operator()() { return __number(d()); }
    template<typename __number>
	__number operator()(__number n) { return __number(n * d()); }
};

//--------------------------------------------------
// SLDAGibbsSampler
//-------------------------------------------------- 
template<class __SLDAModel, typename __float = long double, typename __count = unsigned int>
class SLDAGibbsSampler {
public:
    typedef typename __SLDAModel::size_type size_type;
    typedef typename __SLDAModel::doc_type doc_type;
    typedef typename __SLDAModel::sent_type sent_type;
    typedef typename __SLDAModel::topic_type topic_type;
    typedef typename __SLDAModel::word_type word_type;
    typedef __SLDAModel slda_type;
    typedef __float prob_type;
    typedef __count count_type;

private:
    slda_type& m;
    prob_type alpha, beta, delta;

    SLDAGibbsSampler(SLDAGibbsSampler&) {} // No copy

    // Counts for the current model
    count_type* n_dx; // document-s_topic, M by S
    count_type* n_d_; // document-*, M
    count_type* n_xz; // s_topic-w_topic, S by T
    count_type* n_x_; // s_topic-*, S
    count_type* n_zw; // topic-word, T by W
    count_type* n_z_; // topic-*, T

    prob_type* sample_x; // cdf for p(x|x*), S + 1
    prob_type* sample_z; // cdf for p(z|z*), T + 1

    // Model initialization
    void alloc(doc_type m, topic_type s, topic_type t, word_type w) {
	n_dx = new count_type[m * s];
	n_d_ = new count_type[m];
	n_xz = new count_type[s * t];
	n_x_ = new count_type[s];
	n_zw = new count_type[t * w];
	n_z_ = new count_type[t];
	fill_n(n_dx, m * s, 0);
	fill_n(n_d_, m, 0);
	fill_n(n_xz, s * t, 0);
	fill_n(n_x_, s, 0);
	fill_n(n_zw, t * w, 0);
	fill_n(n_z_, t, 0);

	sample_x = new prob_type[s + 1];
	sample_z = new prob_type[t + 1];
    }

    void free() {
	delete[] n_dx; 
	delete[] n_d_; 
	delete[] n_xz; 
	delete[] n_x_; 
	delete[] n_zw; 
	delete[] n_z_; 

	delete[] sample_x;
	delete[] sample_z;
    }

    void collectCount(slda_type& m) {
	for (size_type i = 0; i < m.N; ++i) {
	    const topic_type& x = m.x[i];
	    const topic_type& z = m.z[i];
	    const word_type& w = m.w[i];
	    ++n_xz[x * m.T + z];
	    ++n_x_[x];
	    ++n_zw[z * m.W + w];
	    ++n_z_[z];
	}

	for (size_type i = 0; i < m.N; ) {
	    const doc_type& d = m.d[i];
	    const sent_type& s = m.s[i];
	    const topic_type& x = m.x[i];
	    ++n_dx[d * m.S + x];
	    ++n_d_[d];

	    while (m.s[++i] == s) ;
	}
    }

    void collectCount(slda_type& m, slda_type& refm) {
	collectCount(m);
	for (size_type i = 0; i < refm.N; ++i) {
	    const topic_type& x = refm.x[i];
	    const topic_type& z = refm.z[i];
	    const word_type& w = refm.w[i];
	    ++n_xz[x * m.T + z];
	    ++n_x_[x];
	    ++n_zw[z * m.W + w];
	    ++n_z_[z];
	}
    }

public:
    // One-model initialization
    SLDAGibbsSampler(slda_type& m, prob_type alpha, prob_type beta, prob_type delta): 
	m(m), alpha(alpha), beta(beta), delta(delta) 
    {
	this->alpha = (alpha == 0.0? 50.0/m.T: alpha);
	alloc(m.M, m.S, m.T, m.W);
	collectCount(m);
    }

    // Two-model initialization (`test')
    SLDAGibbsSampler(slda_type& m, slda_type& refm, prob_type alpha, prob_type beta, prob_type delta):
	m(m), alpha(alpha), beta(beta), delta(delta)
    {
	this->alpha = (alpha == 0.0? 50.0/m.T: alpha); 
	alloc(m.M, m.S, m.T, m.W);
	collectCount(m, refm);
    }

    void update(RNG<prob_type>& rng, prob_type* ppl = 0) {
	const prob_type W_beta = m.W * beta;
	const prob_type T_delta = m.T * delta;
	const prob_type S_alpha = m.S * alpha;

	if (ppl) *ppl = 0.0;

	// Pass 1: Iterate over words
	for (size_type i = 0; i < m.N; ++i) {
	    const topic_type& x = m.x[i];
	    const topic_type& z = m.z[i]; // We might change this value later (non-const)
	    const word_type& w = m.w[i];

	    // 1) Decrement.  Like tuple (x,z,w) never exists.
	    const size_type xT = x * m.T;

	    --n_xz[xT + z];
	    --n_x_[x];
	    --n_zw[z * m.W + w];
	    --n_z_[z];

	    // 2) Collect conditional probabilities
	    prob_type prob_of_sum = 0.0;
	    prob_type f2d = n_x_[x] + T_delta;
	    for (topic_type zz = 0; zz < m.T; ++zz) {
		prob_type f2n = n_xz[xT + zz] + delta;
		prob_type f1n = n_zw[zz * m.W + w] + beta;
		prob_type f1d = n_z_[zz] + W_beta; 
		sample_z[zz] = (prob_of_sum += (f1n / f1d) * (f2n / f2d)); // Also obtain CDF
	    }

	    // NOTE: Sometimes old-fashion way works best
	    //if (ppl) *ppl += log(prob_of_sum); 

	    // 3) Toss the dice and get the sample
	    prob_type toss = rng(sample_z[m.T - 1]);

	    topic_type zz = 0;
	    while (toss > sample_z[zz]) ++zz;
	    m.z[i] = zz; // Gotcha!

	    // 4) Increment
	    ++n_xz[xT + zz];
	    ++n_x_[x];
	    ++n_zw[zz * m.W + w];
	    ++n_z_[zz];
	}

	// Pass 2: Iterate over sentences
	for (size_type i = 0; i < m.N; ) {  // Empty incrementor
	    const doc_type& d = m.d[i];
	    const sent_type& s = m.s[i]; // We might change this later
	    const topic_type& x = m.x[i];

	    // Find the end of the sentence
	    size_type i_end = i + 1;
	    while (i_end < m.N && m.s[i_end] == s) ++i_end;

	    const size_type dS = d * m.S;
	    const size_type xT = x * m.T;

	    //--------------------------------------------------
	    // NOTE: This is the double increment-decrement trick!
	    //-------------------------------------------------- 
	    // 1) Decrement all the counts made to words [i, i_end)
	    const size_type s_size = i_end - i;
	    for (size_type j = i; j < i_end; ++j) --n_xz[xT + m.z[j]];
	    n_x_[x] -= s_size;
	    --n_dx[dS + x];
	    --n_d_[d];

	    // 2) Calculate CDF and run the sampler
	    prob_type prob_of_sum = 0.0;
	    prob_type f2d = n_d_[d] + S_alpha;
	    for (topic_type xx = 0; xx < m.S; ++xx) {
		prob_type f2n = n_dx[dS + xx] + alpha;
		prob_type f1 = 1.0;

		// 2.1) Increment the counts back by pretending the sentence got assigned as topic 'xx'
		const size_type xxT = xx * m.T;
		for (size_type j = i; j < i_end; ++j)
		    f1 *= (n_xz[xxT + m.z[j]]++ + delta) / (n_x_[xx]++ + T_delta);

		// 2.2) Decrement
		for (size_type j = i; j < i_end; ++j) --n_xz[xxT + m.z[j]];
		n_x_[xx] -= s_size;

		// 2.3) Get sample
		sample_x[xx] = (prob_of_sum += f1 * (f2n / f2d)); // Also obtain CDF
	    }

	    // if (ppl) *ppl += log(prob_of_sum); 

	    // 3) Toss the dice and get the sample
	    prob_type toss = rng(sample_x[m.S - 1]);

	    topic_type xx = 0;
	    while (toss > sample_x[xx]) ++xx;

	    for (size_type j = i; j < i_end; ++j) m.x[j] = xx; // Gotcha!

	    // 4) Increment
	    const size_type xxT = xx * m.T;
	    for (size_type j = i; j < i_end; ++j) ++n_xz[xxT + m.z[j]];
	    n_x_[xx] += s_size;
	    ++n_dx[dS + xx];
	    ++n_d_[d];

	    // Done with this sentence
	    i += s_size;
	}

	// Compute the perplexity (Can't figure out a better way to do this)
	if (ppl) {
	    *ppl = 0.0;

	    for (size_type i = 0; i < m.N; ++i) {
		const doc_type& d = m.d[i];
		const word_type& w = m.w[i];

		const size_type dS = d * m.S;
		prob_type f_theta_d = n_d_[d] + S_alpha;

		prob_type sum = 0.0;
		for (topic_type xx = 0; xx < m.S; ++xx) {
		    prob_type f_theta_n = n_dx[dS + xx] + alpha;
		    prob_type f_tau_d = n_x_[xx] + T_delta;

		    const size_type xxT = xx * m.T;

		    prob_type partial_sum = 0.0;
		    for (topic_type zz = 0; zz < m.T; ++zz) {
			prob_type f_tau_n = n_xz[xxT + zz] + delta;
			prob_type f_phi_n = n_zw[zz * m.W + w] + beta;
			prob_type f_phi_d = n_z_[zz] + W_beta;

			partial_sum += (f_tau_n / f_tau_d) * (f_phi_n / f_phi_d);
		    }

		    sum += (f_theta_n / f_theta_d) * partial_sum;
		}

		*ppl += log(sum);
	    }

	    *ppl = exp(-*ppl / m.N);
	}
    }

    void outputWordCluster(ostream &o, vector<string>& vocab, unsigned int top_n = 10) {
	const prob_type W_beta = m.W * beta;
	for (topic_type z = 0; z < m.T; ++z) {
	    const size_type zW = z * m.W;
	    const prob_type denom = n_z_[z] + W_beta;

	    typedef pair<word_type, prob_type> item;
	    vector<item> rank;
	    for (word_type w = 0; w < m.W; ++z) rank.push_back(item(w, (n_zw[zW + w] + beta) / denom));

	    if (rank.size() > top_n) nth_element(rank.begin(), rank.begin() + top_n, rank.end(), second_rcmp());
	    rank.erase(rank.begin() + top_n, rank.end());
	    stable_sort(rank.begin(), rank.end(), second_rcmp());
	    o << z << ' ';
	    foreach (const item& r, rank) { o << vocab.at(r.first) << ' '; }
	    o << '\n';
	}
    }

    void outputTheta(ostream& o, vector<string>& docno) {
	//--------------------------------------------------
	// const prob_type T_alpha = m.T * alpha;
	// const doc_type M = m.M;
	// const topic_type T = m.T;
	// bool use_docno = docno.size() == M;
	//-------------------------------------------------- 

	//--------------------------------------------------
	// for (doc_type d = 0; d < M; ++d) {
	//     if (use_docno) o << docno.at(d) << ' ';
	//-------------------------------------------------- 

	//--------------------------------------------------
	//     prob_type denom = T_alpha + n_d_[d];
	//     for (topic_type z = 0; z < T; ++z) o << (alpha + n_dz[d * T + z]) / denom << ' ';
	//     o << "\n";
	// }
	//-------------------------------------------------- 
    }

    void outputPhiTheta(ostream& o, vector<string>& docno) {
	//--------------------------------------------------
	// const prob_type W_beta = m.W * beta;
	// const prob_type T_alpha = m.T * alpha;
	// const doc_type M = m.M;
	// const word_type W = m.W;
	// const topic_type T = m.T;
	// bool use_docno = docno.size() == M;
	//-------------------------------------------------- 

	//--------------------------------------------------
	// for (doc_type d = 0; d < M; ++d) {
	//     if (use_docno) o << docno.at(d) << ' ';
	//-------------------------------------------------- 

	//--------------------------------------------------
	//     prob_type theta_denom = T_alpha + n_d_[d];
	//     for (word_type w = 0; w < W; ++w) {
	// 	prob_type sum = 0;
	// 	for (topic_type z = 0; z < T; ++z) {
	// 	    prob_type theta_dz = (alpha + n_dz[d * T + z]) / theta_denom;
	// 	    prob_type phi_zw = (beta + n_zw[z * W + w]) / (W_beta + n_z_[z]);
	// 	    sum += theta_dz * phi_zw;
	// 	}
	// 	o << sum << ' ';
	//     }
	//     o << "\n";
	// }
	//-------------------------------------------------- 
    }

//--------------------------------------------------
//     Adding note here
//-------------------------------------------------- 
//     void outputTheta(ostream& o) {
// 	__R T_alpha = T * alpha;
// 	o << alpha << ' ' << M << ' ' << T << endl; // The first line: alpha, M, and T
// 	for (__I j = 0; j < M; ++j) {
// 	    __R denom = T_alpha + n_dx[j];
// 	    for (short k = 0; k < T; ++k) o << (alpha + n_dt[j * T + k]) / denom << ' ';
// 	    o << endl;
// 	}
//     }
//     
//     void outputPhi(ostream& o) {
// 	__R W_beta = W * beta;
// 	o << beta << ' ' << T << ' ' << W << endl; // The first line: beta, T, and W
// 	for (short k = 0; k < T; ++k) {
// 	    __R denom = W_beta + n_xt[k];
// 	    for (__I v = 0; v < W; ++v) o << (beta + n_wt[v * T + k]) / denom << ' ';
// 	    o << endl;
// 	}
//     }

    ~SLDAGibbsSampler() { free(); }
};

//--------------------------------------------------
// SLDAModel
//-------------------------------------------------- 
template<typename __size = unsigned int,
    typename __doc = unsigned int, 
    typename __sent = unsigned int,
    typename __topic = unsigned short,
    typename __word = unsigned int>
class SLDAModel {
public:
    typedef __size size_type;
    typedef __doc doc_type;
    typedef __sent sent_type;
    typedef __topic topic_type;
    typedef __word word_type;

private:
    size_type N;
    doc_type M; 
    sent_type L;
    topic_type S;
    topic_type T;
    word_type W; 
    doc_type *d; // sized N, each in the range of [0, M)
    sent_type *s; // size N, each in the range of [0, L)
    topic_type *x; // sized N, each in the range of [0, S)
    topic_type *z; // sized N, each in the range of [0, T)
    word_type *w; // sized N, each in the range of [0, W)

    SLDAModel(SLDAModel& model) {} // No copy

    void alloc(size_type n) { 
	d = new doc_type[n]; 
	s = new sent_type[n]; 
	x = new topic_type[n]; 
	z = new topic_type[n]; 
	w = new word_type[n]; 
    }

    void free() { delete[] d; delete[] s; delete[] x; delete[] z; delete[] w; }

public:
    template<typename prob_type>
	SLDAModel(size_type N, doc_type M, sent_type L, topic_type S, topic_type T, word_type W, vector<word_type>& wseq, 
	    vector<doc_type>& doc, vector<sent_type>& sent, RNG<prob_type>& rng): N(N), M(M), L(L), S(S), T(T), W(W) 
    { 
	alloc(N);

	typename vector<word_type>::const_iterator witer = wseq.begin();
	typename vector<doc_type>::const_iterator diter = doc.begin();
	typename vector<sent_type>::const_iterator siter = sent.begin();
	for (size_type i = 0; i < N; ++i) {
	    d[i] = *diter++;
	    s[i] = *siter++;
	    x[i] = rng(S);
	    z[i] = rng(T);
	    w[i] = *witer++;
	}
    }

    explicit SLDAModel(istream& seq) { 
	seq >> N >> M >> L >> S >> T >> W;
	alloc(N);

	for (size_type i = 0; i < N; ++i) seq >> d[i] >> s[i] >> x[i] >> z[i] >> w[i];
    }

    void save(ostream& seq) {
	seq << N << ' ' << M << ' ' << L << ' ' << S << ' ' << T << ' ' << W << "\n";
	for (size_type i = 0; i < N; ++i) seq << d[i] << ' ' << s[i] << ' ' << x[i] << ' ' << z[i] << ' ' << w[i] << "\n";
    }

    ~SLDAModel() { free(); }

    size_type getN() { return N; }
    doc_type getD() { return M; }
    sent_type getL() { return L; }
    topic_type getS() { return S; }
    topic_type getT() { return T; } 
    word_type getW() { return W; }

    friend class SLDAGibbsSampler< SLDAModel<size_type, doc_type, sent_type, topic_type, word_type> >;
};

//--------------------------------------------------
// Utility
//-------------------------------------------------- 
void bye(Getopt& g) {
    cerr << g;
    exit(0);
}

//--------------------------------------------------
// Main program
//-------------------------------------------------- 
int main(int argc, char** argv) {
    // Getopt stuff
    Getopt g(argc, argv);

    float alpha = 0;
    float beta = 0.01;
    float delta = 0.01;
    unsigned int nwtopic = 100;
    unsigned int nstopic = 100;
    unsigned int niter = 1000;
    unsigned int interval = 1;
    unsigned int sample_num = 1;
    unsigned int sample_lag = 50;
    unsigned int top_n = 100;
    string model = "model.unnamed";
    string output = "none";
    vector<string> arg;

    g   << $("train", "Run in the training mode")
	<< $("test", "Run in the test mode")
	<< $("test-cohesion", "Run in the test mode (for testing cohesion)")
	<< $("dump", "Run in the dump mode")
	<< $(&output, "output", "Specify the output: 'none', 'phi_theta'")
	<< $(&model, "model,m", "Path to the model")
	<< $(&alpha, "alpha,a", "Specify the parameter 'alpha'")
	<< $(&beta, "beta,b", "Specify the parameter 'beta'")
	<< $(&delta, "delta,d", "Specify the parameter 'delta'")
	<< $(&nwtopic, "nwtopic,T", "Specify the number of word topics")
	<< $(&nstopic, "nstopic,S", "Specify the number of sentence topics")
	<< $(&niter, "niter,I", "Specify the number of runs")
	<< $(&interval, "interval", "Specify the interval (#iter) between perplexity reports")
	<< $(&sample_num, "sample-num", "Specify the number of read-outs (after burned-in)")
	<< $(&sample_lag, "sample-lag", "Specify the number of iteration between each read-out")
	<< $(&top_n, "top-n", "Show only top N results")
	<< $(&arg, "arg", "", -1)
	<< $$$("lda {--train|--test|--test-cohesion} [options..] [input-files..]");

    // 0) Prepare utilities
    RNG<long double> rng;

    if (!g["train"] && !g["test"] && !g["test-cohesion"] && !g["dump"]) bye(g);
    if (arg.empty()) arg.push_back("-");

    if (g["train"]) {
	// 1) Scan input
	Vocabulary<> vocab;
	vector<unsigned int> wseq;
	vector<unsigned int> doc;
	vector<unsigned int> sent;
	vector<string> docno;

	unsigned int docno_id = 0;
	unsigned int sentno_id = 0;
	string last_doc;

	foreach (const string& filename, arg) {
	    AutoIn in(filename);
	    string words, word;

	    while (getline(in(), words)) {
		istringstream wstream(words);

		wstream >> word;
		docno.push_back(word);

		strip_after_first(word, ':');
		if (last_doc != word) {
		    if (!last_doc.empty()) ++docno_id;
		    last_doc = word;
		}

		while (wstream >> word) {
		    wseq.push_back(vocab.encode(word));
		    doc.push_back(docno_id);
		    sent.push_back(sentno_id);
		}
		++sentno_id;
	    }
	}
	++docno_id;

	SLDAModel<> training_set(wseq.size(), docno_id, sentno_id, nstopic, nwtopic, vocab.size(), wseq, doc, sent, rng);
	SLDAGibbsSampler<SLDAModel<> > sampler(training_set, alpha, beta, delta);

	// NOTE: Clean up memory space explicitly
	wseq.clear(); 
	doc.clear(); 

	// Run the sampler
	for (unsigned int iter = 1; iter <= niter; ++iter) {
	    if (interval && iter % interval == 0) {
		long double ppl = 0.0;
		sampler.update(rng, &ppl);
		cerr << iter << ' ' << ppl << "\n";
	    }
	    else {
		sampler.update(rng);
		cerr << iter << "\n";
	    }
	}

	// Save the model
	namespace fs = boost::filesystem;
	fs::create_directory(model);
	if (!fs::exists(model)) bye("Cannot create the model directory " + model);

	fs::path basedir(model);

	{
	    cerr << "Write vocabulary\n";
	    fs::ofstream vocab_out(basedir / "vocab");
	    vocab.save(vocab_out);

	    cerr << "Write model\n";
	    fs::ofstream model_out(basedir / "model");
	    training_set.save(model_out);

	    cerr << "Write docno\n";
	    fs::ofstream docno_out(basedir / "docno");
	    foreach (const string& no, docno) { docno_out << no << '\n'; }
	}

	if (output == "phi_theta") {
	    cerr << "Output phi_theta\n";
	    sampler.outputPhiTheta(cout, docno);
	}
	else if (output == "theta") {
	    cerr << "Output theta\n";
	    sampler.outputTheta(cout, docno);
	}
    }
    else if (g["test"]) {
	// 1) Scan input
	namespace fs = boost::filesystem;
	if (!fs::exists(model)) bye("Cannot open the model directory " + model);
	fs::path basedir(model);

	fs::ifstream vocab_in(basedir / "vocab");
	Vocabulary<> vocab(vocab_in);

	vector<unsigned int> wseq;
	vector<unsigned int> doc;
	vector<unsigned int> sent;
	vector<string> docno;

	unsigned int docno_id = 0;
	unsigned int sentno_id = 0;
	string last_doc;

	foreach (const string& filename, arg) {
	    AutoIn in(filename);
	    string words, word;

	    while (getline(in(), words)) {
		istringstream wstream(words);

		wstream >> word;
		docno.push_back(word);

		strip_after_first(word, ':');
		if (last_doc != word) {
		    if (!last_doc.empty()) ++docno_id;
		    last_doc = word;
		}

		while (wstream >> word) {
		    wseq.push_back(vocab.encode(word));
		    doc.push_back(docno_id);
		    sent.push_back(sentno_id);
		}
		++sentno_id;
	    }
	}
	++docno_id;

	fs::ifstream model_in(basedir / "model");
	SLDAModel<> training_set(model_in);
	SLDAModel<> test_set(wseq.size(), docno_id, sentno_id, training_set.getS(), training_set.getT(), vocab.size(), wseq, doc, sent, rng);
	SLDAGibbsSampler<SLDAModel<> > sampler(test_set, training_set, alpha, beta, delta);

	// NOTE: Clean up memory space explicitly
	wseq.clear(); 
	doc.clear(); 
	sent.clear();

	// Run the sampler
	for (unsigned int iter = 1; iter <= niter; ++iter) {
	    if (interval && iter % interval == 0) {
		long double ppl = 0.0;
		sampler.update(rng, &ppl);
		cerr << iter << ' ' << ppl << "\n";
	    }
	    else {
		sampler.update(rng);
		cerr << iter << "\n";
	    }
	}

	if (output == "phi_theta") {
	    cerr << "Output phi_theta\n";
	    sampler.outputPhiTheta(cout, docno);
	}
	else if (output == "theta") {
	    cerr << "Output theta\n";
	    sampler.outputTheta(cout, docno);
	}
    }
    else if (g["test-cohesion"]) {
	// 1) Scan input
	namespace fs = boost::filesystem;
	if (!fs::exists(model)) bye("Cannot open the model directory " + model);
	fs::path basedir(model);

	fs::ifstream model_in(basedir / "model");
	SLDAModel<> training_set(model_in);

	// 2) Read the whole test collection in
	vector<string> sentences;
	vector<string> sentno;
	vector<string> docno;

	foreach (const string& filename, arg) {
	    AutoIn in(filename);
	    string line;

	    while (getline(in(), line)) {
		istringstream iss(line);
		string word;

		iss >> word;
		sentno.push_back(word); // Save the sentence no.
		strip_after_first(word, ':');
		docno.push_back(word); // Save the docno

		strip_before_first(line, ' ');
		sentences.push_back(line);
	    }
	}

	unsigned int size = docno.size();
	for (unsigned int i = 0; i < size; ++i) {
	    string current_docno = docno[i];
	    string current_sentno = sentno[i];

	    // Stuff that I seemed to write a thousand time...
	    typedef pair<string, float> item;
	    vector<item> rank;

	    // Locate the current document in [head, tail)
	    unsigned head = i, tail = i;
	    while (head > 0 && docno[head - 1] == current_docno) --head;
	    while (tail < size && docno[tail] == current_docno) ++tail;

	    for (unsigned int j = 0; j < size; ++j) {
		if (j != i && j >= head && j < tail) continue;

		// Ready
		fs::ifstream vocab_in(basedir / "vocab");
		Vocabulary<> vocab(vocab_in);

		vector<unsigned int> wseq;
		vector<unsigned int> doc;
		vector<unsigned int> sent;
		unsigned int sentid = 0;

		// Steady...
		for (unsigned int k = head; k < tail; ++k) {
		    string& sentence = (k == i? sentences[j]: sentences[k]);
		    istringstream iss(sentence);
		    string word;

		    while (iss >> word) {
			wseq.push_back(vocab.encode(word));
			doc.push_back(0);
			sent.push_back(sentid);
		    }
		    ++sentid;
		}

		// Go!
		SLDAModel<> test_set(wseq.size(), 1, sentid, training_set.getS(), training_set.getT(), vocab.size(), wseq, doc, sent, rng);
		SLDAGibbsSampler<SLDAModel<> > sampler(test_set, training_set, alpha, beta, delta);

		wseq.clear();
		doc.clear();
		sent.clear();

		// Run the sampler
		long double ppl = 0.0;
		for (unsigned int iter = 1; iter <= niter; ++iter) sampler.update(rng);

		// Take 20 consecutive samples and make average
		long double ppl_sum = 0.0;
		for (unsigned int n = 0; n < 20; ++n) {
		    sampler.update(rng, &ppl);
		    ppl_sum += ppl;
		}

		ppl = ppl_sum / 20;

		cerr << sentno[i] << ' ' << sentno[j] << ' ' << ppl << '\n';
		rank.push_back(item(sentno[j], ppl));
	    }

	    if (rank.size() > top_n) 
		nth_element(rank.begin(), rank.begin() + top_n, rank.end(), second_cmp());

	    rank.erase(rank.begin() + top_n, rank.end());
	    stable_sort(rank.begin(), rank.end(), second_cmp());

	    foreach (const item& r, rank) {
		cout << current_sentno << ' ' << r.first << ' ' << r.second << '\n';
	    }
	}
    }
    else if (g["dump"]) {
	// 1) Scan input
	namespace fs = boost::filesystem;
	if (!fs::exists(model)) bye("Cannot open the model directory " + model);
	fs::path basedir(model);

	fs::ifstream model_in(basedir / "model");
	SLDAModel<> training_set(model_in);
	SLDAGibbsSampler<SLDAModel<> > sampler(training_set, alpha, beta, delta);

	vector<string> vocab;
	fs::ifstream vocab_in(basedir / "vocab");
	copy(istream_iterator<string>(vocab_in),
		istream_iterator<string>(), back_inserter(vocab));

	if (output == "word-cluster") {
	    sampler.outputWordCluster(cout, vocab, top_n);
	}
	//--------------------------------------------------
	// else if (output == "topic-cluster") {
	//     sampler.outputTopicCluster(cout);
	// }
	// else if (output == "sentence-topic-cluster") {
	//     sampler.outputSentenceTopicCluster(cout);
	// }
	//-------------------------------------------------- 
	else {
	    cerr << "No such format";
	}
    }

    cerr << "Done\n";
    return 0;
}

