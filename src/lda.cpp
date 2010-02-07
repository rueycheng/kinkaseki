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
// LDAGibbsSampler
//-------------------------------------------------- 
template<class __LDAModel, typename __float = float, typename __count = unsigned int>
class LDAGibbsSampler {
public:
    typedef typename __LDAModel::size_type size_type;
    typedef typename __LDAModel::doc_type doc_type;
    typedef typename __LDAModel::topic_type topic_type;
    typedef typename __LDAModel::word_type word_type;
    typedef __LDAModel lda_type;
    typedef __float prob_type;
    typedef __count count_type;

private:
    struct Count {
	count_type* n_dz;
	count_type* n_d_;
	count_type* n_zw;
	count_type* n_z_;

	Count(doc_type m, topic_type t, word_type w) {
	    n_dz = new count_type[m * t];
	    n_d_ = new count_type[m];
	    n_zw = new count_type[t * w];
	    n_z_ = new count_type[t];
	    fill_n(n_dz, m * t, 0);
	    fill_n(n_d_, m, 0);
	    fill_n(n_zw, t * w, 0);
	    fill_n(n_z_, t, 0);
	}

	~Count() {
	    delete[] n_dz; 
	    delete[] n_d_; 
	    delete[] n_zw; 
	    delete[] n_z_; 
	}
    };

    lda_type& m;
    prob_type alpha, beta;

    LDAGibbsSampler(LDAGibbsSampler&) {} // No copy

    // Counts for the current model
    count_type* n_dz; // document-topic, M by T
    count_type* n_d_; // document-*, M
    count_type* n_zw; // topic-word, T by W
    count_type* n_z_; // topic-*, T

    // Counts for the reference model (if any)
    word_type Wref;
    count_type* nn_zw;
    count_type* nn_z_;

    prob_type* sample; // cdf for p(t|t*,w,w*), T + 1

    // Model initialization
    void alloc(doc_type m, topic_type t, word_type w, word_type wref = 0) {
	n_dz = new count_type[m * t];
	n_d_ = new count_type[m];
	n_zw = new count_type[t * w];
	n_z_ = new count_type[t];
	fill_n(n_dz, m * t, 0);
	fill_n(n_d_, m, 0);
	fill_n(n_zw, t * w, 0);
	fill_n(n_z_, t, 0);

	if (wref) {
	    Wref = wref;
	    nn_zw = new count_type[t * wref];
	    nn_z_ = new count_type[t];
	    fill_n(nn_zw, t * wref, 0);
	    fill_n(nn_z_, t, 0);
	} else {
	    Wref = 0;
	    nn_zw = 0; // No ref
	    nn_z_ = 0; // No ref
	}

	sample = new prob_type[t + 1];
    }

    void free() {
	delete[] n_dz; 
	delete[] n_d_; 
	delete[] n_zw; 
	delete[] n_z_; 
	delete[] nn_zw; 
	delete[] nn_z_; 
	delete[] sample;
    }

    void collectCount(lda_type& m) {
	for (size_type i = 0; i < m.N; ++i) {
	    const doc_type& d = m.d[i];
	    const topic_type& z = m.z[i];
	    const word_type& w = m.w[i];
	    ++n_dz[d * m.T + z];
	    ++n_d_[d];
	    ++n_zw[z * m.W + w];
	    ++n_z_[z];
	}
    }

    void collectRefCount(lda_type& refm) {
	for (size_type i = 0; i < refm.N; ++i) {
	    const topic_type& z = refm.z[i];
	    const word_type& w = refm.w[i];
	    ++nn_zw[z * refm.W + w];
	    ++nn_z_[z];
	}
    }

public:
    // One-model initialization
    LDAGibbsSampler(lda_type& m, prob_type alpha, prob_type beta): 
	m(m), alpha(alpha), beta(beta) 
    {
	this->alpha = (alpha == 0.0? 50.0/m.T: alpha);
	alloc(m.M, m.T, m.W);
	collectCount(m);
    }

    // With reference model (`test')
    LDAGibbsSampler(lda_type& m, lda_type& refm, prob_type alpha, prob_type beta):
	m(m), alpha(alpha), beta(beta)
    {
	this->alpha = (alpha == 0.0? 50.0/m.T: alpha); 
	alloc(m.M, m.T, m.W, refm.W);
	collectCount(m);
	collectRefCount(refm);
    }

    void update(RNG<>& rng, prob_type* ppl = 0) {
	const prob_type W_beta = m.W * beta;
	const prob_type T_alpha = m.T * alpha;

	if (ppl) *ppl = 0.0;

	for (size_type i = 0; i < m.N; ++i) {
	    const doc_type& d = m.d[i];
	    topic_type& z = m.z[i]; // We might change this value later (non-const)
	    const word_type& w = m.w[i];

	    // 1) Decrement.  Like tuple (d,z,w) never exists.
	    const size_type dT = d * m.T;

	    --n_dz[dT + z];
	    --n_d_[d];
	    --n_zw[z * m.W + w];
	    --n_z_[z];

	    // 2) Collect conditional probabilities
	    prob_type prob_of_sum = 0.0;
	    prob_type f2d = n_d_[d] + T_alpha;
	    for (topic_type zz = 0; zz < m.T; ++zz) {
		prob_type f2n = n_dz[dT + zz] + alpha;
		prob_type f1n = n_zw[zz * m.W + w] + beta; //+ ((nn_zw && w < Wref)? nn_zw[zz * Wref + w]: 0);
		prob_type f1d = n_z_[zz] + W_beta; // + (nn_z_? nn_z_[zz]: 0); 
		sample[zz] = (prob_of_sum += (f1n / f1d) * (f2n / f2d)); // Also obtain CDF
	    }

	    // NOTE: Sometimes old-fashion way works best
	    if (ppl) *ppl += log(prob_of_sum); 

	    // 3) Toss the dice and get the sample
	    prob_type toss = rng(sample[m.T - 1]);

	    topic_type zz = 0;
	    while (toss > sample[zz]) ++zz;
	    m.z[i] = zz; // Gotcha!

	    // 4) Increment
	    ++n_dz[dT + zz];
	    ++n_d_[d];
	    ++n_zw[zz * m.W + w];
	    ++n_z_[zz];
	}

	if (ppl) *ppl = exp(-(*ppl) / m.N);
    }

    void outputTheta(ostream& o, vector<string>& docno) {
	const prob_type T_alpha = m.T * alpha;
	const doc_type M = m.M;
	const topic_type T = m.T;
	bool use_docno = docno.size() == M;

	for (doc_type d = 0; d < M; ++d) {
	    if (use_docno) o << docno.at(d) << ' ';

	    prob_type denom = T_alpha + n_d_[d];
	    for (topic_type z = 0; z < T; ++z) o << (alpha + n_dz[d * T + z]) / denom << ' ';
	    o << "\n";
	}
    }

    void outputPhiTheta(ostream& o, vector<string>& docno) {
	const prob_type W_beta = m.W * beta;
	const prob_type T_alpha = m.T * alpha;
	const doc_type M = m.M;
	const word_type W = m.W;
	const topic_type T = m.T;
	bool use_docno = docno.size() == M;

	for (doc_type d = 0; d < M; ++d) {
	    if (use_docno) o << docno.at(d) << ' ';

	    prob_type theta_denom = T_alpha + n_d_[d];
	    for (word_type w = 0; w < W; ++w) {
		prob_type sum = 0;
		for (topic_type z = 0; z < T; ++z) {
		    prob_type theta_dz = (alpha + n_dz[d * T + z]) / theta_denom;
		    prob_type phi_zw = (beta + n_zw[z * W + w]) / (W_beta + n_z_[z]);
		    sum += theta_dz * phi_zw;
		}
		o << sum << ' ';
	    }
	    o << "\n";
	}
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

    ~LDAGibbsSampler() { free(); }
};

//--------------------------------------------------
// LDAModel
//-------------------------------------------------- 
template<typename __size = unsigned int,
    typename __doc = unsigned int, 
    typename __topic = unsigned short,
    typename __word = unsigned int>
class LDAModel {
public:
    typedef __size size_type;
    typedef __doc doc_type;
    typedef __topic topic_type;
    typedef __word word_type;

private:
    size_type N;
    doc_type M; 
    topic_type T;
    word_type W; 
    doc_type *d; // sized N, each in the range of [0, M)
    topic_type *z; // sized N, each in the range of [0, T)
    word_type *w; // sized N, each in the range of [0, W)

    LDAModel(LDAModel& model) {} // No copy
    void alloc(size_type n) { d = new doc_type[n]; z = new topic_type[n]; w = new word_type[n]; }
    void free() { delete[] d; delete[] z; delete[] w; }

public:
    LDAModel(size_type N, doc_type M, topic_type T, word_type W, vector<word_type>& wseq, 
	    vector<doc_type>& doc, RNG<>& rng): N(N), M(M), T(T), W(W) 
    { 
	alloc(N);

	typename vector<word_type>::const_iterator witer = wseq.begin();
	typename vector<doc_type>::const_iterator diter = doc.begin();
	for (size_type i = 0; i < N; ++i) {
	    d[i] = *diter++;
	    z[i] = rng(T);
	    w[i] = *witer++;
	}
    }

    explicit LDAModel(istream& seq) { 
	seq >> N >> M >> T >> W;
	alloc(N);

	for (size_type i = 0; i < N; ++i) seq >> d[i] >> z[i] >> w[i];
    }

    void save(ostream& seq) {
	seq << N << ' ' << M << ' ' << T << ' ' << W << "\n";
	for (size_type i = 0; i < N; ++i) seq << d[i] << ' ' << z[i] << ' ' << w[i] << "\n";
    }

    ~LDAModel() { free(); }

    size_type getN() { return N; }
    doc_type getD() { return M; }
    topic_type getT() { return T; } 
    word_type getW() { return W; }

    friend class LDAGibbsSampler< LDAModel<size_type, doc_type, topic_type, word_type> >;
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
    unsigned int ntopic = 100;
    unsigned int niter = 1000;
    unsigned int interval = 1;
    unsigned int sample_num = 1;
    unsigned int sample_lag = 50;
    string model = "model.unnamed";
    string output = "none";
    vector<string> arg;

    g   << $("train", "Run in the training mode")
	<< $("test", "Run in the test mode")
	<< $("use-docno", "Treat the first token of each input line as the docno")
	<< $(&output, "output", "Specify the output: 'none', 'phi_theta'")
	<< $(&model, "model,m", "Path to the model")
	<< $(&alpha, "alpha,a", "Specify the parameter 'alpha'")
	<< $(&beta, "beta,b", "Specify the parameter 'beta'")
	<< $(&ntopic, "ntopic,T", "Specify the number of word topics")
	<< $(&niter, "niter,I", "Specify the number of runs")
	<< $(&interval, "interval", "Specify the interval (#iter) between perplexity reports")
	<< $(&sample_num, "sample-num", "Specify the number of read-outs (after burned-in)")
	<< $(&sample_lag, "sample-lag", "Specify the number of iteration between each read-out")
	<< $(&arg, "arg", "", -1)
	<< $$$("lda {--train|--test} [options..] [input-files..]");

    // 0) Prepare utilities
    RNG<> rng;

    if (!g["train"] && !g["test"]) bye(g);
    if (g["train"] && g["test"]) bye("Specify either --train or --test.  Not both.");

    if (arg.empty()) arg.push_back("-");

    if (g["train"]) {
	// 1) Scan input
	Vocabulary<> vocab;
	vector<unsigned int> wseq;
	vector<unsigned int> doc;
	vector<string> docno; // Only used when '--use-docno'
	unsigned int docno_id = 0;

	foreach (const string& filename, arg) {
	    AutoIn in(filename);
	    string words, word;

	    while (getline(in(), words)) {
		istringstream wstream(words);

		if (g["use-docno"]) {
		    wstream >> word;
		    docno.push_back(word);
		}

		while (wstream >> word) {
		    wseq.push_back(vocab.encode(word));
		    doc.push_back(docno_id);
		}
		++docno_id;
	    }
	}

	LDAModel<> training_set(wseq.size(), docno_id, ntopic, vocab.size(), wseq, doc, rng);
	LDAGibbsSampler<LDAModel<> > sampler(training_set, alpha, beta);

	// NOTE: Clean up memory space explicitly
	wseq.clear(); 
	doc.clear(); 

	// Run the sampler
	for (unsigned int iter = 1; iter <= niter; ++iter) {
	    if (interval && iter % interval == 0) {
		float ppl = 0.0;
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

	    if (g["use-docno"]) {
		cerr << "Write docno\n";
		fs::ofstream docno_out(basedir / "docno");
		foreach (const string& no, docno) { docno_out << no << '\n'; }
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
    else if (g["test"]) {
	// 1) Scan input
	namespace fs = boost::filesystem;
	if (!fs::exists(model)) bye("Cannot open the model directory " + model);
	fs::path basedir(model);

	fs::ifstream vocab_in(basedir / "vocab");
	Vocabulary<> vocab(vocab_in);

	vector<unsigned int> wseq;
	vector<unsigned int> doc;
	vector<string> docno; // Only used when '--use-docno'
	unsigned int docno_id = 0;

	foreach (const string& filename, arg) {
	    AutoIn in(filename);
	    string words, word;

	    while (getline(in(), words)) {
		istringstream wstream(words);

		if (g["use-docno"]) {
		    wstream >> word;
		    docno.push_back(word);
		}

		while (wstream >> word) {
		    wseq.push_back(vocab.encode(word));
		    doc.push_back(docno_id);
		}
		++docno_id;
	    }
	}

	fs::ifstream model_in(basedir / "model");
	LDAModel<> training_set(model_in);
	LDAModel<> test_set(wseq.size(), docno_id, training_set.getT(), vocab.size(), wseq, doc, rng);
	LDAGibbsSampler<LDAModel<> > sampler(test_set, training_set, alpha, beta);

	// NOTE: Clean up memory space explicitly
	if (true) {
	    wseq.clear(); 
	    doc.clear(); 
	}

	// Run the sampler
	for (unsigned int iter = 1; iter <= niter; ++iter) {
	    if (interval && iter % interval == 0) {
		float ppl = 0.0;
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

    cerr << "Done\n";
    return 0;
}

