#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <algorithm>
#include <iterator>
#include <ctime>
#include <limits>
#include <boost/format.hpp>
#include <boost/random.hpp>
#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>

#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH

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

#include "magicbox.h"

using namespace std;
using namespace magicbox;

class MissingModelPath {};
class OutOfMemory {};
class Crap {};

//--------------------------------------------------
// RandomSequence
//-------------------------------------------------- 
template<typename __R = double, class __Generator = boost::mt19937, 
    class __Distribution = boost::uniform_01<__Generator, __R> >
class RandomSequence {
    __Distribution d;

public:
    RandomSequence(uint32_t seed = std::time(0)): d(__Generator(seed)) {}
    __R operator()() { return d(); }
    __R operator()(__R f) { return f * d(); }
    unsigned int operator()(unsigned int n) { return (int)(n * d()); }
};

//--------------------------------------------------
// SLDAGibbsSampler
//-------------------------------------------------- 
template<typename __R = float, typename __I = unsigned int>
class SLDAGibbsSampler {
    __R alpha, beta, delta;
    __I N, L, M, W;
    short S, T;
    const __I *ws, *ww; 
    short *wt;
    const __I *sd;
    short *st; // The values do not stay constant

    __I *n_wt, *n_xt, *n_st, *n_sx, *n_ds, *n_dx;
    RandomSequence<> seq;
    RandomSequence<long double> seq2;
    double *distro; 
    long double *distro2;

public:
    SLDAGibbsSampler(__R alpha, __R beta, __R delta,
	    __I N, __I M, __I L, __I W, short S, short T, 
	    const __I ws[], const __I ww[], short wt[], const __I sd[], short st[]): 
	alpha(alpha), beta(beta), delta(delta), N(N), L(L), M(M), W(W), S(S), T(T), 
	ws(ws), ww(ww), wt(wt), sd(sd), st(st)
    {
	n_wt = new __I[W * T];
	n_xt = new __I[T];
	n_st = new __I[S * T];
	n_sx = new __I[S];
	n_ds = new __I[M * S];
	n_dx = new __I[M];

	distro = new double[T];
	distro2 = new long double[S];

	if (!n_wt || !n_xt || !n_st || !n_sx || !n_ds || !n_dx || !distro) throw OutOfMemory();

	fill_n(n_wt, W * T, 0);
	fill_n(n_xt, T, 0);
	fill_n(n_st, S * T, 0);
	fill_n(n_sx, S, 0);
	fill_n(n_ds, M * S, 0);
	fill_n(n_dx, M, 0);

	for (__I i = 0; i < N; ++i) {
	    ++n_wt[ww[i] * T + wt[i]];
	    ++n_xt[wt[i]];
	    ++n_st[st[ws[i]] * T + wt[i]];
	    ++n_sx[st[ws[i]]];
	    ++n_ds[sd[ws[i]] * S + st[ws[i]]];
	    ++n_dx[sd[ws[i]]];
	}
    }

    ~SLDAGibbsSampler() {
	delete[] n_wt;
	delete[] n_xt;
	delete[] n_st;
	delete[] n_sx;
	delete[] n_ds;
	delete[] n_dx;
	delete[] distro;
	delete[] distro2;
    }

    void update() {
	__R W_beta = W * beta;
	__R T_delta = T * delta;
	__R S_alpha = S * alpha;

	// Pass 1: Iterate over words
	for (__I i = 0; i < N; ++i) {
	    __I w = ww[i];
	    short k = wt[i], l = st[ws[i]];

	    // 0) Stuff
	    __I wT = w * T;
	    __I lT = l * T;

	    // 1) Decrement.  Imagine that (wd_i, wt_i, ww_i) never exists...
	    --n_wt[wT + k];
	    --n_xt[k];
	    --n_st[lT + k];
	    --n_sx[l];

	    // 2) Collect sampling probabilities
	    __R f2d = T_delta + n_sx[l];
	    for (short kk = 0; kk < T; ++kk) {
		__R f1n = beta + n_wt[wT + kk];
	        __R f1d = W_beta + n_xt[kk];
		__R f2n = delta + n_st[lT + kk];
		distro[kk] = (f1n / f1d) * (f2n / f2d);
	    }

	    for (short kk = 1; kk < T; ++kk) distro[kk] += distro[kk - 1]; // Obtain CDF
	    double toss = seq(distro[T - 1]); // Toss the dice

	    short kk = 0;
	    while (toss > distro[kk]) ++kk; // Sample!

	    // 3) Increment.
	    wt[i] = kk;

	    ++n_wt[wT + kk];
	    ++n_xt[kk];
	    ++n_st[lT + kk];
	    ++n_sx[l];
	}

	// Pass 2: Iterate over sentences
	__I i = 0;
	while (i < N) {
	    __I s = ws[i];
	    __I d = sd[s];
	    short l = st[s];

	    // 0) Locate all the words in the same sentence
	    __I end = i + 1;
	    while (end < N && ws[end] == s) ++end;

	    __I dS = d * S;
	    __I lT = l * T;

	    // 1) Decrement all the counts made between word [i, end).
	    for (__I j = i; j < end; ++j) --n_st[lT + wt[j]];
	    n_sx[l] -= (end - i);
	    n_ds[dS + l] -= (end - i);
	    n_dx[d] -= (end - i);

	    // 2) Calcualte CDF and run the sampler
	    __R f2d = S_alpha + n_dx[d];
	    for (short ll = 0; ll < S; ++ll) {
		__I llT = ll * T;
		for (short k = 0; k < T; ++k) distro[k] = n_st[llT + k];
		__I total = n_sx[ll];

		double logf1 = 0.0; 
		for (__I j = i; j < end; ++j) {
		    __I count = distro[wt[j]];
		    if (count == 0) count = 1;
		    if (total == 0) total = 1;
		    logf1 += log(count) - log(total);

		    ++distro[wt[j]];
		    ++total;
		}

		__R f2n = alpha + n_ds[dS + ll];
		distro2[ll] = logf1 + log(f2n) - log(f2d);
	    }

	    __R dmax = numeric_limits<__R>::min(), dmin = numeric_limits<__R>::max();
	    for (short ll = 0; ll < S; ++ll) {
		if (distro2[ll] > dmax) dmax = distro2[ll];
		if (distro2[ll] < dmin) dmin = distro2[ll];
	    }

	    if (dmax - dmin > 11356 * 2) throw Crap();

	    long double sum = 0;
	    for (short ll = 0; ll < S; ++ll) {
		distro2[ll] = exp(distro2[ll] - dmin);
		sum += distro2[ll];
	    }

	    for (short ll = 0; ll < S; ++ll) distro2[ll] /= sum;
	    for (short ll = 1; ll < S; ++ll) distro2[ll] += distro2[ll - 1]; // Obtain CDF
	    double toss = seq2(distro2[S - 1]); // Toss the dice
	    short ll = 0;
	    while (toss > distro2[ll]) ++ll; // Sample!

	    // 3) Increment all the counts back
	    st[s] = ll;

	    __I llT = ll * T;
	    for (__I j = i; j < end; ++j) ++n_st[llT + wt[j]];
	    n_sx[ll] += (end - i);
	    n_ds[dS + ll] += (end - i);
	    n_dx[d] += (end - i);

	    // 4) Proceed to the next unprocessed word
	    i = end;
	}
    }

    __R perplexity() {
	__R ppl = 0.0;

	for (__I i = 0; i < N; ++i) {
	    __I w = ww[i];
	    __I d = sd[ws[i]];

	    // 0) Stuff
	    __I wT = w * T;
	    __I dS = d * S;

	    __R S_alpha = S * alpha;
	    __R T_delta = T * delta;
	    __R W_beta = W * beta;

	    __R sum = 0.0;
	    for (short ll = 0; ll < S; ++ll) {
		__R f_theta = (alpha + n_ds[dS + ll]) / (S_alpha + n_dx[d]);
		__I llT = ll * T;

		__R partial_sum = 0.0;
		for (short kk = 0; kk < T; ++kk) {
		    __R f_tau = (delta + n_st[llT + kk]) / (T_delta + n_sx[ll]);
		    __R f_phi = (beta + n_wt[wT + kk])/ (W_beta + n_xt[kk]);
		    partial_sum += f_tau * f_phi;
		}

		sum += f_theta * partial_sum;
	    }

	    ppl += log(sum);
	}

	return exp(-ppl / N);
    }

    void outputStream(ostream& o) {
	o << N << endl; // The first line: N
	for (__I i = 0; i < N; ++i) 
	    o << sd[ws[i]] << ' ' << st[ws[i]] << ' ' << ws[i] << ' ' << ww[i] << ' ' << wt[i] << endl;
    }

    void outputTheta(ostream& o) {
	__R S_alpha = S * alpha;
	o << alpha << ' ' << M << ' ' << S << endl; // The first line: alpha, M, and S
	for (__I j = 0; j < M; ++j) {
	    __R denom = S_alpha + n_dx[j];
	    for (short l = 0; l < S; ++l) o << (alpha + n_ds[j * S + l]) / denom << ' ';
	    o << endl;
	}
    }

    void outputTau(ostream& o) {
	__R T_delta = T * delta;
	o << delta << ' ' << S << ' ' << T << endl; // The first line: delta, S, and T
	for (short l = 0; l < S; ++l) {
	    __R denom = T_delta + n_sx[l];
	    for (short k = 0; k < T; ++k) o << (delta + n_st[l * T + k]) / denom << ' ';
	    o << endl;
	}
    }
    
    void outputPhi(ostream& o) {
	__R W_beta = W * beta;
	o << beta << ' ' << T << ' ' << W << endl; // The first line: beta, T, and W
	for (short k = 0; k < T; ++k) {
	    __R denom = W_beta + n_xt[k];
	    for (__I w = 0; w < W; ++w) o << (beta + n_wt[w * T + k]) / denom << ' ';
	    o << endl;
	}
    }
};

//--------------------------------------------------
// Vocabulary
//-------------------------------------------------- 
class Vocabulary {
    unordered_map<string, unsigned int> vocab;
    vector<string> rvocab;

public:
    Vocabulary() { rvocab.push_back("<unk>"); }
    unsigned int size() { return rvocab.size(); }
    string decode(unsigned int code) { return rvocab.at(code); }
    unsigned int encode(string& word) {
	if (vocab.find(word) != vocab.end()) return vocab[word];
	rvocab.push_back(word);
	return vocab[word] = rvocab.size() - 1;
    }

    void output(ostream& o) { foreach (const string& word, rvocab) { o << word << endl; } }
};

//--------------------------------------------------
// Main program
//-------------------------------------------------- 
int main(int argc, char** argv) {
    // Getopt stuff
    Getopt g(argc, argv);

    float alpha = 0;
    float beta = 0.01;
    float delta = 0;
    unsigned int nstopic = 100;
    unsigned int ntopic = 100;
    unsigned int nrun = 1000;
    unsigned int ntop = 100;
    string model = "";
    vector<string> input;

    g   << $("verbose", "Show verbose output")
	<< $(&model, "model,m", "Path to the output model")
	<< $(&alpha, "alpha,a", "Specify the parameter 'alpha'")
	<< $(&beta, "beta,b", "Specify the parameter 'beta'")
	<< $(&delta, "delta,d", "Specify the parameter 'delta'")
	<< $(&nstopic, "nstopic,S", "Specify the number of sentence topics")
	<< $(&ntopic, "ntopic,T", "Specify the number of word topics")
	<< $(&nrun, "nrun,R", "Specify the number of runs")
	<< $(&ntop, "show-top", "Show only top S results")
	<< $(&input, "input", "", -1)
	<< $$$("lda [options..] [input-file..]");

    // Validate
    if (model.empty()) throw MissingModelPath();
    if (input.empty()) input.push_back("-");

    // 1) Scan through all the texts.  Translate words into numbers.
    typedef pair<unsigned int, unsigned int> T2;
    vector<T2> word2sent;
    vector<T2> sent2doc;

    unsigned int docid = 0;
    unsigned int sentid = 0;
    Vocabulary vocab; // Word ID starts from 1 (0 stands for OOV)

    foreach (const string& input_path, input) {
	ifstream inf;
	if (input_path != "-") inf.open(input_path.c_str());
	istream& in = (input_path == "-")? cin: inf;

	string line, word;
	while (getline(in, line)) {
	    istringstream iss(line);
	    bool has_word = false;

	    while (iss >> word) {
		has_word = true;
		unsigned id = vocab.encode(word);
		word2sent.push_back(make_pair(sentid, id));
	    }

	    if (has_word) sent2doc.push_back(make_pair(docid, sentid++));
	    else ++docid; 
	}
    }

    // 2) Get parameters
    unsigned int nword = vocab.size();
    unsigned int nsent = sentid;
    unsigned int ndoc = docid;
    unsigned int size = word2sent.size();

    cerr << "nword: " << nword << endl;
    cerr << "nsent: " << nsent << endl;
    cerr << "ndoc: " << ndoc << endl;
    cerr << "size: " << size << endl;

    if (alpha == 0) alpha = 50.0 / nstopic;
    if (delta == 0) delta = 50.0 / ntopic;

    // 3) Grab new arrays.  Fill in everything.
    unsigned int* w_sent = new unsigned int[size];
    unsigned int* w_word = new unsigned int[size];
    short* w_topic = new short[size];
    unsigned int* s_doc = new unsigned int[nsent];
    short* s_topic = new short[nsent]; 

    if (!w_sent || !w_word || !w_topic ||!s_doc ||!s_topic) throw OutOfMemory();

    RandomSequence<> seq;

    {
	unsigned int i = 0;
	foreach (const T2& p, word2sent) {
	    w_sent[i] = p.first;
	    w_word[i] = p.second;
	    w_topic[i] = seq(ntopic); // Start with arbitrary topic assignments
	    ++i;
	}

	word2sent.clear(); // Release memory

	foreach (const T2& p, sent2doc) {
	    s_doc[p.second] = p.first;
	    s_topic[p.second] = seq(nstopic); // Random start
	}

	sent2doc.clear();
    }

    // 4) Get the gibbs sampler up and running
    SLDAGibbsSampler<> gibbs(alpha, beta, delta,  // Symmetric priors
	    size, ndoc, nsent, nword, nstopic, ntopic, 
	    w_sent, w_word, w_topic, s_doc, s_topic);

    for (unsigned int iter = 0; iter < nrun; ++iter) {
	float ppl = gibbs.perplexity();
	cerr << "Iteration #" << iter + 1 << " (perplexity: " << ppl << ")" << endl;
	gibbs.update();
    }

    // 5) Save the model
    namespace fs = boost::filesystem;

    fs::create_directory(model);
    if (!fs::exists(model)) throw Crap();

    fs::path basedir(model);

    { 
	cerr << "Write to " << basedir / "vocab" << endl;
	fs::ofstream vocab_f(basedir / "vocab"); 
	vocab.output(vocab_f); 
    }

    { 
	cerr << "Write to " << basedir / "stream" << endl;
	fs::ofstream stream_f(basedir / "stream"); 
	gibbs.outputStream(stream_f); 
    }

    { 
	cerr << "Write to " << basedir / "theta" << endl;
	fs::ofstream theta_f(basedir / "theta"); 
	gibbs.outputTheta(theta_f); 
    }

    { 
	cerr << "Write to " << basedir / "tau" << endl;
	fs::ofstream tau_f(basedir / "tau"); 
	gibbs.outputTau(tau_f); 
    }

    { 
	cerr << "Write to " << basedir / "phi" << endl;
	fs::ofstream phi_f(basedir / "phi"); 
	gibbs.outputPhi(phi_f); 
    }

    // 6) House-cleansing
    delete[] w_sent;
    delete[] w_word;
    delete[] w_topic;
    delete[] s_doc;
    delete[] s_topic;

    return 0;
}

