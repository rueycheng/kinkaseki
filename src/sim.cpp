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
// Similarity functions
//-------------------------------------------------- 
float kl(const vector<float>& q, const vector<float>& d) {
    float result = 0;

    vector<float>::const_iterator qq = q.begin();
    vector<float>::const_iterator qend = q.end();
    vector<float>::const_iterator dd = d.begin();
    vector<float>::const_iterator dend = d.end();

    while (qq != qend && dd != dend) {
	float p = *qq++;
	result += p * log(p / *dd++);
    }

    return result;
}

//--------------------------------------------------
// Main program
//-------------------------------------------------- 
int main(int argc, char** argv) {

    // Getopt
    vector<string> input;
    unsigned int top_n = 100;

    Getopt g(argc, argv);
    g   << $("kl", "Calculate KL-divergence")
	<< $("cosine", "Calculate cosine similarity")
	<< $("in-memory", "Cache all the representation in memory")
	<< $("invert-role", "Invert the role of q and d in KL-divergence")
	<< $(&top_n, "top-n", "Show only top N results")
	<< $(&input, "input", "", -1)
	<< $$$("files..");

    if (input.size() != 2) throw 0;

    namespace fs = boost::filesystem;

    AutoIn query_in(input.at(1));
    string query_line;
    while (getline(query_in(), query_line)) {
	string qno;
	vector<float> qv;

	istringstream iss(query_line);
	iss >> qno;
	copy(istream_iterator<float>(iss),
		istream_iterator<float>(), back_inserter(qv));
	
	typedef pair<string, float> item;
	vector<item> rank;

	AutoIn doc_in(input.at(0));
	string doc_line;
	while (getline(doc_in(), doc_line)) {
	    string dno;
	    vector<float> dv;

	    istringstream iss2(doc_line);
	    iss2 >> dno;
	    copy(istream_iterator<float>(iss2),
		    istream_iterator<float>(), back_inserter(dv));

	    float score = g["invert-role"]? kl(dv, qv): kl(qv, dv);
	    rank.push_back(item(dno, score));
	}

	if (rank.size() > top_n) {
	    nth_element(rank.begin(), rank.begin() + top_n, rank.end(), second_cmp());
	    rank.erase(rank.begin() + top_n, rank.end());
	}

	stable_sort(rank.begin(), rank.end(), second_cmp());
	foreach (const item& r, rank) {
	    cout << qno << ' ' << r.first << ' ' << r.second << '\n';
	}
    }

    return 0;
}
