#ifndef EXTERNAL_MERGER_H
#define EXTERNAL_MERGER_H

#include <algorithm>
#include <fstream>
#include <iostream>
#include <iterator>
#include <string>
#include <vector>
#include <boost/lexical_cast.hpp>

namespace util {

// ------------------------------------------------------------
/// @brief ExternalMerger
///
/// @tparam T
/// @tparam Less
// ------------------------------------------------------------
template<typename T, typename Less = std::less<T>>
class ExternalMerger {
protected:
    std::string prefix;
    std::vector<std::string> runFiles;

public:
    class Result {
	std::string file;
	std::ifstream* in; 
	std::istream_iterator<T> iter;
	bool autoDelete;

    public:
	Result(const std::string& file, bool autoDelete = false): 
	    file(file), in(new std::ifstream(file)), autoDelete(autoDelete) 
	{ }

	~Result() {
	    delete in;
	    if (autoDelete) ::unlink(file.c_str());
	}

	std::istream_iterator<T> iterator() {
	    return std::istream_iterator<T>(*in);
	}

	std::istream_iterator<T> end() {
	    return std::istream_iterator<T>();
	}
    };

    struct Comparator {
	Less comp;

	bool operator()(const std::istream_iterator<T>& lhs, const std::istream_iterator<T>& rhs) {
	    return comp(*lhs, *rhs);
	}
    };

    ExternalMerger(const std::string& prefix):
	prefix(prefix) {}

    ~ExternalMerger() {
	for (int i = 0; i < runFiles.size(); ++i) ::unlink(runFiles[i].c_str());
    }

    template<typename Iterator> void addRun(Iterator first, Iterator last) {
	std::string file = "/tmp/ExternalMerger-" 
	    + prefix + "-" + boost::lexical_cast<std::string>(runFiles.size());
	std::ofstream out(file);
	std::copy(first, last, std::ostream_iterator<T>(out, ""));

	runFiles.push_back(file);
    }

    void mergeToFile(const std::string& file) {
	if (runFiles.size() == 1) {
	    std::ifstream in(runFiles[0]);
	    std::ofstream out(file);

	    out << in.rdbuf();
	}
	else {
	    std::ofstream out(file);
	    std::ostream_iterator<T> result(out, "");

	    std::vector<std::ifstream*> pool;
	    std::vector<std::istream_iterator<T>> runs;

	    for (int i = 0; i < runFiles.size(); ++i) {
		pool.push_back(new std::ifstream(runFiles[i]));
		runs.push_back(std::istream_iterator<T>(*pool[i]));
	    }

	    Comparator comp;
	    std::istream_iterator<T> endOfRun;
	    std::make_heap(runs.begin(), runs.end(), comp);
	    while (!runs.empty()) {
		std::pop_heap(runs.begin(), runs.end(), comp);
		std::istream_iterator<T>& current = runs.back();

		*result++ = *current++;
		if (current != endOfRun) 
		    std::push_heap(runs.begin(), runs.end(), comp);
		else
		    runs.pop_back();
	    }

	    for (int i = 0; i < runFiles.size(); ++i) delete pool[i];
	}
    }

    Result merge() {
	std::string file = "/tmp/ExternalMerger-" + prefix + "-result";
	mergeToFile(file);
	return Result(file);
    }
};

}

#endif
