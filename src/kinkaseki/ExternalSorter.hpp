#ifndef EXTERNAL_SORTER_H
#define EXTERNAL_SORTER_H

#include <algorithm>
#include <vector>

#include "ExternalMerger.hpp"

namespace util {

// ------------------------------------------------------------
/// @brief ExternalSorter
///
/// @tparam T
/// @tparam Less
// ------------------------------------------------------------
template<typename T, typename Less = std::less<T>> 
class ExternalSorter {
public:
    typedef ExternalMerger<T, Less> Merger;
    typedef std::vector<T> Container;

protected:
    Merger merger;
    Container container;
    int maxSize;

    void sortAndFlush() {
	if (container.empty()) return;

	std::stable_sort(container.begin(), container.end(), Less());
	// merger.addRun(container.rbegin(), container.rend());  // in reverse order
	merger.addRun(container.begin(), container.end());  // in reverse order
	container.clear();
    }

public:
    typedef typename Merger::Result Result;

    ExternalSorter(int memoryLimit, const std::string& prefix):
	maxSize(memoryLimit / sizeof(T)), merger(prefix) {}

    void add(const T& element) {
	container.push_back(element);
	if (container.size() >= maxSize) sortAndFlush();
    }

    Result run() {
	sortAndFlush();
	return merger.merge();
    }
};

}

#endif
