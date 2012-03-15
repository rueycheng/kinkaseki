#ifndef KINKASEKI_EXTERNAL_SORTER_H
#define KINKASEKI_EXTERNAL_SORTER_H

#include <algorithm>
#include <vector>

#include "kinkaseki/ExternalMerger.hpp"

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

	std::stable_sort(container.begin(), container.end());
	merger.addRun(container.rbegin(), container.rend());  // in reverse order
	container.clear();
    }

public:
    typedef typename Merger::Result Result;

    ExternalSorter(int memoryLimit) {
	maxSize = memoryLimit / sizeof(T);
    }

    void add(const T& element) {
	container.push_back(element);
	if (container.size() >= maxSize) sortAndFlush();
    }

    Result run() {
	sortAndFlush();
	return merger.merge();
    }
};

#endif
