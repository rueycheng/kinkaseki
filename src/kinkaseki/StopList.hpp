#ifndef KINKASEKI_STOPLIST_H
#define KINKASEKI_STOPLIST_H

#include <iosfwd>
#include <unordered_set>

namespace kinkaseki {

template<typename Tp, 
    typename Container = std::unordered_set<Tp> > // FIXME
class StopList_Impl {
protected:
    Container map;

public:
    StopList_Impl() {}

    StopList_Impl(const typename Tp::value_type* array[], int n) {
	insert(&array[0], &array[n]);
    }

    template<typename Iterator> 
    StopList_Impl(Iterator first, Iterator last) {
	insert(first, last);
    }

    template<typename Iterator>
    void insert(Iterator first, Iterator last) {
	while (first != last) map.insert(*first++);
    }

    bool match(const Tp& v) {
	typename Container::iterator iter = map.find(v);
	return iter != map.end();
    }
};

typedef StopList_Impl<std::string> StopList; 

}
#endif
