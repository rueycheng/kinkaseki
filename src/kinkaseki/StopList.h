#ifndef KINKASEKI_STOPLIST_H
#define KINKASEKI_STOPLIST_H

#include <iosfwd>
#include <set>

namespace kinkaseki {

class StopList {
public:
    typedef std::string value_type;
    typedef std::set<std::string> container_type;

protected:
    container_type map;

public:
    StopList() {}

    template<typename Tp>
    StopList(Tp a[]) {
	int size = sizeof(a) / sizeof(Tp);
	StopList(a, a + size);
    }

    template<typename Iterator> 
    StopList(Iterator first, Iterator last) {
	while (first != last) map.insert(*first++);
    }

    bool has(const value_type& v) {
	container_type::iterator iter = map.find(v);
	return iter != map.end();
    }
};

}
#endif
