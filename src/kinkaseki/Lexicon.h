#ifndef KINKASEKI_LEXICON_H
#define KINKASEKI_LEXICON_H

#include <iosfwd>
#include <vector>

namespace kinkaseki {

template<typename Tp, 
    typename Integer = int,
    typename Container = unordered_map<Tp, Integer>, // FIXME
    typename Sequence = std::vector<Tp> >
class Lexicon_Impl {
protected:
    Container map;
    Sequence seq;

public:
    Lexicon_Impl() {}

    typename Sequence::size_type size() { 
	return seq.size();
    }

    Integer encode(const Tp& key) {
	typename Container::iterator iter = map.find(key);
	if (iter != map.end()) return iter->second;
	
	Integer value = seq.size();
	seq.push_back(key);
	map.insert(make_pair(key, value));
	return value;
    }

    Tp decode(const Integer& value) {
	return seq[value];
    }
};

typedef Lexicon_Impl<std::string> Lexicon;

}

#endif
