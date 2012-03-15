#ifndef KINKASEKI_FUNCTORS_H
#define KINKASEKI_FUNCTORS_H

namespace kinkaseki {

//--------------------------------------------------
// Choose1st
//-------------------------------------------------- 
template<typename Pair, typename Predicate>
struct Choose1st {
    Predicate pred;

    bool operator()(const Pair& x, const Pair& y) {
	return pred(x.first, y.first);
    }
};

template<typename Pair, typename Predicate>
Choose1st<Pair, Predicate> choose1st(Predicate pred) {
    return Choose1st<Pair, Predicate>();
}

//--------------------------------------------------
// Choose2nd
//-------------------------------------------------- 
template<typename Pair, typename Predicate>
struct Choose2nd {
    Predicate pred;

    bool operator()(const Pair& x, const Pair& y) {
	return pred(x.second, y.second);
    }
};

template<typename Pair, typename Predicate>
Choose2nd<Pair, Predicate> choose2nd(Predicate pred) {
    return Choose2nd<Pair, Predicate>();
}

}

#endif
