#ifndef KINKASEKI_PRIORITY_QUEUE_H
#define KINKASEKI_PRIORITY_QUEUE_H

#include <vector>
#include <unordered_map>

namespace kinkaseki {

//--------------------------------------------------
// Home-made priority queue
//-------------------------------------------------- 
template<typename _Tp, typename _Sequence = std::vector<_Tp>, 
    typename _Map = std::unordered_map<_Tp, typename _Sequence::distance_type>,
    typename _Compare = std::less<typename _Sequence::value_type> >
class PriorityQueue {
public:
    typedef typename _Sequence::value_type value_type;
    typedef typename _Sequence::reference reference;
    typedef typename _Sequence::const_reference const_reference;
    typedef typename _Sequence::size_type size_type;
    typedef _Sequence container_type;

protected:
    _Sequence c;
    _Map m;
    _Compare comp;

    typedef typename _Sequence::iterator _Iterator;
    typedef typename _Sequence::difference_type _Distance;
    typedef typename _Sequence::value_type _Value;

    void __push_heap(_Iterator __first, _Distance __holeIndex,
		     _Distance __topIndex, _Tp __value)
    {
	_Distance __parent = (__holeIndex - 1) / 2;
	while (__holeIndex > __topIndex
	       && __comp(*(__first + __parent), __value)) {
	    *(__first + __holeIndex) = *(__first + __parent);
	    __holeIndex = __parent;
	    __parent = (__holeIndex - 1) / 2;
	}

	*(__first + __holeIndex) = __value;
    }

    void push_heap(_Iterator __first, _Iterator __last)
    {
	_Value __value = *(__last - 1);
	__push_heap(__first, _Distance((__last - __first) - 1),
		    _Distance(0), __value);
    }

    void __pop_heap(_Iterator __first, _Iterator __last, _Iterator __result)
    {
	_Value __value = *__result;
	*__result = *__first;
	__adjust_heap(__first, _Distance(0), _Distance(__last - __first), __value);
    }

    void pop_heap(_Iterator __first, _Iterator __last)
    {
	--__last;
	__pop_heap(__first, __last, __last);
    }


    void __adjust_heap(_Iterator __first, _Distance __holeIndex,
		     _Distance __len, _Tp __value)
    {
	const _Distance __topIndex = __holeIndex;
	_Distance __secondChild = __holeIndex;
	while (__secondChild < (__len - 1) / 2) {
	    __secondChild = 2 * (__secondChild + 1);
	    if (comp(*(__first + __secondChild),
		     *(__first + (__secondChild - 1))))
		__secondChild--;
	    *(__first + __holeIndex) = *(__first + __secondChild);
	    __holeIndex = __secondChild;
	}

	if ((__len & 1) == 0 && __secondChild == (__len - 2) / 2) {
	    __secondChild = 2 * (__secondChild + 1);
	    *(__first + __holeIndex) = *(__first + __secondChild - 1);
	    __holeIndex = __secondChild - 1;
	}

	__push_heap(__first, __holeIndex, __topIndex, __value);
    }

    void make_heap(_Iterator __first, _Iterator __last)
    {
	if (__last - __first < 2) return;

	const _Distance __len = __last - __first;
	_Distance __parent = (__len - 2) / 2;
	while (true) {
	    _Value __value = *(__first + __parent);
	    __adjust_heap(__first, __parent, __len, __value);
	    if (__parent == 0) return;
	    __parent--;
	}
    }

public:
    explicit PriorityQueue(const _Compare& __comp = _Compare(),
			   const _Sequence& __s = _Sequence(),
			   const _Map& __m = _Map())
    : c(__s), m(__m), comp(__comp) 
    { 
	make_heap(c.begin(), c.end(), comp); 
    }

    template<typename _InputIterator> 
    PriorityQueue(_InputIterator __first, _InputIterator __last,
		  const _Compare& __comp = _Compare(), 
		  const _Sequence& __s = _Sequence(),
		  const _Map& __m = _Map()) 
    : c(__s), m(__m), comp(__comp) 
    {
	c.insert(c.end(), __first, __last);
	make_heap(c.begin(), c.end(), comp);
    }

    bool empty() const { return c.empty(); }
    size_type size() const { return c.size(); }
    const_reference top() const { return c.front(); }

    void push(const value_type& __x) {
	c.push_back(__x);
	push_heap(c.begin(), c.end(), comp);
    }

    void pop() {
	pop_heap(c.begin(), c.end(), comp);
	c.pop_back();
    }
};

}

#endif
