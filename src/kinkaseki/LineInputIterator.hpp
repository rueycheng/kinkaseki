#ifndef KINKASEKI_LINE_INPUT_ITERATOR_H
#define KINKASEKI_LINE_INPUT_ITERATOR_H

#include <iostream>
#include <string>
#include <boost/iterator/iterator_facade.hpp>

namespace kinkaseki {

class LineInputIterator:
    public boost::iterator_facade<LineInputIterator, 
	const std::string, boost::forward_traversal_tag>
{
    std::istream* _in;
    std::string _value;
    bool _ok;

public:
    LineInputIterator(): 
	_in(0), _value(), _ok(false) {}

    explicit LineInputIterator(std::istream& in): 
	_in(&in), _value(), _ok(false) { getline(); }

private:
    friend class boost::iterator_core_access;

    void getline() {
	_ok = _in && *_in;
	if (_ok) 
	    _ok = std::getline(*_in, _value);
    }

    void increment() {
	getline();
    }

    const std::string& dereference() const {
	return _value;
    }

    bool equal(const LineInputIterator& other) const {
	return (_ok == other._ok) && (!_ok || _in == other._in);
    }
};

}

#endif
