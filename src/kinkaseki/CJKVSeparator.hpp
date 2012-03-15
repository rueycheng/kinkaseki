#ifndef KINKASEKI_CJKV_SEPARATOR_H
#define KINKASEKI_CJKV_SEPARATOR_H

#include <boost/tokenizer.hpp>
#include <boost/token_functions.hpp>

namespace kinkaseki {

class CJKVSeparator {
    template<typename InputIterator>
    short codesize(InputIterator next, InputIterator end) {
	if (next == end) return 0;
	if ((*next & 0x80) == 0x00) return 1;
	if ((*next & 0xE0) == 0xC0) 
	    return (next + 1 != end && (*(next + 1) & 0xC0) == 0x80)? 2: -2;
	if ((*next & 0xF0) == 0xE0)
	    return ((next + 1 != end && (*(next + 1) & 0xC0) == 0x80) && 
		    (next + 2 != end && (*(next + 2) & 0xC0) == 0x80))? 3: -3;
	if ((*next & 0xF8) == 0xF0)
	    return ((next + 1 != end && (*(next + 1) & 0xC0) == 0x80) && 
		    (next + 2 != end && (*(next + 2) & 0xC0) == 0x80) && 
		    (next + 3 != end && (*(next + 3) & 0xC0) == 0x80))? 4: -4;
	return -1;
    }

public:
    void reset() {}

    template<typename InputIterator, typename Token>
    bool operator()(InputIterator& next, InputIterator end, Token& tok) {
	using boost::tokenizer_detail::get_iterator_category;
	using boost::tokenizer_detail::assign_or_plus_equal;

	typedef get_iterator_category<InputIterator> getter;
	typedef assign_or_plus_equal<
	    typename getter::iterator_category> assigner;

	// while (next != end && (isspace(*next) || ispunct(*next))) ++next;
	while (next != end && isspace(*next)) ++next;
	if (next == end) return false;

	short size; 
	while ((size = codesize(next, end)) < 0) next += -size;
	if (size == 0) return false;

	InputIterator start = next;
	assigner::clear(tok); // NOTE: HACK!

	if (size == 1 && !ispunct(*start)) {
	    do {
		++next;
	    } while (codesize(next, end) == 1 && !isspace(*next) && !ispunct(*next));
	}
	else next += size;

        assigner::assign(start, next, tok); // NOTE: Hack!
	return true;
    }
};

}

#endif
