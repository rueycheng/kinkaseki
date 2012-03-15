#ifndef COMMON_H
#define COMMON_H

//--------------------------------------------------
// #include "expatpp.h"
// #include "expatpp_callback.h"
//-------------------------------------------------- 

#include <boost/algorithm/string.hpp>
#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/format.hpp>
#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH
#include <boost/random.hpp>
#include <boost/tokenizer.hpp>
#include <boost/token_functions.hpp>

//--------------------------------------------------
// Utilities
//-------------------------------------------------- 
struct pack {
    const char* target;
    short size;

    template<typename T> pack(const T& data) {
	target = reinterpret_cast<const char*>(&data);
	size = sizeof(T);
    }
};

struct unpack {
    char* target;
    short size;

    template<typename T> unpack(T& data) {
	target = reinterpret_cast<char*>(&data);
	size = sizeof(T);
    }
};

namespace std {
    ostream& operator<<(ostream& out, const pack& p) { return out.write(p.target, p.size); }
    istream& operator>>(istream& in, const unpack& p) { return in.read(p.target, p.size); }
}

//--------------------------------------------------
// Error reporting
//-------------------------------------------------- 
void say(const std::string& message) {
    std::cout << message << "\n";
}

void bye(const std::string& message) {
    std::cerr << message << "\n";
    exit(0);
}

void die(const std::string& message) {
    std::cerr << message << "\n";
    exit(1);
}

//--------------------------------------------------
// CJKVSeparator
//-------------------------------------------------- 
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

	typedef boost::tokenizer_detail::assign_or_plus_equal<
	    typename boost::tokenizer_detail::get_iterator_category<
	    InputIterator>::iterator_category
	    > assigner;

	while (next != end && (isspace(*next) || ispunct(*next))) ++next;
	if (next == end) return false;

	short size; 
	while ((size = codesize(next, end)) < 0) next += -size;
	if (size == 0) return false;

	InputIterator start = next;
	assigner::clear(tok); // NOTE: HACK!

	if (size == 1) {
	    do {
		++next;
	    } while (codesize(next, end) == 1 && !isspace(*next) && !ispunct(*next));
	}
	else next += size;

        assigner::assign(start, next, tok); // NOTE: Hack!
	return true;
    }
};

typedef boost::tokenizer<CJKVSeparator> CJKVTokenizer;

//--------------------------------------------------
// //--------------------------------------------------
// // TagContext
// //-------------------------------------------------- 
// class TagContext {
//     std::string context;
//     std::vector<int> offset;
// 
// public:
//     TagContext() {}
//     TagContext(const char* path) {
// 	context.assign(path);
// 	if (path[0] != '/') throw;
// 	for (int i = 0; path[i]; i++) 
// 	    if (path[i] == '/') offset.push_back(i + 1);
//     }
// 
//     bool operator==(const TagContext& t) { return context == t.context; }
//     bool operator!=(const TagContext& t) { return context != t.context; }
//     bool operator==(const char* path) { return context == path; }
//     bool operator!=(const char* path) { return context != path; }
// 
//     void enter(const char* comp) {
// 	offset.push_back(context.size());
// 	context = context + '/' + comp; // instead of `+='
//     }
// 
//     void leave() {
// 	if (offset.size()) {
// 	    context.erase(context.begin() + offset.back(), context.end());
// 	    offset.pop_back();
// 	}
//     }
// 
//     bool empty() const { return context.empty(); }
//     const char* path() const { return context.c_str(); }
//     const char* top() const { 
// 	if (offset.size()) return context.c_str() + offset.back() + 1;
// 	else return context.c_str();
//     }
// };
// 
// //--------------------------------------------------
// // ContextualSAXHandler
// //-------------------------------------------------- 
// class ContextualSAXHandler: public expatpp::callback::SAXHandler {
//     TagContext context;
// 
// public:
//     typedef std::map<std::string, std::string> attribute_map;
// 
//     void startElement(const expatpp::Char* name, const expatpp::Char** atts) { 
// 	if (strcmp(name, "root") == 0 && context.empty()) return;
// 	context.enter(name);
// 	enterContext(context.path(), name, make_map(atts));
//     }
// 
//     void endElement(const expatpp::Char* name) {
// 	if (context.path() == context.top()) return;
// 	leaveContext(context.path(), name);
// 	context.leave();
//     }
// 
//     void characterData(const expatpp::Char* s, int len) {
// 	if (context.empty()) return;
// 	text(context.path(), std::string(s, len));
//     }
// 
//     virtual void enterContext(std::string context, std::string name, attribute_map attr) {}
//     virtual void leaveContext(std::string context, std::string name) {}
//     virtual void text(std::string context, std::string text) {}
// 
// private:
//     inline attribute_map make_map(const expatpp::Char** atts) {
// 	attribute_map m;
// 
// 	while (*atts) {
// 	    m[*atts] = *(atts + 1);
// 	    atts += 2;
// 	}
// 	return m;
//     }
// };
//-------------------------------------------------- 

//--------------------------------------------------
// Procedures
//-------------------------------------------------- 

//--------------------------------------------------
// Strip
//-------------------------------------------------- 
// after: [i, end) 
// before [begin, i)
void strip_after_first(std::string& s, const char c) { 
    size_t cpos = s.find(c);
    if (cpos != std::string::npos) s.erase(cpos); 
}

void strip_before_first(std::string& s, const char c) { 
    size_t cpos = s.find(c);
    if (cpos != std::string::npos) s.erase(0, cpos); 
}

//--------------------------------------------------
// Comparator for std::pair
//-------------------------------------------------- 
struct first_cmp {
    template<typename T, typename S>
    bool operator()(const std::pair<T,S>& l, const std::pair<T,S>& r) {
	return l.first < r.first;
    }
};

struct first_rcmp {
    template<typename T, typename S>
    bool operator()(const std::pair<T,S>& l, const std::pair<T,S>& r) {
	return r.first < l.first;
    }
};

struct second_cmp {
    template<typename T, typename S>
    bool operator()(const std::pair<T,S>& l, const std::pair<T,S>& r) {
	return l.second < r.second;
    }
};

struct second_rcmp {
    template<typename T, typename S>
    bool operator()(const std::pair<T,S>& l, const std::pair<T,S>& r) {
	return r.second < l.second;
    }
};

//--------------------------------------------------
// Exception
//-------------------------------------------------- 
struct InvalidFormat {};
struct NotSupportedFunction {};
struct MissingArgument {};


#endif
