#ifndef MAGICBOX_H
#define MAGICBOX_H

#include <iostream>
#include <string>
#include <boost/program_options.hpp>

namespace magicbox {
    //--------------------------------------------------
    // (Plain) Option and TypedOption
    //-------------------------------------------------- 
    template<typename T> struct TypedOption {
	T* _t;
	const char* _spec;
	const char* _desc;
	int _pos;
    };

    struct Option {
	const char* _spec;
	const char* _desc;
    };

    struct RunOption {
	const char* _usage;
    };

    //--------------------------------------------------
    // This could come in handy: `$'
    //-------------------------------------------------- 
    template<typename T> TypedOption<T> 
    $(T* t, const char* spec, const char* desc, int pos = 0) {
	TypedOption<T> to = { t, spec, desc, pos }; return to;
    }

    Option $(const char* spec, const char* desc) {
	Option o = { spec, desc }; return o;
    }

    RunOption $$$(const char* usage = 0) { 
	RunOption ro = { usage }; return ro;
    }

    //--------------------------------------------------
    // Getopt
    //-------------------------------------------------- 
    namespace po = boost::program_options;

    class Getopt {
	int _argc;
	char** _argv;
	po::options_description* _visible;
	po::options_description* _hidden;
	po::positional_options_description* _pos;
	po::variables_map* _m;
	std::string _basename;
	std::string _usage;

	std::string basename(const char* filename) {
	    std::string name(filename);
	    size_t sep = name.rfind('/');
	    if (sep != std::string::npos) name.erase(0, sep + 1);
	    return name;
	}

    public:
	Getopt(int argc, char** argv) {
	    _argc = argc;
	    _argv = argv;
	    _visible = new po::options_description("Options");
	    _hidden = new po::options_description;
	    _pos = new po::positional_options_description;
	    _m = new po::variables_map;
	    _basename = basename(argv[0]);
	    _usage = _basename + " [option..]";

	    _visible->add_options()("help", "Show this help screen");
	}

	~Getopt() {
	    delete _visible;
	    delete _hidden;
	    delete _pos;
	    delete _m;
	}

	template<typename T> Getopt& operator<<(const TypedOption<T>& o) {
	    if (o._pos == 0) {
		if (o._t) _visible->add_options()(o._spec, po::value<T>(o._t), o._desc);
		else _visible->add_options()(o._spec, po::value<T>(), o._desc);
	    } else {
		_hidden->add_options()(o._spec, po::value<T>(o._t), o._desc);
		_pos->add(o._spec, o._pos);
	    }

	    return *this;
	}

	Getopt& operator<<(const Option& o) {
	    _visible->add_options()(o._spec, o._desc);
	    return *this;
	}

	Getopt& operator<<(const RunOption& o) {
	    try {
		po::options_description all;
		all.add(*_visible).add(*_hidden);
		po::store(po::command_line_parser(_argc, _argv).options(all).positional(*_pos).run(), *_m);
		po::notify(*_m);
	    } catch (po::invalid_option_value e) {}

	    if (o._usage) _usage = o._usage;
	    if (_m->count("help")) {
		print(std::cerr);
		exit(0);
	    }

	    return *this;
	}

	bool operator[](const char* key) { return _m->count(key); }

	std::ostream& print(std::ostream& o) {
	    o << "Usage: " << _usage << std::endl << std::endl;
	    return o << *_visible << std::endl; 
	}
    };

    //--------------------------------------------------
    // AutoStream
    //-------------------------------------------------- 
    template<class __Implementation, class __Interface, __Interface& __StandardStream> 
    class AutoStream {
	__Implementation stream;

    public:
	AutoStream<__Implementation, __Interface, __StandardStream>(std::string filename) {
	    if (filename != "-") stream.open(filename.c_str());
	}

	__Interface& operator()() {
	    return stream.is_open()? stream: __StandardStream;
	}
    };

    typedef AutoStream<std::ifstream, std::istream, std::cin> AutoIn;
    typedef AutoStream<std::ofstream, std::ostream, std::cout> AutoOut;
}

namespace std {
    std::ostream& operator<<(std::ostream& o, magicbox::Getopt& g) { return g.print(o); }
}

#endif
