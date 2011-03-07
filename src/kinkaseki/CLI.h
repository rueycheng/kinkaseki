#ifndef KINKASEKI_CLI_H
#define KINKASEKI_CLI_H

#include <iosfwd>
#include <string>
#include <vector>
#include <boost/program_options.hpp>
#include <boost/lexical_cast.hpp>

namespace kinkaseki {

namespace po = boost::program_options;

//--------------------------------------------------
// CLI
//-------------------------------------------------- 
class CLI {
    po::options_description general;
    po::options_description visible;
    po::options_description hidden;
    po::positional_options_description positional;
    po::variables_map vm;

protected:
    // Get the base part of a filename
    std::string basename(const char* filename) {
	std::string name(filename);
	size_t sep = name.rfind('/');
	if (sep != std::string::npos) name.erase(0, sep + 1);
	return name;
    }

public:
    CLI(): general("Options") {
	general.add_options()
	    ("help", "Show this help screen")
	    ("verbose", "Show verbose output");
    }

    // The most basic form of bind, supporting only set/unset operations
    // (No type specfication)
    CLI& bind(const char* spec, const char* desc = "") {
	visible.add_options()(spec, desc);
	return *this;
    }

    // Advanced form of bind, allowing user to bind an option to a variable
    template<typename T> CLI& bind(T& var, const char* spec, const char* desc = "") {
	std::string po_desc(desc);
	std::string value = boost::lexical_cast<std::string>(var);

	if (po_desc.size() > 0) po_desc += "\n";
	po_desc += "Defaults '" + value + "'";

	visible.add_options()(spec, po::value<T>(&var), po_desc.c_str());
	return *this;
    }

    // Commit the settings via boost framework
    std::vector<std::string> commit(int argc, char** argv, std::ostream& out = std::cerr) {
	std::string usage = "Usage: " + basename(argv[0]) + " [option..]";
	std::vector<std::string> passthrough;

	try {
	    po::options_description all;
	    all.add(general).add(visible).add(hidden);

	    po::parsed_options parsed = 
		po::command_line_parser(argc, argv).options(all).allow_unregistered().run();
	    po::store(parsed, vm);
	    po::notify(vm);

	    passthrough = 
		po::collect_unrecognized(parsed.options, po::include_positional);
	} catch (po::invalid_option_value e) {}

	if (vm.count("help")) {
	    out << usage << "\n";
	    out << general << "\n";
	    out << visible;
	    exit(0);
	}

	return passthrough;
    }

    // Tiny wrapper over variable_maps::count().  Comes in handy
    bool operator[](const char* key) {
	return vm.count(key);
    }
};

}

#endif
