#include <iostream>
#include <boost/foreach.hpp>
#include "compat.h"
#include "kinkaseki/CLI.h"

int main(int argc, char** argv) {
    kinkaseki::CLI cli(argc, argv);

    cli
	.bind("verbose", "Show verbose output")
	.setSynopsis("Segment input texts using the glue algorithm\n")
	.setTexts(
	    "  No concrete example so far.\n"
	);

    std::vector<std::string> args = cli.parse();

    BOOST_FOREACH (std::string& arg, args) {
	std::cout << arg << "\n";
    }

    return 0;
}
