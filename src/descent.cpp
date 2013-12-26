#include <iostream>
#include <vector>

#include "kinkaseki/CLI.hpp"
#include "kinkaseki/Lexicon.hpp"
#include "kinkaseki/LineReader.hpp"

int main(int argc, char** argv)
{
    using namespace std;

    kinkaseki::CLI cli(argc, argv);
    vector<string> args = cli.parse();

    typedef int Unigram;
    typedef pair<int, int> Bigram;

    // Shameless copy from ubt.cpp
    kinkaseki::Lexicon lexicon;
    const Unigram UNK = lexicon.encode("");
    const Unigram EOL = lexicon.encode("\n");

    vector<Unigram> text;
    kinkaseki::LineReader reader;

    while (istream& line = reader.getline(cin)) {
	string token;
	while (line >> token) 
	    text.push_back(lexicon.encode(token));
	text.push_back(EOL);
    }


    return 0;
}
