#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <cmath>
#include <numeric>

#include "magicbox.h"
#include "common.h"
#include "compat.h"

using namespace std;
using namespace magicbox;

//--------------------------------------------------
// TRECDocument
//-------------------------------------------------- 
struct TRECDocument {
    string docno;
    string sysno;
    vector<string> contents;
    vector<string> facets;

    void clear() {
	docno.clear();
	sysno.clear();
	contents.clear();
	facets.clear();
    }

    void normalize() {
	boost::trim(docno);
	boost::trim(sysno);
	foreach (string& content, contents) { boost::trim(content); }
	foreach (string& facet, facets) { boost::trim(facet); }
    }
};

namespace std {
    ostream& operator<<(ostream& o, const TRECDocument& doc) {
	unsigned int id = 0;

	foreach (const string& content, doc.contents) {
	    o << boost::format("%s#%02d") % doc.docno % ++id;
	    if (!doc.sysno.empty()) o << '@' << doc.sysno;
	    o << ' ' << content << "\n";
	}

	o << boost::format("%s:facet") % doc.docno;
	foreach (const string& facet, doc.facets) { o << ' ' << facet; }
	o << "\n";

	return o;
    }
}

//--------------------------------------------------
// TRECReader
//-------------------------------------------------- 
class TRECReader: public ContextualSAXHandler {
    TRECDocument doc;
    vector<TRECDocument> collection;

    istream& in;
    ostream& err;

    bool sentStart;
    bool sentEnd;

public:
    TRECReader(istream& in, ostream& err = std::cerr): in(in), err(err), sentStart(false), sentEnd(false) {}

    void enterContext(string context, string name, attribute_map attr) {
	boost::to_lower(context);
	if (context == "/doc") {
	    doc.docno = attr["docid"];
	    doc.sysno = attr["sysid"];
	}
	else if (context == "/doc/facet") doc.facets.push_back("");
	else if (context == "/doc/text/p") doc.contents.push_back("");
    }

    void leaveContext(string context, string name) {
	boost::to_lower(context);
	if (context == "/doc") {
	    doc.normalize();
	    collection.push_back(doc);
	    doc.clear();
	}
    }

    void text(string context, string text) {
	boost::to_lower(context);
	replace(text.begin(), text.end(), '\n', ' ');

	if (context == "/doc/docno") doc.docno += text;
	else if (context == "/doc/facet") doc.facets.back() += text;
	else if (context == ("/doc/text/p")) doc.contents.back() += text;
    }

    const vector<TRECDocument>& getCollection() { return collection; }
    vector<TRECDocument>::size_type size() { return collection.size(); }
    void clear() { collection.clear(); }

    bool findNewDocument() {
	char buf[BUFSIZE];

	string start = "<?xml version=\"1.0\" encoding=\"UTF-8\"?><root>";
	string end = "</root>";

	if (sentEnd) return false;

	if (!sentStart) {
	    parser.parse(start.c_str(), start.size(), false);
	    sentStart = true;
	}

	vector<TRECDocument>::size_type alreadyFound = collection.size();
	while (alreadyFound == collection.size() && in.read(buf, BUFSIZE))
	    parser.parse(buf, in.gcount(), false) || carp(err, buf);

	if (alreadyFound - collection.size() > 0) return true;

	parser.parse(buf, in.gcount(), false);

	if (!sentEnd) {
	    parser.parse(end.c_str(), end.size(), true);
	    sentEnd = true;
	}

	return alreadyFound - collection.size() > 0;
    }
};

//--------------------------------------------------
// Main program
//-------------------------------------------------- 
int main(int argc, char** argv) {

    // Getopt
    vector<string> input;

    Getopt g(argc, argv);
    g   << $(&input, "input", "", -1)
	<< $$$("files..");

    if (input.empty()) input.push_back("-");

    namespace fs = boost::filesystem;

    {
	foreach (const string& filename, input) {
	    AutoIn in(filename);
	    TRECReader reader(in());

	    while (reader.findNewDocument()) {
		foreach (const TRECDocument& doc, reader.getCollection()) { cout << doc << "\n"; }
		reader.clear();
	    }
	}
    }

    return 0;
}
