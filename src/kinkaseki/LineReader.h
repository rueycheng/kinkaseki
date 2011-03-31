#ifndef KINKASEKI_LINEREADER_H
#define KINKASEKI_LINEREADER_H

#include <iosfwd>

namespace kinkaseki {

class LineReader {
    std::istringstream line;
    std::string buf;

public:
    std::istream& getline(std::istream& in) {
	if (std::getline(in, buf)) {
	    line.str(buf);
	    line.clear();
	    return line;
	}
	else return in;
    }
};

}

#endif
