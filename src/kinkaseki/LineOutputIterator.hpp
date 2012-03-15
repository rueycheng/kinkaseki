#ifndef KINKASEKI_LINE_OUTPUT_ITERATOR_H
#define KINKASEKI_LINE_OUTPUT_ITERATOR_H

namespace kinkaseki {

struct LineOutputIterator:
    public std::ostream_iterator<std::string>
{
    LineOutputIterator(std::ostream& out):
	std::ostream_iterator<std::string>(out, "\n") {}
};

}

#endif
