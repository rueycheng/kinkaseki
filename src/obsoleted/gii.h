#ifndef GII_H
#define GII_H

//--------------------------------------------------
// Vocabulary
//-------------------------------------------------- 
template<typename WordID, class GrowPolicy>
class BasicVocabulary {
public:
    typedef WordID data_type;
    typedef GrowPolicy policy_type;
    typedef unordered_map<string, WordID> lookup_type;
    typedef vector<string> sequence_type;

    //--------------------------------------------------
    // Policy classes
    //-------------------------------------------------- 
    struct NoGrow {
	NoGrow(lookup_type&, sequence_type&) {}
	void grow(const string& word) {}
    };

    struct AutoGrow {
	lookup_type& l;
	sequence_type& s;

	AutoGrow(lookup_type& l, sequence_type& s): l(l), s(s) {}
	void grow(const string& word) { s.push_back(word); l[word] = s.size() - 1; }
    };

private:
    lookup_type lookup;
    sequence_type seq;
    policy_type pol;

    BasicVocabulary(): pol(lookup, seq) { seq.push_back("<unk>"); }

    template<typename Iterator> BasicVocabulary(Iterator first, Iterator last): BasicVocabulary() {
	copy(first, last, back_inserter(seq));
	data_type id = 1;
	while (first != last) lookup[*first++] = id++;
    }

    WordID size() { return seq.size() - 1; }

    WordID operator[](const string& word) {
	if (lookup.find(word) == string::npos) pol.grow(word);
	if (lookup.find(word) == string::npos) return 0;
	else return lookup[word];
    }
};

#endif
