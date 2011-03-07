#ifndef PORTER_H
#define PORTER_H
struct stemmer {
   char * b;       /* buffer for word to be stemmed */
   int k;          /* offset to the end of the string */
   int j;          /* a general offset into the string */
};

struct stemmer * create_stemmer(void);
void free_stemmer(struct stemmer * z);
int stem(struct stemmer * z, char * b, int k);
#endif
