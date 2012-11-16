

#ifndef clustgun_kmeriterator_h
#define clustgun_kmeriterator_h


#include <iostream>
#include <string>
#include <cmath>


using namespace std;

typedef short aminoacid;


string string_int_2_kmer(int kmer_code, int kmerlength, const char * aminoacid_int2ASCII, int aminoacid_count);

class KmerIterator {
public:
	int kmerlength;
	const char * sequence;
	int seqlen;
	int kmer_start_pos; // refers to the first aa of the k-mer
	int code;
	int coded_length;
	
	const char* aminoacid_int2ASCII;
	aminoacid * aminoacid_ASCII2int;
	int aminoacid_count;
	
	KmerIterator(const char* seq, int seqlen, int startpos, int kmerlength, const char* aminoacid_int2ASCII, aminoacid* aminoacid_ASCII2int, int aminoacid_count);
	bool nextKmer();
	void reset(int startpos);
	
};

#endif
