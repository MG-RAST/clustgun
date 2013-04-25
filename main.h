//
//  Header.h
//  clustgun
//
//  Created by Wolfgang Gerlach on 5/9/12.
//  Copyright (c) 2012 University of Chicago. All rights reserved.
//

#ifndef clustgun_main_h
#define clustgun_main_h


// makes NDEBUG default:
#ifndef DEBUG
#define NDEBUG
#endif

#define USE_HashTableSimple

//#define USE_FLUSH_NEEDED //  avoiding flush made it 30% slower, just do the flush! 


#include <cstdlib>
#include <cstring>
#include <assert.h>
#include <map>
#include <set>
#include <cmath>

#include <iostream>
#include <fstream>
#include <algorithm>

#ifndef USE_HashTableSimple
#include <sparsehash/sparse_hash_map>
#endif

#include "olist.hpp"
#include "hat.hpp"
#include "kmer_iterator.hpp"
#include "fasta_parser.hpp"
#include "binarypath.hpp"
#include "read_blosum.hpp"

#include <omp.h>
#include "basic_tools_omp.hpp" 

#include <time.h>

#include <boost/program_options/detail/config_file.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>

#include <boost/iostreams/tee.hpp>
#include <boost/iostreams/stream.hpp>

#include <boost/dynamic_bitset.hpp>


#include "git_ref.h"


using namespace std;

#ifndef USE_HashTableSimple
using google::sparse_hash_map; 
#endif

namespace bio = boost::iostreams;
using bio::tee_device;
using bio::stream;



#define GCC_VERSION (__GNUC__ * 10000 \
	+ __GNUC_MINOR__ * 100 \
	+ __GNUC_PATCHLEVEL__)


#ifndef USE_HashTableSimple
#ifdef __APPLE__
using tr1::hash;// ext::hash;  // or __gnu_cxx::hash, or maybe tr1::hash, depending on your OS
#else
#if GCC_VERSION >= 40300
#include <unordered_map>
using std::unordered_map; // hash set was removed: http://gcc.gnu.org/gcc-4.3/changes.html
#else
using  __gnu_cxx::hash;
#endif
#endif
#endif

typedef short aminoacid;

typedef tee_device<ostream, ofstream> TeeDevice;
typedef stream<TeeDevice> TeeStream;

//stream<tee_device<ostream, ofstream > > my_split;
TeeStream  log_stream;


template <class T>
class mybasicstring;

typedef mybasicstring<char> mystring;


template <class T1, class T2, class T3> class triplet;




const int aminoacid_count = 21;
const char aminoacid_int2ASCII[aminoacid_count+1] = "ARNDCEQGHILKMFPSTWYV*"; // 20 AA according to wikipedia.. ;) plus stop codon "*"
																				// warning X can be represented either as -1 or aminoacid_count,
																				// normal values go from 0 to aminoacid_count-1
aminoacid aminoacid_ASCII2int[256];


// k-mer filtering, not used currently
const int low_abundance_threshold=5; // for kmers in the whole set of reads


// overlap detection:
int kmerlength;
int cluster_kmer_overlap_threshold;


// overlap verification:
int min_overlap_length;
//const double min_overlap_fraction  = 0.2; // an overlap length of 20% of the length of the shorter sequence is required
//int windowLength; // length of sliding window
//int windowScoreThreshold; // total BLOSUM score required for each window
double avgScoreThreshold; // average BLOSUM score required for an overlap
string blosum_file;
string prefixname;

//const int limit_input_reads = -1; // -1
//const int maxclustercount = 5000000;

const int max_protein_length = 65000;
//const int minimal_input_sequence_length = 10;





int thread_count;
int thread_count_requested;

//max_thread_count defined in basic_tools_omp.hpp

#ifdef DEBUG_DEADLOCK
int thread_status[max_thread_count]; //= {'0','0','0','0','0', '0','0','0','0','0'};
#endif
int reads_per_thread[max_thread_count]; //= {0,0,0,0,0 ,0,0,0,0,0};


typedef OList<int, short> cluster_member_list; // (Read, Offset)
typedef OList_rwlock<int, short> kmer_appearance_list; // (Cluster, Offset)


string string_int_2_kmer(int kmer_code);

bool fexists(string filename);





struct eqint
{
	//bool operator()(const char* s1, const char* s2) const
	//{
	//   return (s1 == s2) || (s1 && s2 && strcmp(s1, s2) == 0);
	//}
	bool operator()(int s1, int s2) const
	{
		return (s1 == s2);
	}
};

#ifdef DEBUG_DEADLOCK

inline void thread_update (int this_thread_id, int line) {
	#pragma omp critical(thread_status)
	{
	thread_status[this_thread_id]=line; // TODO not safe I think, should be critical section
	}
}
#define MACRO_THREAD_UPDATE(a,b) thread_update(a,b);

#else
#define MACRO_THREAD_UPDATE(a,b)
#endif

#ifdef USE_HashTableSimple
typedef class HashTableSimple HashTable;
#else
typedef class HashTableSparse HashTable;
#endif

class HashTableSimple  {
private:
	kmer_appearance_list * _hash;
	//vector<kmer_appearance_list > * _hash;
	size_t size;
public:
	
	//typedef vector<kmer_appearance_list * >::iterator vec_iterator;
	//typedef HashTableSimple_iterator iterator;
	//typedef kmer_appearance_list* iterator;
	
	//HashTableSimple() : ReaderWriterLock(thread_count) {
	HashTableSimple()  {
		size = (int) pow((double) aminoacid_count, (int) kmerlength);
		try {
			//_hash = new vector<kmer_appearance_list >(size,  kmer_appearance_list(thread_count));
			_hash = new kmer_appearance_list[size];
			//std::fill(*_hash, (*_hash)+size, kmer_appearance_list(thread_count));
			
			// init kmer_appearance_list(thread_count)
		} catch (bad_alloc& ba) {
			cerr << "error: (HashTableSimple) bad_alloc caught: " << ba.what() << endl;
			cerr << "parameter k was probably too big" << endl;
			exit(1);
		}
		
		for (size_t i = 0 ; i < size; ++i ) {
			_hash[i]=kmer_appearance_list(thread_count);
		}
		
	}
	
	~HashTableSimple() {
		delete [] _hash;
	}
	
	kmer_appearance_list* operator[] (size_t x) {
		return &(_hash[x]);
		//return &((*_hash)[x]);
	}
	
	kmer_appearance_list* at(size_t x) {
		
		if (x >= size) {
			cerr << "error: (HashTableSimple) x >= size" << endl;
			exit(1);
		}
		
		//return &(_hash->at(x));
		//return &((*_hash)[x]);
		return &(_hash[x]);
	}
	
	// should be optimized away
	inline void set_reader_lock(){};
	inline void unset_reader_lock(){};
	
	inline void set_writer_lock(){};
	inline void unset_writer_lock(){};
	
};

#ifndef USE_HashTableSimple
class HashTableSparse : public ReaderWriterLock, public sparse_hash_map<int, kmer_appearance_list * , hash<int>, eqint> {
	
		
	public:
	 
	HashTableSparse() : ReaderWriterLock(thread_count) {}
	
	
	
};
#endif

class Clustgun {
public:
	bool list_all_members;
	bool sort_input_seq;
	bool avgcov;
	string outputfile;
	string listfile;
	
	Clustgun() {
		list_all_members = false;
		sort_input_seq = false;
		avgcov = false;
	}
	
	void cluster(string inputfile);
	
	
};


template <typename T>
string NumberToString ( T Number )
{
	stringstream ss;
	ss << Number;
	return ss.str();
}


bool fexists(string filename)
{
	ifstream ifile(filename.c_str());
	return ifile;
}






template <class T>
class mybasicstring  {

private:
	T * data;
	size_t data_len;
	
public:
	mybasicstring(size_t seq_len){
		this->data_len = seq_len;
		this->data =new T[seq_len+1];
		this->data[this->data_len] = '\0';
	}
	
	mybasicstring(basic_string<T> * str) {
		this->data_len = str->length();
		data = new T[this->data_len+1];
		strcpy(data, str->c_str());
		data[this->data_len] = '\0';
	}
	

	mybasicstring() : data((T *) NULL), data_len(0) {}
	
	
	
	mybasicstring(char * seq, size_t seq_len){
		this->data = seq;
		this->data_len = seq_len;
	}
	
	
	
	mybasicstring(basic_string<T> str) {
		this->data_len = str.length();
		data = new T[this->data_len+1];
		strcpy(data, str.c_str());
		data[this->data_len] = '\0';
	}
	
	void delete_data () {
		if (data != NULL) {
			delete [] data;
		}
	}
	
	
	const size_t length(){
		return this->data_len;
	}
	
	const T * c_str() {
#ifdef DEBUG
		if (this->data == (T *) NULL) {
			cerr << "error: data == NULL" << endl;
			exit(1);
		}
#endif
		return data;
	}
	
	T& operator[] (size_t __pos)  {
		return data[__pos];
	}
	
	T& at (size_t __pos)  {
		if (__pos >= data_len) {
			cerr << "error: array out of bound, need to throw exeception.." << endl;
			exit(1);
		}
		return data[__pos];
	}
	

	
};



template <class T>
triplet<bool, T, int> majority_vote(vector<T> * vector_of_offsets, int number_of_elements) {
	// do a majority vote on the offsets, requires confirmation?	
	
	if (number_of_elements == 0) {
		return triplet<bool, T, int>(false, (T)0, 0);
	}
	
	int majority_vote_counter = 0;
	T current_vote_candidate= T();
	T e= T();
	for (int i = 0; i < number_of_elements; ++i) {
		#ifdef DEBUG
		e = vector_of_offsets->at(i);
		#else
		e = (*vector_of_offsets)[i];
		#endif
		
		//cout << "e: " << e<< endl;
		if (majority_vote_counter == 0) {
			current_vote_candidate = e; 
			majority_vote_counter = 1;
		} else {
			if (current_vote_candidate == e) {
				++majority_vote_counter;
			} else {
				--majority_vote_counter;
			}
		}
	}
	
	bool majority_found = false;
	
	if (majority_vote_counter > (number_of_elements/2)) {
		
		//fine
		majority_found = true;
	} else {
		// need to count
		int candidate_count = 0;
		for (int i = 0; i <= number_of_elements-1; ++i) {
			#ifdef DEBUG
			if (e == vector_of_offsets->at(i)) {
			#else
			if (e == (*vector_of_offsets)[i]) {
			#endif
				++candidate_count;
			}
		}
		if (candidate_count > (number_of_elements/2)+1 ) {
			majority_found = true;
		}
	}
	//if ((int) current_vote_candidate == 0) {
	//	cerr << "urks." << endl;
//		exit(1);
//	}
	
		
	// warning: majority_vote_counter is only a lower bound on the count of matching kmers!
	return triplet<bool, T, int>(majority_found,current_vote_candidate, majority_vote_counter);
}




#endif
