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

#include <cstdlib>
#include <cstring>
#include <assert.h>
#include <map>
#include <set>
#include <cmath>

#include <iostream>
#include <fstream>
#include <algorithm>

#include <sparsehash/sparse_hash_map>

#include "olist.hpp"
#include "hat.hpp"
#include "kmer_iterator.hpp"
#include "fasta_parser.hpp"
#include "binarypath.hpp"
#include "read_blosum.hpp"


//#include <iomanip>
#ifdef TIME
#include <time.h>
#endif

#include <boost/program_options/detail/config_file.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>




using namespace std;
using google::sparse_hash_map; 

#define GCC_VERSION (__GNUC__ * 10000 \
	+ __GNUC_MINOR__ * 100 \
	+ __GNUC_PATCHLEVEL__)


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


typedef short aminoacid;

template <class T1, class T2, class T3> class triplet;

const int aminoacid_count = 21;
const char aminoacid_int2ASCII[aminoacid_count+1] = "ARNDCEQGHILKMFPSTWYV*"; // 20 AA according to wikipedia.. ;) plus stop codon "*"
aminoacid aminoacid_ASCII2int[256];


// k-mer filtering, not used currently
const int low_abundance_threshold=5; // for kmers in the whole set of reads


// overlap detection:
const int kmerlength = 5;
int cluster_kmer_overlap_threshold = 5;// 7


// overlap verification:
const int min_overlap_length = 10;
const double min_overlap_fraction  = 0.2; // an overlap length of 20% of the length of the shorter sequence is required
const int windowLength = 10; // length of sliding window
const int windowScoreThreshold = 0; // total BLOSUM score required for each window
const int avgScoreThreshold = 3; // average BLOSUM score required for an overlap
string blosum_file = "BLOSUM62";

string prefixname = "cluster";

//const int limit_input_reads = -1; // -1
//const int maxclustercount = 5000000;

const int max_protein_length = 65000;
const int minimal_input_sequence_length = 10;

const int hat_increase_steps = 1048576;





typedef OList<int, short> cluster_member_list; // (Read, Offset)
typedef OList<int, short> kmer_appearance_list; // (Cluster, Offset)



string string_int_2_kmer(int kmer_code);

bool fexists(string filename);


class Clustgun {
public:
	bool list_all_members;
	bool sort_input_seq;
	bool avgcov;
	string outputfile;
	
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
triplet<bool, T, int> majority_vote(vector<T> * vector_of_offsets, int number_of_elements) {
	// do a majority vote on the offsets, requires confirmation?	
	
	if (number_of_elements == 0) {
		return triplet<bool, T, int>(false, (T)0, 0);
	}
	
	int majority_vote_counter = 0;
	T current_vote_candidate;
	T e;
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
