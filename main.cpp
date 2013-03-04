//
//  main.cpp
//  clustgun
//
//  Created by Wolfgang Gerlach on 5/9/12.
//  Copyright (c) 2012 University of Chicago. All rights reserved.
//


#include "main.h"


#ifdef DEBUG_DEADLOCK
void printThreadStatus(){
	
	cerr << "reporting threads is " << omp_get_thread_num() << endl;
	
#pragma omp critical(thread_status)
	{
		
		for (int tt = 0; tt < 10 ; tt++) {
			cerr << thread_status[tt] << ",";
			
		}
	}
#pragma omp critical(reads_per_thread)
	{
		
		cerr << endl;
		for (int tt = 0; tt < 10 ; tt++) {
			cerr << reads_per_thread[tt] << ",";
			
		}
		cerr << endl;
	}

	
}
#endif

struct eqstrXX
{
	//bool operator()(const char* s1, const char* s2) const
	//{
	 //   return (s1 == s2) || (s1 && s2 && strcmp(s1, s2) == 0);
	//}
	bool operator()(string s1, string s2) const
	{
		return (! s1.compare(s2));
	}
};

struct pair_second_arg {
	bool operator() (pair<int, int> i,pair<int, int> j) { return (i.second>j.second);}
} myobject;




string string_int_2_kmer(int kmer_code) {
	
  
	string kmerstring(kmerlength, '.');

	
	for (int i = 0 ; i < kmerlength-1; ++i){
		//cout << "decode: " << (kmer_code % aminoacid_count) << endl;
		kmerstring[i] = aminoacid_int2ASCII[kmer_code % aminoacid_count];
		kmer_code =  kmer_code / aminoacid_count;
	}
	//cout << "decode: " << (kmer_code % aminoacid_count) << endl;
	kmerstring[kmerlength-1] = aminoacid_int2ASCII[kmer_code % aminoacid_count];
	
	//cout << "got: " <<kmerstring << endl;
	
	return kmerstring;
}




// vector_of_offsets contains output !
int getOffsets(const char * sequence, int seq_len, vector<short> * vector_of_offsets, HashTable * cluster_kmer_hash, int max_cluster) {
	
	#ifndef USE_HashTableSimple
	HashTable::iterator cluster_kmer_hash_it;
	#endif
	
	KmerIterator * mykmerit;
	#ifdef DEBUG
	try {
#endif
		mykmerit = new KmerIterator(sequence, seq_len ,0, kmerlength, aminoacid_int2ASCII, aminoacid_ASCII2int, aminoacid_count);
#ifdef DEBUG
	} catch (bad_alloc& ba) {
		cerr << "error: (getOffsets) bad_alloc caught: " << ba.what() << endl;
		exit(1);
	}
#endif
	int number_of_overlapping_kmers_seen = 0;
	
	//cout << "search: " << *sequence << endl;
	kmer_appearance_list::iterator mylist_it;
	
	while (mykmerit->nextKmer()) {
		int code = mykmerit->code;
		short readpos = mykmerit->kmer_start_pos;
		
		#ifndef USE_HashTableSimple
		
		cluster_kmer_hash_it = cluster_kmer_hash->find(code);
		
		if (cluster_kmer_hash_it != cluster_kmer_hash->end()) {
		#endif
			#ifdef USE_HashTableSimple
			kmer_appearance_list * mylist = (*cluster_kmer_hash)[code];
			#else
			kmer_appearance_list * mylist = cluster_kmer_hash_it->second;
			#endif
			//cout << "---- "<< endl;
			//mylist->resetIterator();
			mylist->set_reader_lock();
			for (mylist_it = mylist->begin(); mylist_it != mylist->end(); ++mylist_it) {
			//while (mylist->nextElement()) {
				int clusterhit = mylist_it.getFirst();
							
				if (clusterhit == max_cluster) {
					number_of_overlapping_kmers_seen++;
					//cout << "match: " << string_int_2_kmer(code) << endl;
					short clusterpos = mylist_it.getSecond();
					
					//cout << "readpos: " << readpos << endl;
					//cout << "clusterpos: " << clusterpos << endl;
					//cout << "offset: " << clusterpos-readpos << endl;
					//	cout << "number_of_overlapping_kmers_seen: " << number_of_overlapping_kmers_seen << endl;
					
					#ifdef DEBUG
					try {
					vector_of_offsets->at(number_of_overlapping_kmers_seen-1)=clusterpos-readpos;
					} catch (out_of_range) {
						cerr << "error: (vector_of_offsets->at(number_of_overlapping_kmers_seen-1)=clusterpos-readpos) out_of_range" << endl;
						cerr << "number_of_overlapping_kmers_seen-1: " << number_of_overlapping_kmers_seen-1 << endl;
						cerr << "vector_of_offsets->size(): " <<vector_of_offsets->size() << endl;
						exit(1);
					}
					#else
					(*vector_of_offsets)[number_of_overlapping_kmers_seen-1]=clusterpos-readpos;
					#endif
				}
				
			}	   
			mylist->unset_reader_lock();
		#ifndef USE_HashTableSimple	
		} 
		#endif
		
		
	} // end while ( kmer iteration )
	delete mykmerit;
	return number_of_overlapping_kmers_seen; // vector_of_offsets contains output !
}

bool sortbyPairSecond( const pair<int, short>& i, const pair<int, short>& j ) {
    return j.second > i.second;
}

//bool sortbyPairSecondLength( const pair<string * , string * >& i, const pair<string * , string * >& j ) {
//    return j.second->length() > i.second->length();
//}


	
mystring * computeConsensus(cluster_member_list * mymemberlist, // cluster already locked
						  HashedArrayTree<pair<mystring, mystring > > * inputSequences,
						  vector<char> * alignment_column,
						  vector<bool> * aminoacid_occurence,
						  vector<int> * aminoacid_scores,
						  short * score_matrix,
						  int ** cluster_aminoacid_counts	) {
	
	
	// total procedure requires two passes of member list
	
	
	cluster_member_list::iterator mylist_it;
	
	
	int new_consensus_length = 0; // without base offset adjustment !
	int base_offset = INT_MAX;
	
	// need to adjust base offset (left-most member)
	// need to find right most member (offset + length)
	for (mylist_it = mymemberlist->begin(); mylist_it != mymemberlist->end(); ++mylist_it) {
		
		
		
		#ifdef DEBUG
		int read_id = mylist_it.getFirst();
		mystring * read_sequence = &(inputSequences->at(read_id).second);
		#else
		mystring * read_sequence = &((*inputSequences)[mylist_it.getFirst()].second);
		#endif
		
		
		short pos =  mylist_it.getSecond();
		
		if ((int)pos + (int)read_sequence->length() > new_consensus_length) {
			new_consensus_length = pos + read_sequence->length();
		}
		
		// make base_offset small
		if (pos < base_offset) {
			base_offset = pos;
		}
		
	}
	// adjust new_consensus_length
	new_consensus_length-=base_offset;
	
	#ifdef DEBUG
	if (new_consensus_length > max_protein_length) {
		cerr << "error: new_consensus_length > max_protein_length" << endl;
		exit(1);
	}
	#endif
	
	// init cluster_aminoacid_counts
	for (int i = 0 ; i < new_consensus_length ; ++i ) {
		int * blabla = cluster_aminoacid_counts[i];
		for (int j = 0; j < aminoacid_count ; ++j  ) {
			blabla[j]=0;
		}
	}
	
	
	
//	#ifdef DEBUG
//	if (sorted_members.size() == 0) {
//		cerr << "error: sorted_members.length() == 0" << endl;
//		exit(1);
//	}
//	#endif
	
	//int base_offset = sorted_members[0].second;
	//cout << "base_offset: " << base_offset << endl;
	mystring * new_consensus_seq;
	#ifdef DEBUG
	try {
	#endif
		new_consensus_seq = new mystring(new_consensus_length);
	#ifdef DEBUG
	} catch (bad_alloc& ba) {
		cerr << "error: (computeConsensus) bad_alloc caught: " << ba.what() << endl;
		exit(1);
	}
	#endif
	
	#ifdef DEBUG
	for (int i= 0 ; i < new_consensus_length; ++i) {
		
		(*new_consensus_seq)[i]='-'; // later we can detect something has not been written
	}
	#endif
	
	
	// update the real memberlist
	// update offsets (leftmost member has now offset zero)
	// AND put aminoacid counts into matrix
	for (mylist_it = mymemberlist->begin(); mylist_it != mymemberlist->end(); ++mylist_it) {
		//mylist_it.getSecond() -= base_offset;
		
		int read_id = mylist_it.getFirst();
		int offset = mylist_it.getSecond() - base_offset;
		mylist_it.getSecond() = offset;

		mystring * read_sequence = &((*inputSequences)[read_id].second);
		
		
		// go through aminoacids in member:
		for (size_t i = 0; i < read_sequence->length(); ++i ) {
			#ifdef DEBUG
			
			char aa_char = read_sequence->at(i);
			aminoacid aa = aminoacid_ASCII2int[(int)aa_char];
			
			
			if (offset+(int)i > new_consensus_length || offset+i  < 0 ) {
				cerr << "error: offset+i > new_consensus_length" << endl;
				cerr << "offset: " << offset  << endl;
				cerr << "i: " << i << endl;
				cerr << "new_consensus_length: " << new_consensus_length << endl;
				cerr << "max_protein_length: " << max_protein_length << endl;
				exit(1);
			}
			
			int * aa_array = cluster_aminoacid_counts[offset+i];
			
			if (aa >= aminoacid_count) {
				cerr << "aa >= aminoacid_count || aa < 0" << endl;
				cerr << (int) aa << endl;
				cerr << (int) aa_char << " " << aa_char << endl;
				exit(1);
			}
			if (aa > 0 ) { // known amino acid
				aa_array[aa]++;
			}
			#else
			aminoacid aa = aminoacid_ASCII2int[(int)(*read_sequence)[i]];
			if (aa > 0 ) { // known amino acid
				cluster_aminoacid_counts[offset+i][aa]++;
			}
			#endif
		}
		
	}
	
	
	// go through each column to find consenss
	for (int col_index = 0 ; col_index < new_consensus_length ; ++col_index ) {
		
		
		int * col = cluster_aminoacid_counts[col_index];
		
		// do a majority vote:
		aminoacid max_aa = 0;
		int max_aa_count = 0;
		int total_aa_count = 0;
		
		for (aminoacid this_aa = 0; this_aa < aminoacid_count ; ++this_aa  ) {
			int this_aa_count = col[this_aa];
			total_aa_count += this_aa_count;
			if (this_aa_count > max_aa_count) {
				max_aa = this_aa;
				max_aa_count = this_aa_count;
			}
		}
		
		if ((double)max_aa_count > (double)total_aa_count*0.5) {
			// majority found
			char c = aminoacid_int2ASCII[max_aa];
			(*new_consensus_seq)[col_index] = c;
			continue; // goto next column
		}
		
		// majority not found, try sum-of-pairs strategy
		//int max_score= 0;
		//max_aa = 0;
		
		int best_candiate_score=INT_MIN;
		aminoacid best_candidate=-1; // X unknown
		for (aminoacid candidate_aa = 0; candidate_aa < aminoacid_count; ++candidate_aa) {
			//(*aminoacid_occurence)[aa] = false;
			//(*aminoacid_scores)[aa] = 0;
			
			
			if ( col[candidate_aa] != 0 ) {
				int candidate_score = 0;
				int candidate_aa_256 = 256*candidate_aa;
				for (aminoacid this_aa = 0; this_aa < aminoacid_count; ++this_aa) {
					if ( col[this_aa] > 0 ) {
						candidate_score += ( col[this_aa] * score_matrix[this_aa+candidate_aa_256]  );
						
					}
				}
				
				if (candidate_score > best_candiate_score) {
					best_candiate_score = candidate_score;
					best_candidate = candidate_aa;
				}
			}
			
		}
		if (best_candidate > 0) {
			char c = aminoacid_int2ASCII[best_candidate];
			(*new_consensus_seq)[col_index] = c;
		} else {
			(*new_consensus_seq)[col_index] = 'X'; // do not like this, but can happen.
		}
	}

	#ifdef DEBUG
	for (int i= 0 ; i < new_consensus_length; ++i) {
			
		if ( (*new_consensus_seq)[i] == '-' ) {
			cerr << "error: (*new_consensus_seq)[i] == '-'" << endl;
			cerr << new_consensus_seq->c_str() << endl;
			exit(1);
		}
	}
	#endif
		
	return new_consensus_seq;
}

void addConsensusSequence(mystring * newconsensus, int cluster, HashTable * cluster_kmer_hash) {

//#ifdef DEBUG
//#pragma omp critical(cerr)
//	{
//		cerr << "thread: " << omp_get_thread_num() << " add cluster: " << cluster << " sequence: " << sequence << endl;
//	}
//#endif
	
	const char * sequence = newconsensus->c_str();
	int seq_len = newconsensus->length();
	
	#ifndef USE_HashTableSimple
	HashTable::iterator cluster_kmer_hash_it;
	#endif
	
	KmerIterator * mykmer_insertion_it;
	#ifdef DEBUG
	try {
#endif
		mykmer_insertion_it = new KmerIterator(sequence, seq_len, 0, kmerlength, aminoacid_int2ASCII, aminoacid_ASCII2int, aminoacid_count);
		#ifdef DEBUG
	} catch (bad_alloc& ba) {
		cerr << "error: (addConsensusSequence) bad_alloc caught: " << ba.what() << endl;
		exit(1);
	}
#endif
	//mykmer_insertion_it->reset(0);

	while (mykmer_insertion_it->nextKmer()) {
		int code = mykmer_insertion_it->code;
		int pos = mykmer_insertion_it->kmer_start_pos;
		
		
		
		#ifdef DEBUG
		// check if kmer is really existing
		
		string kmer = string_int_2_kmer(code, kmerlength, aminoacid_int2ASCII, aminoacid_count);

		char * ptrToSubString;
		ptrToSubString = strstr(sequence,kmer.c_str());
		
		if (ptrToSubString == NULL) {
			
			
			cerr << "kmer not found! " << endl;
			exit(1);
		}
		#endif
		
		
		kmer_appearance_list * mylist;
		
		#ifndef USE_HashTableSimple
		cluster_kmer_hash->set_reader_lock();
		
		cluster_kmer_hash_it = cluster_kmer_hash->find(code);
		#endif
		
		#ifdef USE_HashTableSimple

		mylist = (*cluster_kmer_hash)[code];
	
		#else
		
		if (cluster_kmer_hash_it != cluster_kmer_hash->end()) {
		
			//found list
			
			mylist = cluster_kmer_hash_it->second;
			
			cluster_kmer_hash->unset_reader_lock();
			
			//cout << "append to old list: " << string_int_2_kmer(code) << endl;
			//exit(0);
			
		} else {
			//did not find list
			#ifdef DEBUG
			try {
			#endif
				mylist = new kmer_appearance_list(thread_count);
			#ifdef DEBUG
			} catch (bad_alloc& ba) {
				cerr << "error: (new kmer_appearance_list) bad_alloc caught: " << ba.what() << endl;
				exit(1);
			}		
			#endif
			
			
			
			cluster_kmer_hash->unset_reader_lock(); // to avoid dead lock
			
			cluster_kmer_hash->set_writer_lock();
			
			(*cluster_kmer_hash)[code]=mylist; // that line is independent of which hash table I use..
			
			cluster_kmer_hash->unset_writer_lock();
			
			
			
			
		}
		#endif
		mylist->set_writer_lock();
		mylist->append(cluster,pos);
		mylist->unset_writer_lock();
		
		//if (code == 3566)
		//	cout << "insertX: "<< clusterhit << " " << pos << endl;
		//if (code == 1579170)
		//	cout << "insert: "<< clusterhit << " " << pos << endl;
	}// edn while

	delete mykmer_insertion_it;
}

// removes kmers from index
void removeConsensusSequence(mystring * consensus_old, int cluster, HashTable * cluster_kmer_hash) {

	//#ifdef DEBUG
	//#pragma omp critical(cerr)
	//{
	//cerr << "thread: " << omp_get_thread_num() << " remove cluster: " << cluster << " consensus_old: " << consensus_old->c_str() << endl;
	//}
	//#endif
	
	#ifndef USE_HashTableSimple
	HashTable::iterator cluster_kmer_hash_it;
	#endif
	
	KmerIterator * mykmer_deletion_it;
	#ifdef DEBUG
	try {
#endif
		mykmer_deletion_it = new KmerIterator(consensus_old->c_str(), consensus_old->length(), 0, kmerlength, aminoacid_int2ASCII, aminoacid_ASCII2int, aminoacid_count);
		#ifdef DEBUG
	} catch (bad_alloc& ba) {
		cerr << "error: (removeConsensusSequence) bad_alloc caught: " << ba.what() << endl;
		exit(1);
	}	
#endif
	//mykmer_deletion_it->reset(0);
	
	//if (cluster == 606) {
	//	cout << "delete cluster 606   --------------------------------------------" << endl;
		
	//}
	
	while(mykmer_deletion_it->nextKmer()) {
		
		int code = mykmer_deletion_it->code;
		//int pos = mykmer_deletion_it->kmer_start_pos;
		
		//cout << "searchdel: " << code << " in  " << clusterhit << endl;
		
		//if (cluster == 606 && code == 2644139) {cout << "delete: " << string_int_2_kmer(code) << endl;}
		
		

		kmer_appearance_list * mylist;
		kmer_appearance_list::iterator mylist_it;

		#ifndef USE_HashTableSimple
		cluster_kmer_hash->set_reader_lock();
		#endif
		
		#ifdef USE_HashTableSimple
		if ((*cluster_kmer_hash)[code]!=NULL) {
		#else
		cluster_kmer_hash_it = cluster_kmer_hash->find(code);
		if (cluster_kmer_hash_it != cluster_kmer_hash->end()) { // kmer-list exists
		#endif
			
			#ifdef USE_HashTableSimple
			mylist = (*cluster_kmer_hash)[code];
			#else
			mylist = cluster_kmer_hash_it->second;
			#endif
			
			#ifndef USE_HashTableSimple
			cluster_kmer_hash->unset_reader_lock();
			#endif
			//mylist->resetIterator();
			#ifdef DEBUG
			bool found_kmer = false;
			#endif
			
			//while (mylist->nextElement()) {
			
			//cerr << "start for loop" << endl;
			mylist->set_writer_lock();
			for (mylist_it = mylist->begin(); mylist_it != mylist->end(); ++mylist_it) {
				//if (cluster == 606 && code == 2644139) {cout << "mylist->currentArrayPosition: " << mylist->currentArrayPosition << " mylist->lastArrayPosition: " << mylist->lastArrayPosition << " mylist->getFirst(): "<< mylist->getFirst() << endl; }
				
				//cerr << mylist_it.getFirst() << " " << cluster << endl;
				if (mylist_it.getFirst() == cluster) {
					//if (cluster == 606 && code == 2644139) {mylist->print();}
					#ifdef DEBUG
					found_kmer = true;
					#endif
					mylist->erase(mylist_it);
					//if (cluster == 606 && code == 2644139) {cout << "after erase element" << endl; mylist->print();}
					//if (mylist->getLength() > 1) {
					//exit(0);
					//}
					break;
				}
				
			}
			mylist->unset_writer_lock();
			
			#ifdef DEBUG
			if (! found_kmer) {
				cerr << "error: did not find kmer I wanted to delete! cluster: " << cluster << endl << " kmer: " << code << " \"" << string_int_2_kmer(code) << "\""<< endl;
				cerr << "string: " << consensus_old->c_str() << endl;
				mylist->print();
				exit(1);
			}
			#endif
			//cout << "append to old list: " << string_int_2_kmer(code) << endl;
			//exit(0);
			
		} else {
			cerr << "error: kmer for deletion not found! cluster "<< cluster << " kmer:" << string_int_2_kmer(code) << endl;
			cerr << "string: " << consensus_old->c_str() << endl;
			exit(1);
			//cout << "add new list" << endl;
		}
		//mylist->append(last_cluster,pos);
	} // end while
	
	delete mykmer_deletion_it;

	return;
}

void testConsensusSequence(mystring * consensus_old, int cluster, HashTable * cluster_kmer_hash, int pos_in_source) {
	
	
	
	#ifndef USE_HashTableSimple
	HashTable::iterator cluster_kmer_hash_it;
	#endif
	KmerIterator * mykmer_deletion_it;
#ifdef DEBUG
	try {
#endif
		mykmer_deletion_it = new KmerIterator(consensus_old->c_str(), consensus_old->length(), 0, kmerlength, aminoacid_int2ASCII, aminoacid_ASCII2int, aminoacid_count);
#ifdef DEBUG
	} catch (bad_alloc& ba) {
		cerr << "error: (removeConsensusSequence) bad_alloc caught: " << ba.what() << endl;
		exit(1);
	}
#endif
	
	while(mykmer_deletion_it->nextKmer()) {
		
		int code = mykmer_deletion_it->code;
				
		
		
		kmer_appearance_list * mylist;
		kmer_appearance_list::iterator mylist_it;
		
		#ifndef USE_HashTableSimple
		cluster_kmer_hash->set_reader_lock();
		#endif
		
		#ifdef USE_HashTableSimple
		if (true) {
			mylist = (*cluster_kmer_hash)[code] ;
		#else
		cluster_kmer_hash_it = cluster_kmer_hash->find(code);

		if (cluster_kmer_hash_it != cluster_kmer_hash->end()) { // kmer-list exists
			
			mylist = cluster_kmer_hash_it->second;
		#endif
			#ifndef USE_HashTableSimple
			cluster_kmer_hash->unset_reader_lock();
			#endif
#ifdef DEBUG
			bool found_kmer = false;
#endif
			
			
			mylist->set_reader_lock();
			for (mylist_it = mylist->begin(); mylist_it != mylist->end(); ++mylist_it) {
				
				if (mylist_it.getFirst() == cluster) {
					
#ifdef DEBUG
					found_kmer = true;
#endif
					
					
					break;
				}
				
			}
			mylist->unset_reader_lock();
			
#ifdef DEBUG
			if (! found_kmer) {

				#pragma omp critical(cerr)
				{

				cerr << "error(testConsensusSequence): did not find kmer! cluster: " << cluster << endl << " kmer: " << code << " \"" << string_int_2_kmer(code) << "\""<< endl;
				cerr << "string: " << consensus_old->c_str() << " len:"<< consensus_old->length() << endl;
				cerr << "pos_in_source: " << pos_in_source << endl;
					printThreadStatus();
				}
				mylist->print();
				exit(1);
				
			}
#endif
			
			
		} else {
#pragma omp critical(cerr)
			{
			cerr << "error(testConsensusSequence): kmer not found! cluster "<< cluster << " kmer:" << string_int_2_kmer(code) << endl;
			cerr << "string: " << consensus_old->c_str() << " len:"<< consensus_old->length() << endl;
			cerr << "pos_in_source: " << pos_in_source << endl;
				#ifdef DEBUG_DEADLOCK
				printThreadStatus();
				#endif
			exit(1);
			}
		}
		
	} // end while
	
	delete mykmer_deletion_it;
	return;
}


bool computeSequenceOverlap(int offset, mystring * a, mystring * b, short * score_matrix, int majority_vote_count) {
	int start_i = max(0, offset);
	int start_j = start_i - offset;
	
	//cout << "i: " << i << endl;
	//cout << "j: " << j << endl;
	
	if (false) {
		//int shift = 0;
		//if (j<0) {
		//shift = -j;
		//}
		cout << string(start_j, ' ') << a->c_str() << endl;
		cout << string(start_i, ' ') << b->c_str() << endl;
	}
	
	//if (j<0) {
	//	exit(0);	
	//}
	
	
	
	int a_len = a->length();
	int b_len = b->length();
	
	int overlap_length = min((a_len - start_i), (b_len - start_j));
	//cout << "overlap_length: " << overlap_length << endl; 
	
	
	if (overlap_length < min_overlap_length) {
		return false;
	}
	
	//if majority_vote_count indicates high similarity
	if (majority_vote_count >= (int)(overlap_length-kmerlength+1)*0.8) {
		return true;
	}
	
	
	//int min_required_length = (int) ((double)min(a_len, b_len)*(double)min_overlap_fraction);
	//if (overlap_length < min_required_length) {
		//cout << "failed ==============================" << endl;
	//	return false;
	//}
	
	
	
	
	//exit(0);
	int total_score = 0;
	int tot_len = 0;
	
	int aa_score;
	
	//int currentWinLen = 0;
	//int windowScore = 0;
	
	
	//int start_i = i;
	int end_i = start_i+overlap_length;
	for (int i = start_i; i< end_i ; ++i) {
		tot_len++;
		int j = i - offset;
		
		//cout << "i " << (*a)[i] << " "  << (*b)[j] << " " << i << " " << j  << endl;
		
		
		
#ifdef DEBUG
		//int index;
		char a_char;
		char b_char;
		try {
				
			a_char = a->at(i);
			b_char = b->at(j);
		} catch (out_of_range& oor) {
			cerr << "Out of Range error:(index = a->at(i)+256 * b->at(j)) " << oor.what() << endl;
			exit(1);
		}
		
		aa_score = getScoreSave(a_char, b_char, score_matrix );
#else
		aa_score = getScore((*a)[i],(*b)[j], score_matrix );
#endif	
		
		//aa_score = score_matrix[index];
		
		total_score += aa_score;
		//windowScore += aa_score;
		
//		if (currentWinLen < windowLength) {
//			//cout << "_windows score is : " << windowScore << endl;
//			currentWinLen++;
//		} else {
//			
//			//remove first aa score, of aa that has left window
//			#ifdef DEBUG
//			char a_char;
//			char b_char;
//			try {
//				a_char = a->at(i-windowLength);
//				b_char = b->at(j-windowLength);
//			} catch (out_of_range& oor) {
//				cerr << "Out of Range error:(index = a->at(i-windowLength)+256 * b->at(j-windowLength)) " << oor.what() << endl;
//				exit(1);
//			}
//			
//			//if (index >= 256*256) {
//			//	cerr << "B) index >= 256*256: " << index << endl;
//		//		exit(1);
//			//}
//			//if (index < 0) {
//			//	cerr << "B) index < 0: " << index << endl;
//			//	exit(1);
//			//}
//			aa_score = getScoreSave(a_char, b_char, score_matrix );
//			#else
//			aa_score = getScore( (*a)[i-windowLength], (*b)[j-windowLength], score_matrix );
//			#endif
//			
//			
//			windowScore -= aa_score;
//			
//			
//			//cout << "windows score is : " << windowScore << endl;
//			if (windowScore < windowScoreThreshold) {
//				//cerr << "windowScore is too small" << endl;
//				return false;
//				//exit(1);
//				
//			}
//		}
				
				
		
		//cout << "score: " << score << endl;
		
		
		
	}
	
	#ifdef DEBUG
	if (tot_len <= 0) {
		cerr << "error: tot_len <= 0 : " << tot_len << endl;
		exit(1);
	}
	#endif
	
	double avg_score = (double)total_score/(double)tot_len;
	//cout << "avg_score: " << avg_score << endl;
	if (avg_score < avgScoreThreshold) {
		//cout << "avgScore is too small ---------------------------------" << endl;
		//exit(1);

		return false;
	
	}
		
	return true;
}


bool clusterNeedsUpdate(int size) {
	if (size <= 15) {
	//if ((size <= 15) || (size % 5 == 0)) {
		return true;
	}
	return false;
}


void sortInputSequences(HashedArrayTree<pair<mystring *, mystring * > > * inputSequences) {
	cerr << "sort sequences..." << endl;
	//sort(inputSequences->begin(), inputSequences->end(), sortbyPairSecondLength); // would be nice, but I have no iterators...
	
	int seq_count = inputSequences->size();
	
	//for (int x = 0; x < seq_count; ++x) {
	//	cout <<(*inputSequences)[x].second->length() << endl;
	//	
	//}
	//cout << " --- " << endl;
	
	for (int i = seq_count-1; i >= 0; --i) {
		//cout << " i = " << i << endl;
		for ( int j = 0; j < i; ++j){
			//cout << " j = " << j << endl;
			if ( (*inputSequences)[j].second->length() < (*inputSequences)[j+1].second->length() ) {
				//swap((*inputSequences)[i], (*inputSequences)[i+1]);
				pair<mystring * , mystring * > pp = (*inputSequences)[j];
				(*inputSequences)[j] = (*inputSequences)[j+1];
				(*inputSequences)[j+1] = pp;
				
			}
			
			//if ( (*inputSequences)[j].second->length() < (*inputSequences)[j+1].second->length() ) {
			//	cerr << i << endl;
			//	exit(1);
			//}
		}
		//for (int x = 0; x < seq_count; ++x) {
		//	cout <<(*inputSequences)[x].second->length() << endl;
		
		//}
		//cout << " --- " << endl;
		
	}
	
	cerr << "... sorting done." << endl;
	
	
	//for (int i = 0; i < seq_count-1; ++i) {
	//	cout <<(*inputSequences)[i].second->length() << endl;
	//	if ( (*inputSequences)[i].second->length() < (*inputSequences)[i+1].second->length() ) {
	//		cerr << "error: at " << i << " " << (*inputSequences)[i+1].second->length() << endl;
	//		exit(1);
	//	}
	//}

	
}


void Clustgun::cluster(string inputfile) {
	double vm, rss;

	//log_stream << "huhu" << endl;	
  
	log_stream<< "clustgun starts...\n";

	string date = exec("date");
	log_stream << date << endl;
	
	
	
	
	
	time_t begin, end; 
	time(&begin);
	
	#ifdef TIME
	timespec kmersearch_total = {0, 0};
	timespec kmersearch_start; 
	timespec kmersearch_end; 
	
	timespec validation_total = {0, 0};
	timespec validation_start;
	timespec validation_end;
	
	timespec overlap_total = {0, 0};
	timespec overlap_start;
	timespec overlap_end;
	#endif
	
	
	// create mappings of amino acid ASCII to int and back.
	
	log_stream << "initialize alphabet mapping... \"" << aminoacid_int2ASCII << "\""<< endl;

	// init with -1
	std::fill_n(aminoacid_ASCII2int, 256, -1);
	
	for (aminoacid i = 0 ; i<aminoacid_count; ++i) {
		char aa = aminoacid_int2ASCII[i];
		//cout << aa << " " << i << endl;
		aminoacid_ASCII2int[(int) aa]=i;
		aa = tolower ( aa );
		aminoacid_ASCII2int[(int) aa]=i;
		//cout << aa << " " << i << endl;
	}
  
	
	//for (int i = 0 ; i<255; ++i) {
	//	cout << i << ": " << (char) i  << " " << aminoacid_ASCII2int[i] << endl;
	//}
	
	
	// search BLOSUM file
	string blosum_file_confirmed;
	bool matrix_exists = fexists(blosum_file);
	
	if (matrix_exists) {
		
		blosum_file_confirmed = blosum_file;
	} else {
		
		size_t found=blosum_file.find('/');
		if (found != string::npos) {
			// already seems to an abosulte path name
			cerr << "error: (a) blosum file not found: " << blosum_file << endl;
			exit(1);
		} 
		
		// try to find blosum file in directory of the binary
		string newfilename = getbinarypath();
		newfilename.append("/");
		newfilename.append(blosum_file);
		//cout << "try: " << newfilename << endl;
		matrix_exists = fexists(newfilename);
		
		if (!matrix_exists) {
			cerr <<  "error: (b) blosum file not found: " << blosum_file << endl;
			exit(1);
		}
		blosum_file_confirmed = newfilename;
	}
	
	short * score_matrix = readBLOSUM(blosum_file_confirmed.c_str());
	
	
	// ------------------------------------------------------------
	// ----------------------START---------------------------------
	// ------------------------------------------------------------
	
	//cout << argc << endl;
	
	
	string file (inputfile);
	//cout << file << endl;
	//exit(0);
	
	//string file ("/Users/wolfganggerlach/protassembly/protassembly/protassembly/test.fas");
	//string * zcat_command = new string("/usr/bin/gzcat");	  // zcat expects file with .Z extension !?			   
   // string cat_command = ("cat /Users/wolfganggerlach/protassembly/protassembly/protassembly/test.fas");
	
#ifdef __APPLE__
	string * zcat_command = new string("/usr/bin/gzcat"); 
#else
	string * zcat_command = new string("/bin/zcat"); 
#endif


	
	#ifndef USE_HashTableSimple
	sparse_hash_map<int, int, hash<int>, eqint> kmer_count_hash;
	sparse_hash_map<int, int, hash<int>, eqint>::iterator it;
	#endif
	
	FASTA_Parser * fasta_parser;
	
	#ifndef USE_HashTableSimple
	int total_substrings=0;
	#endif
	
	size_t total_read_count = 0;
	string descr;
	string * fasta_sequence;

	// ----------------------------------------
	// read sequences into memory
	
	HashedArrayTreeString * inputSequencesData = new HashedArrayTreeString(23); // 2^23 = 8MB chunks
	
	
	HashedArrayTree<pair<mystring, mystring > > * inputSequences;
	#ifdef DEBUG
	try {
	#endif
		inputSequences = new HashedArrayTree<pair<mystring, mystring > >(false, 20); // 20 for 2^20=1MB chunks
		#ifdef DEBUG
	} catch (bad_alloc& ba) {
		cerr << "error: (inputSequences) bad_alloc caught: " << ba.what() << endl;
		exit(1);
	}
	#endif
	
	#ifdef DEBUG
	inputSequences->name = string("inputSequences");
	#endif
	
	fasta_parser = new FASTA_Parser(file, true, zcat_command);
	//string my_sequence;
	int ignore_short_seq = 0;
	while (fasta_parser->getNextDescriptionLine(descr)) {
		
		//cerr << "descr: " << descr<< endl;
		fasta_sequence = fasta_parser->getSequence();
		//my_sequence = *sequence;
		
		
		
		if (descr.length() > 0 && (int) fasta_sequence->length() > min_overlap_length) { // won't make sense to keep sequences short as the min_overlap_length
			//cout << "huhu: " << *sequence << endl;
			
			total_read_count++;
			
			if (total_read_count < 10) {
				
				size_t found;
				int dna_amount = 0;
				found=fasta_sequence->find_first_of("ACGTNactgn");
				while (found != string::npos)
				{
					dna_amount++;
					found = fasta_sequence->find_first_of("ACGTNactgn",found+1);
				}
				
				if (dna_amount > 0.9*(fasta_sequence->length()) ) {
					
					cerr << "sequence: " << *fasta_sequence << endl;
					cerr << "#bases: " << dna_amount << " length: " << fasta_sequence->length() << endl;
					
					cerr << "error: this looks like DNA, clustgun can assemble only protein sequences" << endl;
					exit(1);
				}
			}
			
			char * descr_pointer = inputSequencesData->addSequence(descr.c_str(), descr.length());
			
			char * seq_pointer = inputSequencesData->addSequence(fasta_sequence->c_str(), fasta_sequence->length());
			
			// idea:
			// use HATSequence->addData to store concatenation of sequence and description
			// point only to sequence
			// I may store the length of the sequences in fixed length array afterwards.  Or better store length in first 1 or 2 bytes in front of sequence.
			
			
			#ifdef DEBUG
			try {
			#endif
				
				inputSequences->push_back(pair<mystring, mystring >(mystring(descr_pointer, descr.length()), mystring(seq_pointer, fasta_sequence->length())));
				
			#ifdef DEBUG
			} catch (bad_alloc& ba) {
				cerr << "error: (inputSeuqences->push_back) bad_alloc caught: " << ba.what() << endl;
				exit(1);
			}
			#endif
			
			
			
			//cout << "pushed: " << inputSequences->size() << endl;	
			//cout << *sequence << endl;
			//exit(0);
							
			   
			if (total_read_count % 1000000 == 0 ) {
				cerr << "total_read_count: " << total_read_count << endl;
			}
						   
			
		} else {
			ignore_short_seq++;
			//cerr << "warning: sequence was not accepted..." << endl;
			
		}
		delete fasta_sequence;
		
		//if (total_read_count >= limit_input_reads && limit_input_reads != -1) {
		//	cerr << "WARNING: numer of reads for debugging purposes limited !!!!!!" << endl;
		//	break;
		//}
	}
	
	if (ignore_short_seq > 0 ) {
		log_stream << "warning: ignored " << ignore_short_seq << " sequences, because they were shorter than the mininmal overlap length." << endl;
	}
	
	if (total_read_count % 1000000 != 0 || total_read_count == 0 ) {
	cerr << "total_read_count: " << total_read_count << endl;
	}

	
	//cout << "------------- " << *(inputSequences->at(0).second) << endl;
// count k-mers:	
	
//	KmerIterator * mykmerit = new KmerIterator(sequence, 0, kmerlength);
//	
//	while (mykmerit->nextKmer()) {
//		int code = mykmerit->code;
//		total_substrings++;
//		kmer_count_hash[code]++;
//		
//	}
//	
//	delete mykmerit;	
	
	
	delete fasta_parser;
	fasta_parser = 0;
	
	
	// SORT input sequences
	
	//if (sort_input_seq) {
	//	sortInputSequences(inputSequences);
	//}
	
	
	
	
	#ifndef USE_HashTableSimple
	bool kmerabundance = false;
	
	if (kmerabundance) {
		vector< pair<int, int> > * array;
		vector< pair<int, int> >::iterator array_it;
		
		int sum;
		int totsum;
		
		#ifdef DEBUG
		try {
		#endif
			array = new vector< pair<int, int> >;
			#ifdef DEBUG
		} catch (bad_alloc& ba) {
			cerr << "error: (array vector) bad_alloc caught: " << ba.what() << endl;
			exit(1);
		}
		#endif
		// print hash table / put in array
		for ( it=kmer_count_hash.begin() ; it != kmer_count_hash.end(); it++ ) {
			if ((*it).second > low_abundance_threshold) {
				// cout << string_int_2_kmer((*it).first) << " => " << (*it).second << endl;
				#ifdef DEBUG
				try {
				#endif
					array->push_back(pair<int, int >((*it).first, (*it).second));
				#ifdef DEBUG
				} catch (bad_alloc& ba) {
					cerr << "error: (array->push_back) bad_alloc caught: " << ba.what() << endl;
					exit(1);
				}	
				#endif
			}
		}
		
		sort (array->begin(), array->end(), myobject);
		
		//for ( array_it=array->begin() ; array_it != array->end(); array_it++ ) {
		
		//	   cout << (*array_it).first << " => " << (*array_it).second << endl;
		
		//}
		
		cout << "array->size(): " << array->size() << endl;
		
		cout << "total_substrings: " << total_substrings << endl;
		
		// print kmer-abundance profile
		sum = 0;
		totsum = 0;
		int kmer_abundance_stepsize = 1000;
		cout << "kmer abundance profile:" << endl;
		for (int i=1; i<array->size(); ++i) {
			totsum += array->at(i).second;
			if (i % kmer_abundance_stepsize == 0) {
				cout << i << "\t" << sum/kmer_abundance_stepsize << endl;
				sum = 0;
			} else {
				sum += array->at(i).second;
				
			}
			
		}
		
		cout << "last" << "\t" << sum/((array->size())%kmer_abundance_stepsize) << "\t(last) "<< endl;
		cout << "totsum: " << totsum << endl;
	
	}
	#endif
	
	// -------------------------------------------------------------------------
	
	
	
		
	
	
	
	
	HashTable * cluster_kmer_hash;
	
	
	
	// --------- CLUSTERS ----------- // cluster objects would have been nice, but I am afraid of memory inefficiency
	

	HashedArrayTree_rwlock<mystring * > * cluster_consensus_sequences; //with lock
	
	
	
	HashedArrayTree_rwlock<cluster_member_list * > * cluster_member_lists; //with lock
	
	
	HashedArrayTree_rwlock<omp_lock_t > * cluster_locks; //with lock
	
	
	int last_cluster=-1;

	
	
	
	
	
	
	// -------------------------------------------------------------------------
	
	//iterate through all reads to create clusters and match against existing clusters

	
	
	size_t sequence_count = inputSequences->size();
	
	if (total_read_count != sequence_count) {
		cerr << "total_read_count != sequence_count" << endl;
		exit(1);
	}
	
	if (sequence_count == 0 ) {
		cerr << "error: no sequences found..." << endl;
		exit(1);
	}
	
	//int num_threads = thread_count;
	const int top_n_clusters = 3;
	
	
	int stat_real_cluster_count=0;
	int reads_processed = 0;
	
	const int hat_increase_steps = 1048576;
	
	
	if (thread_count_requested > 0) {
		omp_set_num_threads(thread_count_requested);
	}

	#ifdef DEBUG_DEADLOCK
	for (int i = 0; i < max_thread_count; ++i)
	{
		thread_status[i] =0;
		reads_per_thread[i] = 0;
	}
	#endif
	
	size_t partial_read_count[max_thread_count];
	for (int i = 0; i < max_thread_count; ++i)
	{
		partial_read_count[i] =0;
	}
	
	size_t chunk_start = 0;
	int chunk_size = 20000;
	
	//#####################################################################################
	
	#pragma omp parallel 
	{
	
		int this_thread_id = omp_get_thread_num();
		
		if (this_thread_id == 0) {
			thread_count = omp_get_num_threads();
			
			
			
			int size_t_size_bits = sizeof(size_t) * CHAR_BIT; // usually CHAR_BIT = 8
			
			if (thread_count > size_t_size_bits) {
				cerr << "error: number of threads bigger than number of bits in size_t type." << sizeof(size_t) << endl;
				cerr << "threads: " << thread_count << endl;
				cerr << "bit in size_t: " << size_t_size_bits << endl;
				exit(1);
				
			}
			
			
#ifdef DEBUG
			try {
#endif
				cluster_consensus_sequences = new HashedArrayTree_rwlock<mystring * >(true, thread_count, (size_t) 20, NULL);
#ifdef DEBUG
			} catch (bad_alloc& ba) {
				cerr << "error: (cluster_consensus_sequences) bad_alloc caught: " << ba.what() << endl;
				exit(1);
			}
#endif
			
#ifdef DEBUG
			cluster_consensus_sequences->name = string("cluster_consensus_sequences");
#endif
			
			
#ifdef DEBUG
			try {
#endif
				cluster_member_lists = new HashedArrayTree_rwlock<cluster_member_list * >(true, thread_count, (size_t)20, NULL);
#ifdef DEBUG
			} catch (bad_alloc& ba) {
				cerr << "error: (cluster_member_lists) bad_alloc caught: " << ba.what() << endl;
				exit(1);
			}
#endif
#ifdef DEBUG
			cluster_member_lists->name = string("cluster_member_lists");
#endif
			
			
			
#ifdef DEBUG
			try {
#endif
				cluster_locks = new HashedArrayTree_rwlock<omp_lock_t >(true, thread_count, (size_t) 20);
#ifdef DEBUG
			} catch (bad_alloc& ba) {
				cerr << "error: (cluster_member_lists) bad_alloc caught: " << ba.what() << endl;
				exit(1);
			}
#endif
#ifdef DEBUG
			cluster_locks->name = string("cluster_member_lists");
#endif
			
			cluster_consensus_sequences->reserve(2*hat_increase_steps);
			cluster_member_lists->reserve(2*hat_increase_steps);
			cluster_locks->reserve(2*hat_increase_steps);
			
			
			cluster_kmer_hash = new HashTable();
			
			log_stream << "number of threads actually used: " << thread_count<< endl;

			
			log_stream << "#reads\t#clusters";
			log_stream << "\tseconds";
			log_stream << "\tVM[kb]\tRSS[kb]";
			log_stream << endl;
			
			process_mem_usage(vm, rss);
			log_stream << "0\t0\t0\t" << (int)vm << "\t" << (int)rss << endl;
			
			#pragma omp flush
			
		}
		//int num_threads = omp_get_num_threads();
		//thread_count = num_threads;
		
		
//#pragma omp critical(cerr)
//		{
//		cerr << "real thread_count: " << omp_get_num_threads() << endl;
//		cerr << "thread starts: " << this_thread_id << endl;
//		}
		
		// barrier needed for shared data structures
		#pragma omp barrier
		
		#pragma omp flush
		
		HashedArrayTree<short> * countOfOverlappingKmers = 0; // thread private
		HashedArrayTree<int > * lastreadseen=0;  // thread private

		
		vector<short> * vector_of_offsets;
		vector<char> * alignment_column;
		vector<bool> * aminoacid_occurence;
		vector<int> * aminoacid_scores;
		
		int max_cluster_id_array[top_n_clusters];
		int max_cluster_kmer_count_array[top_n_clusters];
		short max_cluster_offsets[top_n_clusters];
		
		vector_of_offsets = new vector<short>();
		vector_of_offsets->resize(max_protein_length);
		
		alignment_column = new vector<char>();
		alignment_column->resize(1000000);
		
		aminoacid_occurence = new vector<bool>(aminoacid_count, false);
		aminoacid_scores = new vector<int>(aminoacid_count, 0);
		
		
		mystring * sequence = 0;
		
		#ifndef USE_HashTableSimple
		HashTable::iterator cluster_kmer_hash_it;
		#endif
		
		int** cluster_aminoacid_counts = new int*[max_protein_length];
		
		for (int i = 0 ; i <max_protein_length ; ++i ) {
			cluster_aminoacid_counts[i]=new int[aminoacid_count];
		}
				
#ifdef DEBUG
		try {
#endif
			countOfOverlappingKmers = new HashedArrayTree<short>(true, 20, 0); 
#ifdef DEBUG
		} catch (bad_alloc& ba) {
			cerr << "error: (countOfOverlappingKmers) bad_alloc caught: " << ba.what() << endl;
			exit(1);
		}
#endif
		
#ifdef DEBUG
		countOfOverlappingKmers->name = string("countOfOverlappingKmers");
#endif
		
		countOfOverlappingKmers->reserve(5*hat_increase_steps);
		//if ( countOfOverlappingKmers->capacity() < 5*hat_increase_steps ) {
		//	cerr << "ERROR!!!!!!!!!!" << endl;
		//	exit(1);
		//}
		
#ifdef DEBUG
		try {
#endif
			lastreadseen = new HashedArrayTree<int >(true, 20, -1); // this avoids initialization of match counts
#ifdef DEBUG
		} catch (bad_alloc& ba) {
			cerr << "error: (lastreadseen) bad_alloc caught: " << ba.what() << endl;
			exit(1);
		}
#endif
		
#ifdef DEBUG
		lastreadseen->name = string("lastreadseen");
#endif
		
		
		lastreadseen->reserve(5*hat_increase_steps);
		
		size_t first_read_id;
		size_t last_read_id;
		
		while (1) {
		
			bool all_reads_processed =false;
			#pragma omp critical(chunkstart)
			{
				if (chunk_start >= sequence_count) {
					all_reads_processed = true;
				} else {
					first_read_id = chunk_start;
					
					last_read_id = chunk_start+chunk_size-1;
					if (sequence_count - 1 < last_read_id) {
						last_read_id = sequence_count - 1;
					}
						
					chunk_start += chunk_size;
				}
			}
			if (all_reads_processed) {
				break;
			}
			
			
					
		for (size_t read_id = first_read_id ; read_id <= last_read_id; read_id++) { // ############################# READ LOOP ###########
			MACRO_THREAD_UPDATE(this_thread_id, __LINE__)
			
			partial_read_count[this_thread_id] = read_id-first_read_id; // no flush to avoid performance-decrease (there are enough flushes elsewhere.. ;) )
			
			//#pragma omp critical(reads_per_thread)
			//{
			//reads_per_thread[this_thread_id]++;  // too slow!
			//}
			
						
			bool repeat_loop = false;
		
						
						
			
			if (read_id % 1000 == 0) { // TODO not perfect, but I do not know better !
				//cerr << "XXt: " << this_thread_id << " new_cluster: " << new_cluster << " cap: " << countOfOverlappingKmers->capacity() << endl;
				#pragma omp flush(last_cluster)
				int reserve_max = last_cluster+5*hat_increase_steps;
				
				countOfOverlappingKmers->reserve(reserve_max); //privat
				
				lastreadseen->reserve(reserve_max); //privat
				
			}
#pragma omp flush(last_cluster)
			if (last_cluster > (int)countOfOverlappingKmers->capacity()-1000) {
				cerr << "last_cluster > (int)countOfOverlappingKmers->capacity()-1000" << endl;
				cerr << "last_cluster: " << last_cluster << endl;
				cerr << "countOfOverlappingKmers->capacity(): " << countOfOverlappingKmers->capacity() << endl;
				exit(1);
			}
			
			
			//cout << "read_id: " << read_id<< endl;
		
			sequence = &((*inputSequences)[read_id].second);
			KmerIterator * mykmerit;
		#ifdef DEBUG		
			try {
		#endif	
				mykmerit = new KmerIterator(sequence->c_str(), sequence->length(), 0, kmerlength, aminoacid_int2ASCII, aminoacid_ASCII2int, aminoacid_count);
		#ifdef DEBUG
			} catch (bad_alloc& ba) {
				cerr << "error: (new KmerIterator) bad_alloc caught: " << ba.what() << endl;
				exit(1);
			}
		#endif	
			//int max_cluster = -1;
			//int max_cluster_kmer_count = 0;
			
			
			
			for (int i = 0; i < top_n_clusters; ++i) {
				max_cluster_id_array[i]=-1;
				max_cluster_kmer_count_array[i]=0;
				
			}
			
			//cout << "------------- "<< read_id << " " << *(inputSequences->at(0).second) << endl;

			
			// check if read has matches to known clusters, iterate through read kmers
			#ifdef TIME
			#pragma omp master
			{
			clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &kmersearch_start);
			}
			#endif
			kmer_appearance_list::iterator mylist_it;
			
			kmer_appearance_list * mylist;
			// try to get a read lock, otherwise wait..
			#ifndef USE_HashTableSimple
			cluster_kmer_hash->set_reader_lock();
			#endif
			
			while (mykmerit->nextKmer()) {
				int code = mykmerit->code;
				//int pos = mykmerit->kmer_start_pos;
				
				//cout << "here code: " << code<< endl;
				
				#ifdef USE_HashTableSimple
				if (true) {
					mylist = (*cluster_kmer_hash)[code];
				#else
				cluster_kmer_hash_it = cluster_kmer_hash->find(code);
				
				// check if k-mer is indexed:
				if (cluster_kmer_hash_it != cluster_kmer_hash->end()) {
					mylist = cluster_kmer_hash_it->second;
				#endif
					
					//cout << "---- "<< endl;
					//mylist->resetIterator();
					// check all occurrences of that k-mer
					mylist->set_reader_lock();
					for (mylist_it = mylist->begin(); mylist_it != mylist->end(); ++mylist_it) {
					//while (mylist->nextElement()) {
						int clusterhit = mylist_it.getFirst();
						
						
						#ifdef DEBUG
						#pragma omp flush(last_cluster)
						if (clusterhit > last_cluster) {
							cerr << "clusterhit > last_cluster" << endl;
							cerr << "clusterhit: " << clusterhit << endl;
							cerr << "last_cluster: " << last_cluster << endl;
							exit(1);
						}
						
						if (clusterhit >= (int)lastreadseen->capacity()) {
							cerr << "clusterhit > lastreadseen->capacity()" << endl;
							cerr << "clusterhit: " << clusterhit << endl;
							cerr << "lastreadseen->capacity(): " << lastreadseen->capacity() << endl;
							exit(1);
						}
						#endif
						
						
						
						//cerr << "clusterhit: " << clusterhit << endl;
						if ((*lastreadseen)[clusterhit] == (int) read_id) {
						//int blubbbla = lastreadseen->at(clusterhit);
						//int blubbbla2 = (int) read_id;
						//if (blubbbla == blubbbla2) { //  back []
							#ifdef DEBUG
							if (clusterhit > (int) countOfOverlappingKmers->capacity()) {
								cerr << "clusterhit > countOfOverlappingKmers->capacity()" << endl;
								cerr << "clusterhit: " << clusterhit << endl;
								cerr << "countOfOverlappingKmers->capacity(): " << countOfOverlappingKmers->capacity() << endl;
								exit(1);
							}
							#endif						
							short& cluster_overlap_count = (*countOfOverlappingKmers)[clusterhit];
							cluster_overlap_count++;
							
						
							
							//(*countOfOverlappingKmers)[clusterhit]++;
							//cout << "seen before" << endl;
							
							// check for threshold
							if (cluster_overlap_count >= cluster_kmer_overlap_threshold) {
								
								//int search_pos_in_top_list = top_n_clusters-1;
								
								//cout << "----" << endl;
								//for (int i = 0; i < top_n_clusters; ++i) {
									
								//	cout << i << ") " << max_cluster_kmer_count_array[i] << " " << max_cluster_id_array[i] << endl;
								//}
								//exit(0);
								
								// check if count is better than the n-th best cluster seen so far 
								// and check if it is already there
								
								// will find first cluster of same size, or first cluster that is bigger, or top of list
								
								bool cluster_already_in_top_list=false;
								//bool best_cluster= false;
								//bool insert_cluster = false;
								
								int search_pos_in_top_list = 0;
								
								
								// search for same size cluster:
								for ( search_pos_in_top_list = 0;  search_pos_in_top_list < top_n_clusters; ++search_pos_in_top_list) {
									if (max_cluster_id_array[search_pos_in_top_list] == clusterhit) {
										cluster_already_in_top_list = true;
										break;
									}
									
								}
								
								if (cluster_already_in_top_list) {
									// update value:
									
									max_cluster_kmer_count_array[search_pos_in_top_list] = cluster_overlap_count;
									
									// try to move cluster up in list
									while (true) {
										if (search_pos_in_top_list == 0) {
											break;
										}
										
										if (max_cluster_kmer_count_array[search_pos_in_top_list] > max_cluster_kmer_count_array[search_pos_in_top_list-1] ) {
											swap(max_cluster_kmer_count_array[search_pos_in_top_list], max_cluster_kmer_count_array[search_pos_in_top_list-1]);
											swap(max_cluster_id_array[search_pos_in_top_list], max_cluster_id_array[search_pos_in_top_list-1]);
											
											search_pos_in_top_list--;
											
										} else {
											break;
										}
										
									} 
									
									
								} 
															
								
								//cout << "search_pos_in_top_list: " << search_pos_in_top_list << endl;
								//if (search_pos_in_top_list >= top_n_clusters) exit(1);
								
								if (! cluster_already_in_top_list && (cluster_overlap_count > max_cluster_kmer_count_array[top_n_clusters-1])) {
									
									
									
									// insert cluster at insertion point, if inspoint is still in list
									
									
									search_pos_in_top_list = top_n_clusters-1;
									
									max_cluster_kmer_count_array[top_n_clusters-1] = cluster_overlap_count;
									max_cluster_id_array[top_n_clusters-1] = clusterhit;
									
									//cout << "AA " << max_cluster_kmer_count_array[search_pos_in_top_list]<< " " << max_cluster_kmer_count_array[search_pos_in_top_list-1] << endl;
									
									while (max_cluster_kmer_count_array[search_pos_in_top_list] > max_cluster_kmer_count_array[search_pos_in_top_list-1] ) {
										//cout << "AA--" << max_cluster_kmer_count_array[search_pos_in_top_list]<< " " << max_cluster_kmer_count_array[search_pos_in_top_list-1] << endl;
										swap(max_cluster_kmer_count_array[search_pos_in_top_list], max_cluster_kmer_count_array[search_pos_in_top_list-1]);
										swap(max_cluster_id_array[search_pos_in_top_list], max_cluster_id_array[search_pos_in_top_list-1]);
										search_pos_in_top_list--;
										if (search_pos_in_top_list == 0) {
											break;
										}
										//exit(0);
									}
										
																	
								} 
								
								
																
								//cout << "----" << endl;
								//for (int i = 0; i < top_n_clusters; ++i) {
									
								//	cout << i << ") " << max_cluster_kmer_count_array[i] << " " << max_cluster_id_array[i] << endl;
								//}
								#ifdef DEBUG
								if (max_cluster_id_array[0] == max_cluster_id_array[1]) {
									cerr << "error: max_cluster_id_array[0] == max_cluster_id_array[1]" << endl;
									exit(1);	
								} 
								if (max_cluster_id_array[0] == max_cluster_id_array[2]) {
									cerr << "error: ax_cluster_id_array[0] == max_cluster_id_array[2]" << endl;
									exit(1);	
								} 
								#endif
								//if (max_cluster_kmer_count_array[0] > 100000) exit(0);  	
								
								
							}
							
							
						} else {
							#ifdef DEBUG
							if (clusterhit >= (int) lastreadseen->capacity() ){
								cerr << "clusterhit >= lastreadseen->capacity()" << endl;
								cerr << clusterhit << " " << lastreadseen->capacity() << endl;
								exit(1);
							}
							
							if (clusterhit >= (int) countOfOverlappingKmers->capacity() ){
								cerr << "clusterhit >= countOfOverlappingKmers->capacity()" << endl;
								cerr << clusterhit << " " << countOfOverlappingKmers->capacity() << endl;
								exit(1);
							}
							
							#endif
							
							(*lastreadseen)[clusterhit] = read_id;
							(*countOfOverlappingKmers)[clusterhit]=1;
							//cout << "not seen before" << endl;
						}// end if
					} // end for
					mylist->unset_reader_lock();
					
				} // end if
				
				
			} // end while ( kmer iteration )
			#ifndef USE_HashTableSimple
			cluster_kmer_hash->unset_reader_lock();
			#endif
			
			MACRO_THREAD_UPDATE(this_thread_id, __LINE__)
				
			#ifdef TIME
			#pragma omp master
			{
			clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &kmersearch_end);
			kmersearch_total = add_time( overlap_total, diff(kmersearch_start, kmersearch_end));
			}
			#endif
			
			
			delete mykmerit; 	
					
	//		if (read_id > 1791) {
	//			cout << "read_id: " << read_id<< endl;
	//			kmer_appearance_list * mylist;
	//			cluster_kmer_hash_it = cluster_kmer_hash.find(2644139);
	//			if (cluster_kmer_hash_it != cluster_kmer_hash.end()) {
	//				mylist = cluster_kmer_hash_it->second;
	//				if (mylist->start == NULL ) {
	//					cerr << "NULL" << endl;
	//					exit(1);
	//				}
	//				mylist->print();
	//				//if (mylist->start->data_array1[0] != 1601) {
	//				//	cout << "read_id" << read_id << endl; 
	//				//	exit(0);
	//				//}
	//			} else {
	//				cout << "not found: " << read_id << endl;
	//				exit(0);
	//			}				
	//		}
			
			//cout << max_cluster_id_array[0] << endl;
			
			#ifdef TIME
			#pragma omp master
			{
			clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &validation_start);
			}
			#endif
			int validated_overlaps=0;
			
			// validate overlaps
			if (max_cluster_id_array[0] >= 0) { // number of k-mer-hits indicates a hit, needs to be verified.
				
				// (probably) found cluster for the read !
				
				//int seed_sequence_number =  (*cluster_seedread)[max_cluster];
				
				if (false) {
					cout << "yeah --------------------" << endl;
					cout << "read: " << read_id << endl;
					cout << "readseq: " << sequence << endl;
					//cout << "max_cluster: " << max_cluster << endl;
					//cout << "max_cluster_kmer_count: " << max_cluster_kmer_count << endl;
					
					//cout << "cluster seed num: " << seed_sequence_number << endl;
				
					//cout << "cluster seed seq: " << *((*inputSequences)[seed_sequence_number].second) << endl;
				}
				
				
				
				
				// check for each cluster in top-list if the overlap is real:
				
				for (int cluster_it=0; cluster_it < top_n_clusters; ++cluster_it) {
					
					//if (read_id == 414) {
					//	for (int i = 0; i < top_n_clusters; ++i) {
							
					//		cout << i << "] " << max_cluster_kmer_count_array[i] << " " << max_cluster_id_array[i] << endl;
					////	}
					//}
					
					MACRO_THREAD_UPDATE(this_thread_id, __LINE__)
					if (max_cluster_id_array[cluster_it] == -1) {
						break;
					}
					
					// store overlaps in vector_of_offsets
					#ifndef USE_HashTableSimple
					cluster_kmer_hash->set_reader_lock();
					#endif
					int number_of_overlapping_kmers_seen = getOffsets(sequence->c_str(),
																	  sequence->length(),
																	  vector_of_offsets,	// <-- output
																	  cluster_kmer_hash,
																	  max_cluster_id_array[cluster_it]);
					#ifndef USE_HashTableSimple
					cluster_kmer_hash->unset_reader_lock();
					#endif
					
					// majority vote on offsets				
					triplet<bool, short, int> major = majority_vote<short>(vector_of_offsets, number_of_overlapping_kmers_seen);
					bool overlap_found=major.first;
					short majority_offset=major.second;
					
					// check if the number on of k-mers at main diagonal is still sufficient for threshold
					if (false) {
						if (overlap_found) {
							int real_kmer_count = 0;
							//cout << "--" << endl;
							for (int i = 0; i<number_of_overlapping_kmers_seen; ++i) {
								//cout <<(*vector_of_offsets)[i] << endl;
								if ((*vector_of_offsets)[i] == majority_offset) {
									real_kmer_count++;	
									
								}
							}
							
							//cout << "real_kmer_count: " << real_kmer_count << endl;
							if (real_kmer_count < cluster_kmer_overlap_threshold) {
								overlap_found = false;
							}
						}
					}
					
					int cluster;
					
					MACRO_THREAD_UPDATE(this_thread_id, __LINE__)
					while (1) { // not a real loop, just want to be able to leave block.. ;)
						if (not overlap_found) {
							break;
						}
			
						//check overlap:
						cluster = max_cluster_id_array[cluster_it];
						
						// lock cluster 
						omp_lock_t * my_cluster_lock = &(*cluster_locks)[cluster];
						//#pragma omp critical(cerr)
						//{
						//	cerr << "thread: " << this_thread_id << " clusterlock: " << cluster << endl;
						//}
						MACRO_THREAD_UPDATE(this_thread_id, __LINE__)
						ScopedLock lalala(my_cluster_lock);
						//#pragma omp flush
						MACRO_THREAD_UPDATE(this_thread_id, __LINE__)
						//check if it exists
						cluster_member_lists->set_reader_lock();
						
						if ((*cluster_member_lists)[cluster] == NULL ) {
							// cluster does not exist anymore!
							overlap_found = false;
							cluster_member_lists->unset_reader_lock();
							break;
						}
						cluster_member_lists->unset_reader_lock();
						MACRO_THREAD_UPDATE(this_thread_id, __LINE__)
						
						cluster_consensus_sequences->set_reader_lock();
						mystring * consensus_seq;
						#ifdef DEBUG
						try {
						consensus_seq = cluster_consensus_sequences->at(cluster);
						} catch (out_of_range& oor) {
							cerr << "Out of Range error:(consensus_seq = cluster_consensus_sequences->at(cluster)) " << oor.what() << endl;
							exit(1);
						}
						#else
						consensus_seq = (*cluster_consensus_sequences)[cluster];
						#endif
						cluster_consensus_sequences->unset_reader_lock();
						
						if (consensus_seq == NULL) {
							
							#ifdef DEBUG
							cluster_member_list * cluster_memberlist;
							try {
							cluster_memberlist = cluster_member_lists->at(cluster);
							} catch (out_of_range& oor) {
								cerr << "Out of Range error:(cluster_member_list * cluster_memberlist = cluster_member_lists->at(cluster)) " << oor.what() << endl;
								exit(1);
							}
							#else
							cluster_member_list * cluster_memberlist = (*cluster_member_lists)[cluster];
							#endif
														
							
							int cluster_read = cluster_memberlist->start->data_array1[0];
							
							#ifdef DEBUG
							//pair<mystring,mystring> pp = inputSequences->at(cluster_read);
							//string tempstring = string(pp.second.c_str());
							
							
							try {
							//consensus_seq = new mystring(tempstring);
								
							consensus_seq = &(inputSequences->at(cluster_read).second);
							} catch (out_of_range& oor) {
								cerr << "Out of Range error:(consensus_seq = inputSequences->at(cluster_read).second) " << oor.what() << endl;
								exit(1);
							}
							#else
							//consensus_seq = new mystring(string((*inputSequences)[cluster_read].second.c_str()));
							
							consensus_seq = &((*inputSequences)[cluster_read].second);
							#endif
						}
						
						//mystring help_sequence = mystring(sequence->c_str());
						overlap_found = computeSequenceOverlap(majority_offset, consensus_seq, sequence, score_matrix, major.third);
						MACRO_THREAD_UPDATE(this_thread_id, __LINE__)
					
						//exit(1);
						//computeSequenceOverlap(int offset, string * a, string * b, short * score_matrix)
					
					
						if (not overlap_found) {
							break;
						}
						MACRO_THREAD_UPDATE(this_thread_id, __LINE__)
						validated_overlaps++;
						//if (read_id == 414) {
						//	cout<< "::: " << cluster_it  << " " << validated_overlaps << endl;
						//}
						if (cluster_it != (validated_overlaps-1)) {
							max_cluster_id_array[validated_overlaps-1] = max_cluster_id_array[cluster_it];
							max_cluster_kmer_count_array[validated_overlaps-1] = max_cluster_kmer_count_array[cluster_it];
						}
						max_cluster_offsets[validated_overlaps-1] = majority_offset;
						break;
						
					}// end while(1)
					
					if (not overlap_found) {
						max_cluster_id_array[cluster_it] = -1;
						max_cluster_kmer_count_array[cluster_it] = 0;
					}
					
				}
								
				
					
				MACRO_THREAD_UPDATE(this_thread_id, __LINE__)
				
				
				
			} // end if (max_cluster_id_array[0] >= 0)
			
			#ifdef TIME
			#pragma omp master
			{
			clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &validation_end);
			validation_total = add_time( validation_total, diff(validation_start, validation_end));
			//cerr << validation_start << " " << validation_end << endl;
			}
			#endif
			
			
			MACRO_THREAD_UPDATE(this_thread_id, __LINE__)
				
			// ************ process overlaps ************
			//cout << "validated_overlaps"  << endl;
			#ifdef TIME
			#pragma omp master
			{
			clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &overlap_start);
			}
			#endif
				
				// disable merging
				//if (validated_overlaps >=2 ) {
				//	validated_overlaps=1 ;
				//}
				
				
			while (true) {
				if (validated_overlaps == 0) {
					MACRO_THREAD_UPDATE(this_thread_id, __LINE__)
					// did not found matches, insert as new cluster
					#pragma omp atomic
					stat_real_cluster_count++;
					
					omp_lock_t * my_cluster_lock;
					
					int new_cluster;
					
					
					MACRO_THREAD_UPDATE(this_thread_id, __LINE__)
					
					
					#pragma omp critical(last_cluster_update)
					{
						last_cluster++;
						
						new_cluster = last_cluster;
						// make sure the HAT-arrays have enough memory allocated
						// this has the advantage I can use []-operator instead of at()-function
						
						// increase global arrays (only one thread can do this!)
						if ((new_cluster % hat_increase_steps) == 0) { // make sure that other threads have sufficient space!
							int reserve_max = new_cluster+hat_increase_steps+1;
MACRO_THREAD_UPDATE(this_thread_id, __LINE__)
							
							
							cluster_consensus_sequences->set_writer_lock();
							cluster_consensus_sequences->reserve(reserve_max); //shared!
							cluster_consensus_sequences->unset_writer_lock();
MACRO_THREAD_UPDATE(this_thread_id, __LINE__)
							cluster_member_lists->set_writer_lock();
							cluster_member_lists->reserve(reserve_max); //shared!
							cluster_member_lists->unset_writer_lock();
MACRO_THREAD_UPDATE(this_thread_id, __LINE__)
							cluster_locks->set_writer_lock();
							cluster_locks->reserve(reserve_max); //shared!
							cluster_locks->unset_writer_lock();
						}
//MACRO_THREAD_UPDATE(this_thread_id, __LINE__)
						my_cluster_lock = &(*cluster_locks)[new_cluster];
						
MACRO_THREAD_UPDATE(this_thread_id, __LINE__)
						//#pragma omp critical(cerr)
						//{
						//	cerr << "thread: " << this_thread_id << " clusterlock: " << new_cluster << endl;
						//}
						
						
						//#pragma omp flush
						
						
						
						
						
						//omp_set_lock(my_cluster_lock);
					} // end pragma critical
					
						
					
					
					
					omp_init_lock(my_cluster_lock);
					ScopedLock lalala(my_cluster_lock);
					
					
					
MACRO_THREAD_UPDATE(this_thread_id, __LINE__)
					// insert seed as first member:
					
					cluster_member_lists->set_reader_lock();
					#ifdef DEBUG
					
					cluster_member_list * templist;
					
					try {
						templist = cluster_member_lists->at(new_cluster);
					} catch (out_of_range& oor) {
						cerr << "Out of Range error:(templist = cluster_member_lists->at(new_cluster)) " << oor.what() << endl;
						exit(1);
					}
					
					
					cluster_member_list*& mymemberlist = cluster_member_lists->at(new_cluster);
					#else
					cluster_member_list*& mymemberlist = (*cluster_member_lists)[new_cluster]; // alias to a pointer
					#endif
					cluster_member_lists->unset_reader_lock();
MACRO_THREAD_UPDATE(this_thread_id, __LINE__)
					//cout << mymemberlist << endl;
					assert( (mymemberlist == NULL) );
					
					cluster_member_lists->set_writer_lock();
//MACRO_THREAD_UPDATE(this_thread_id, __LINE__)
					#ifdef DEBUG
					try {
					#endif
						mymemberlist = new cluster_member_list(); //modified cluster_member_lists ! (alias to a pointer)
					#ifdef DEBUG
					} catch (bad_alloc& ba) {
						cerr << "error: (new cluster_member_list) bad_alloc caught: " << ba.what() << endl;
						exit(1);
					}
					#endif	
					cluster_member_lists->unset_writer_lock();
					mymemberlist->append(read_id, 0); // cluster lock is enough
					MACRO_THREAD_UPDATE(this_thread_id, __LINE__)
					
					
					addConsensusSequence(sequence, new_cluster, cluster_kmer_hash);
					
					#ifdef DEBUG
					// see if new sequence has been added correctly
					if ( omp_test_lock(my_cluster_lock) ) {
						cerr << "omp_test_lock: 1" << endl;
						exit(1);
					}
					
					testConsensusSequence(sequence, new_cluster, cluster_kmer_hash, 7);
					#endif
					
					MACRO_THREAD_UPDATE(this_thread_id, __LINE__)
					//omp_unset_lock(my_cluster_lock);
					
					//#pragma omp flush
					
					break; // end of cluster creation
					
				} else if (validated_overlaps == 1) { // read has match to exactly one cluster
					MACRO_THREAD_UPDATE(this_thread_id, __LINE__)
					int clusterhit = max_cluster_id_array[0];
					
					// how can I be sure cluster still exists ?
					// if  cluster does not exist anymore, there will be no wrong cluster instead, just empty.
					
					// lock the thing (lock exists!), check if it exists, if not fall back..
					
					omp_lock_t * my_cluster_lock = &(*cluster_locks)[clusterhit];
//#pragma omp critical(cerr)
					//{
					//cerr << "thread: " << this_thread_id << " clusterlock: " << clusterhit << endl;
					//}
					
					ScopedLock lalala(my_cluster_lock);
					
					#ifdef DEBUG
					if ( omp_test_lock(my_cluster_lock) ) {
						cerr << "omp_test_lock: 5" << endl;
						exit(1);
					}
					#endif
					//#pragma omp flush
					
					//omp_set_lock(my_cluster_lock);
					MACRO_THREAD_UPDATE(this_thread_id, __LINE__)
					//check if cluster exists...
					
					cluster_member_lists->set_reader_lock();
					
					if ((*cluster_member_lists)[clusterhit] == NULL ) {
						// cluster does not exist anymore!
						validated_overlaps = 0; // fall-back.
						repeat_loop = true;
						cluster_member_lists->unset_reader_lock();
						//omp_unset_lock(my_cluster_lock);
						break; // break the while(1) loop
					}
					
					// get member list of cluster to add new member
					#ifdef DEBUG
					
					cluster_member_list * templist; // already locked
					try {
						templist = cluster_member_lists->at(clusterhit);
					} catch (out_of_range& oor) {
						cerr << "Out of Range error:(templist = cluster_member_lists->at(clusterhit)) " << oor.what() << endl;
						exit(1);
					}
					MACRO_THREAD_UPDATE(this_thread_id, __LINE__)
					
					cluster_member_list*& mymemberlist = cluster_member_lists->at(clusterhit);
					#else
					cluster_member_list*& mymemberlist = (*cluster_member_lists)[clusterhit]; // alias to a pointer
					#endif	
					cluster_member_lists->unset_reader_lock();
					
					if (mymemberlist == NULL){ // in case cluster had size one, I think
						#ifdef DEBUG
						try {
						#endif
							cluster_member_lists->set_writer_lock();
							mymemberlist = new cluster_member_list();
							cluster_member_lists->unset_writer_lock();
						#ifdef DEBUG
						} catch (bad_alloc& ba) {
							cerr << "error: (new cluster_member_list) bad_alloc caught: " << ba.what() << endl;
							exit(1);
						}
						#endif
					}
					MACRO_THREAD_UPDATE(this_thread_id, __LINE__)
					
					
					
					// find old consensus sequence (read or real consensus)
					cluster_consensus_sequences->set_reader_lock();
					#ifdef DEBUG
					mystring * consensus_old;
					try {
					consensus_old = cluster_consensus_sequences->at(clusterhit);
					} catch (out_of_range& oor) {
						cerr << "Out of Range error:(consensus_old = cluster_consensus_sequences->at(clusterhit)) " << oor.what() << endl;
						exit(1);
					}
					#else
					mystring * consensus_old = (*cluster_consensus_sequences)[clusterhit];
					#endif	
					cluster_consensus_sequences->unset_reader_lock();
					
					
					
					mymemberlist->append(read_id, max_cluster_offsets[0]); // locked by cluster lock
					
		
					MACRO_THREAD_UPDATE(this_thread_id, __LINE__)
					
					if (consensus_old == NULL) {
						
						// consensus seuquence does not exist yet, seed member had been used
						
						int seed_id = mymemberlist->start->data_array1[0];
						
						#ifdef DEBUG
					//ugga = 9;
						try {
						consensus_old = &(inputSequences->at(seed_id).second);
						} catch (out_of_range& oor) {
							cerr << "Out of Range error:(consensus_old = inputSequences->at(seed_id).second) " << oor.what() << endl;
							exit(1);
						}
						#else
						consensus_old = &((*inputSequences)[seed_id].second);
						#endif
						
						//#ifdef DEBUG
						//if (consensus_old->data == NULL) {
						//	cerr << "2consensus_old->data == NULL" << endl;
						//	exit(1);
						//}
						//#endif
					} else {
						//#ifdef DEBUG
						//if (consensus_old->data == NULL) {
						//	cerr << "1consensus_old->data == NULL" << endl;
						//	exit(1);
						//}
						//#endif
					}
					
#ifdef DEBUG
					// see if old sequence had been added correctly
					if ( omp_test_lock(my_cluster_lock) ) {
						cerr << "omp_test_lock: 2" << endl;
						exit(1);
					}
				//	testConsensusSequence(consensus_old, clusterhit, cluster_kmer_hash, ugga);
#endif
				
					
					// check if new member is a perfect full-length overlap
					bool extends_consensus = true;
					
					
					
					int member_len = sequence->length();
								
					//bool perfect_overlap = false;
					if (max_cluster_offsets[0] >= 0 && max_cluster_offsets[0]+member_len <= (int) consensus_old->length() ) { // check for full-length overlap
						extends_consensus = false;
					} 
					
					MACRO_THREAD_UPDATE(this_thread_id, __LINE__)
					int membercount = mymemberlist->getLength();
					//cout << "membercount: " << membercount << endl;
					
					
					if (clusterNeedsUpdate(membercount) || extends_consensus) { // mymemberlist->getLength() >= 100
					//if (true) {
						
						// ============ REMOVE OLD CONSENSUS K-MERS =============
						
						MACRO_THREAD_UPDATE(this_thread_id, __LINE__)
						// remove consensus sequence kmers for cluster "clusterhit"
						removeConsensusSequence(consensus_old, clusterhit, cluster_kmer_hash);
						
						MACRO_THREAD_UPDATE(this_thread_id, __LINE__)
						// delete old consensus, but not seed !
						
						cluster_consensus_sequences->set_writer_lock();
						#ifdef DEBUG
						mystring * tempstr;
						try {
						tempstr = cluster_consensus_sequences->at(clusterhit);
						} catch (out_of_range& oor) {
							cerr << "Out of Range error:(tempstr = cluster_consensus_sequences->at(clusterhit)) " << oor.what() << endl;
							exit(1);
						}
						if (tempstr != NULL) {
							consensus_old->delete_data();
							delete consensus_old;
							cluster_consensus_sequences->at(clusterhit) = NULL;
						}
						#else					
						if ((*cluster_consensus_sequences)[clusterhit] != NULL) {
							consensus_old->delete_data();
							delete consensus_old;
							(*cluster_consensus_sequences)[clusterhit] = NULL;
						}
						#endif
						cluster_consensus_sequences->unset_writer_lock();
						//consensus_old = NULL;
						MACRO_THREAD_UPDATE(this_thread_id, __LINE__)
						
						// ============= SORT, UPDATE OFFSETS AND COMPUTE CONSENSUS ==============
						// copy members from list to array (temporarly)
						
						mystring * new_consensus_seq = computeConsensus(mymemberlist,  // locked by cluster lock
																	  inputSequences, 
																	  alignment_column,
																	  aminoacid_occurence,
																	  aminoacid_scores,
																	  score_matrix,
																	  cluster_aminoacid_counts	);
						cluster_consensus_sequences->set_writer_lock();
						#ifdef DEBUG
						try {
						cluster_consensus_sequences->at(clusterhit) = new_consensus_seq;
						} catch (out_of_range& oor) {
							cerr << "Out of Range error:(cluster_consensus_sequences->at(clusterhit) = new_consensus_seq) " << oor.what() << endl;
							exit(1);
						}
						#else
						(*cluster_consensus_sequences)[clusterhit] = new_consensus_seq;
						#endif
						cluster_consensus_sequences->unset_writer_lock();

						
						
						MACRO_THREAD_UPDATE(this_thread_id, __LINE__)
						
						// insert kmers of new_consensus_seq into hash
						
						addConsensusSequence(new_consensus_seq, clusterhit, cluster_kmer_hash);
						
#ifdef DEBUG
						if ( omp_test_lock(my_cluster_lock) ) {
							cerr << "omp_test_lock: 3" << endl;
							exit(1);
						}
						// see if old sequence has been added correctly
						testConsensusSequence(new_consensus_seq, clusterhit, cluster_kmer_hash, 50);
#endif
						
					}
					
					//omp_unset_lock(my_cluster_lock);
					//#pragma omp flush
					
					
					break;
					
				} else { //if (validated_overlaps >= 2) {
					MACRO_THREAD_UPDATE(this_thread_id, __LINE__)
					//cout << "validated_overlaps >= 2 !!!!!!!!" << endl;
					
					
					//cout << "read_id: " << read_id << endl;
					//for (int i = 0; i < validated_overlaps; ++i) {
						
					//	cout << i << "))) " << max_cluster_kmer_count_array[i] << " " << max_cluster_id_array[i] << endl;
					//	cout << "offset: " << max_cluster_offsets[i] << endl;
					//}
					
					
					
					validated_overlaps=2; // sorry....
					
					
					int previous_master_cluster_index = 0;
					// merge clusters iteratively
					//for (int second_cluster_index = 1; second_cluster_index < validated_overlaps; ++second_cluster_index) {
					int second_cluster_index = 1;
						
					
					int cluster_1 = max_cluster_id_array[previous_master_cluster_index];
					int cluster_2 = max_cluster_id_array[second_cluster_index];
					
					
					// try to lock both clusters
					omp_lock_t * my_cluster_lock_1;
					omp_lock_t * my_cluster_lock_2;
					
					if (cluster_1 < cluster_2) { // this should help avoiding deadlock !
						my_cluster_lock_1 = &(*cluster_locks)[cluster_1];
						my_cluster_lock_2 = &(*cluster_locks)[cluster_2];
					} else {
						my_cluster_lock_1 = &(*cluster_locks)[cluster_2];
						my_cluster_lock_2 = &(*cluster_locks)[cluster_1];
					}
					
					
					//omp_set_lock(my_cluster_lock_1);
					//omp_set_lock(my_cluster_lock_2);
//#pragma omp critical(cerr)
//					{
//						cerr << "thread: " << this_thread_id << " clusterlock: " << cluster_1  << " " << cluster_2 << endl;
//					}
					
					MACRO_THREAD_UPDATE(this_thread_id, __LINE__)
					ScopedLock lalala1(my_cluster_lock_1);
					MACRO_THREAD_UPDATE(this_thread_id, __LINE__)
					ScopedLock lalala2(my_cluster_lock_2);
					//#pragma omp flush
					
					

					
					MACRO_THREAD_UPDATE(this_thread_id, __LINE__)
					// get their member lists:
					cluster_member_lists->set_reader_lock();
					cluster_member_list * cluster_1_memberlist = (*cluster_member_lists)[cluster_1];
					cluster_member_list * cluster_2_memberlist = (*cluster_member_lists)[cluster_2];
					

					// check if the clusters actually still exist
					if (cluster_1_memberlist == NULL ||  cluster_2_memberlist == NULL) {
						// cluster does not exist anymore!
						validated_overlaps = 1; // fall-back. 
						repeat_loop = true;
						cluster_member_lists->unset_reader_lock();
						//omp_unset_lock(my_cluster_lock_1);
						//omp_unset_lock(my_cluster_lock_2);
						break;
					}
					cluster_member_lists->unset_reader_lock();
					
					// get consensus sequences
					cluster_consensus_sequences->set_reader_lock();
					mystring * consensus_1 = (*cluster_consensus_sequences)[cluster_1];
					mystring * consensus_2 = (*cluster_consensus_sequences)[cluster_2];
					cluster_consensus_sequences->unset_reader_lock();
					
					bool consensus_1_isread = false;
					bool consensus_2_isread = false;
					
					if (consensus_1 == NULL) {
						consensus_1_isread = true;
						int cluster_read = cluster_1_memberlist->start->data_array1[0];
						
	//					#ifdef DEBUG
	//					consensus_1 = inputSequences->at(cluster_read).second;
	//					#else
	//					consensus_1 = (*inputSequences)[cluster_read].second;
	//					#endif
						consensus_1 = &(*inputSequences)[cluster_read].second;
					}
										
					if (consensus_2 == NULL) {
						consensus_2_isread = true;
						int cluster_read = cluster_2_memberlist->start->data_array1[0];
						
	//					#ifdef DEBUG
	//					consensus_2 = inputSequences->at(cluster_read).second;
	//					#else
	//					consensus_2 = (*inputSequences)[cluster_read].second;
	//					#endif
						consensus_2 = &(*inputSequences)[cluster_read].second;
					}

#ifdef DEBUG
					// see if old sequence has been added correctly
					testConsensusSequence(consensus_1, cluster_1, cluster_kmer_hash, 51);
					testConsensusSequence(consensus_2, cluster_2, cluster_kmer_hash, 52);
#endif
					
					//cout << cluster_1 << endl;
					//cout << cluster_2 << endl;
					
					//cout << max_cluster_offsets[previous_master_cluster_index] << endl;
					//cout << max_cluster_offsets[second_cluster_index] << endl;
					
					int offset_diff = max_cluster_offsets[previous_master_cluster_index] - max_cluster_offsets[second_cluster_index];
					//cout << "seq: " << *sequence << endl;
					//cout << "o: " << offset_diff << endl;
					
					//cout << consensus_1->c_str() << endl;
					//cout << consensus_2->c_str() << endl;
					bool has_overlap = computeSequenceOverlap(offset_diff, consensus_1, consensus_2, score_matrix, 0);
					
					//cout << (has_overlap?string("true"):string("false")) << endl;
					
					if ( has_overlap) {
					
						MACRO_THREAD_UPDATE(this_thread_id, __LINE__)
						//cout << cluster_1_memberlist->getLength() << endl;
						//cout << cluster_2_memberlist->getLength() << endl;
						
						//cout << *consensus_1 << endl;
						//cout << *consensus_2 << endl;
						//exit(0);
						
						int newread_offset;
						
						
						int slave_offset;
						
						int master_cluster_index;
						int slave_cluster_index;
						
						int master_cluster;
						int slave_cluster;
						
						cluster_member_list * master_cluster_member_list;
						cluster_member_list * slave_cluster_member_list;
						
						//string * master_consensus_sequence;
						//string * slave_consensus_sequence;
						//cout << "B" << endl;
						
						
						if (max_cluster_offsets[previous_master_cluster_index] >= max_cluster_offsets[second_cluster_index]) {
							// first sequence is left of second sequence
							
							master_cluster_index = previous_master_cluster_index;
							slave_cluster_index = second_cluster_index;
							
							cluster_member_lists->set_reader_lock();
							master_cluster_member_list = (*cluster_member_lists)[cluster_1];
							slave_cluster_member_list = (*cluster_member_lists)[cluster_2];
							cluster_member_lists->unset_reader_lock();
							
						} else {
							
							master_cluster_index =  second_cluster_index;
							slave_cluster_index = previous_master_cluster_index;
							
							cluster_member_lists->set_reader_lock();
							master_cluster_member_list = (*cluster_member_lists)[cluster_2];
							slave_cluster_member_list = (*cluster_member_lists)[cluster_1];
							cluster_member_lists->unset_reader_lock();
						}
						
						master_cluster = max_cluster_id_array[master_cluster_index];
						slave_cluster = max_cluster_id_array[slave_cluster_index];
						
						//master_consensus_sequence = consensus_2;
						//slave_consensus_sequence = consensus_1;
						
						slave_offset = max_cluster_offsets[master_cluster_index] - max_cluster_offsets[slave_cluster_index];
						newread_offset = max_cluster_offsets[master_cluster_index];
						
						//cout << "C" << endl;
						
						// update and copy slave members to master
						
						//slave_cluster_member_list->resetIterator();
						cluster_member_list::iterator mylist_it;
						
						//while (slave_cluster_member_list->nextElement()) {
						
						
						for (mylist_it = slave_cluster_member_list->begin(); mylist_it != slave_cluster_member_list->end(); ++mylist_it) {
							master_cluster_member_list->append(mylist_it.getFirst(), mylist_it.getSecond() + slave_offset); // protected by cluster lock
						}
						
						
						
						// append current read as member:
						master_cluster_member_list->append(read_id, newread_offset); // protected by cluster lock
						
						// delete stuff
						delete slave_cluster_member_list;
						cluster_member_lists->set_writer_lock();
						(*cluster_member_lists)[slave_cluster] = NULL;
						cluster_member_lists->unset_writer_lock();
						
						removeConsensusSequence(consensus_1, cluster_1, cluster_kmer_hash);
						removeConsensusSequence(consensus_2, cluster_2, cluster_kmer_hash);
						
						MACRO_THREAD_UPDATE(this_thread_id, __LINE__)
						
	//					if ( (*cluster_consensus_sequences)[ cluster_1] != NULL ) {
	//						delete (*cluster_consensus_sequences)[ cluster_1];
	//						(*cluster_consensus_sequences)[ cluster_1] = NULL;
	//						
	//					}
	//					if ( (*cluster_consensus_sequences)[ cluster_2] != NULL ) {
	//						delete (*cluster_consensus_sequences)[ cluster_2];
	//						(*cluster_consensus_sequences)[ cluster_2] = NULL;
	//						
	//					}
						
						//  I can delete now both: (either old consensus, or copy of input read)
						if (not consensus_1_isread) {
							consensus_1->delete_data();
							delete consensus_1;
						}
						if (not consensus_2_isread) {
							consensus_2->delete_data();
							delete consensus_2;
						}
						cluster_consensus_sequences->set_writer_lock();
						(*cluster_consensus_sequences)[ cluster_1] = NULL;
						(*cluster_consensus_sequences)[ cluster_2] = NULL;
						cluster_consensus_sequences->unset_writer_lock();
						// visualize cluster:
						//master_cluster_member_list->resetIterator();
						//while(master_cluster_member_list->nextElement()) {
						//	string * tempseq = ((*inputSequences)[master_cluster_member_list->getFirst()].second);
						//	cout << string(master_cluster_member_list->getSecond()+20, '_') << *tempseq << endl;
						//}
						
						
						// create new consensus:
						
						
						mystring * new_consensus_seq = computeConsensus(master_cluster_member_list,
																	  inputSequences, 
																	  alignment_column,
																	  aminoacid_occurence,
																	  aminoacid_scores,
																	  score_matrix,
																	  cluster_aminoacid_counts);
						//cout << "D3" << endl;
						
						cluster_consensus_sequences->set_writer_lock();
						#ifdef DEBUG
						
						
						assert( (*cluster_consensus_sequences)[master_cluster]==NULL );
						
						
						try {
						cluster_consensus_sequences->at(master_cluster) = new_consensus_seq;
						} catch (out_of_range& oor) {
							cerr << "Out of Range error:(cluster_consensus_sequences->at(master_cluster) = new_consensus_seq) " << oor.what() << endl;
							exit(1);
						}	
						#else
						(*cluster_consensus_sequences)[master_cluster] = new_consensus_seq;
						#endif
						cluster_consensus_sequences->unset_writer_lock();
						//cout << "write2 " << new_consensus_seq << " at " << master_cluster << endl;
						//cout <<"merged con: " << *new_consensus_seq << endl;
						//exit(0);
						//cout << "E" << endl;
						
						
						// insert kmers of new_consensus_seq into hash
						
						addConsensusSequence(new_consensus_seq, master_cluster, cluster_kmer_hash);
						
						
#ifdef DEBUG
						// see if old sequence has been added correctly
						testConsensusSequence(new_consensus_seq, master_cluster, cluster_kmer_hash, 54);
						
#endif
						
						previous_master_cluster_index = master_cluster_index;
						#pragma omp atomic
						stat_real_cluster_count--;
		
						//#pragma omp flush
						
					} // end if ( has_overlap)
							
					//omp_unset_lock(my_cluster_lock_1);
					//omp_unset_lock(my_cluster_lock_2);
					
					if ( not has_overlap) {
						// did not have real overlap, try with single merging
						validated_overlaps = 1;
						continue;
					}
					
					break; // end of cluster merging 
					
				} // end if validated_overlaps
				
				break; // leave while(1) loop
			
			} // end while(1)
			
			if (repeat_loop) {
				read_id--;
			} 			
			
			
			
		} // end internal for-loop for input reads

			size_t this_chunk_size = last_read_id - first_read_id +1;
			
			#pragma omp atomic
			reads_processed+=this_chunk_size;
			
			if(this_thread_id == 0)
			{
				//printThreadStatus();
#ifdef TIME
				//overlap_end=clock();
				//overlap_total = overlap_end - overlap_start;
				clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &overlap_end);
				overlap_total = add_time( overlap_total, diff(overlap_start, overlap_end));
#endif
				
				//if (reads_processed % 10000 == 0 && reads_processed > 0) {
					//cerr << "count: " << read_id << endl;
					
					time(&end);
				
					// for more accurate read counts
					size_t report_reads_processed = reads_processed;
					for (int tt = 1; tt < thread_count; ++tt) { // can skip tt=0, since this is always zero.
						report_reads_processed += partial_read_count[tt];
					}
					
				
#pragma omp flush(stat_real_cluster_count)
					log_stream << report_reads_processed << "\t" << stat_real_cluster_count ;
					
					log_stream << "\t" << difftime(end, begin);
					
					
					if (true) { //TODO
						double vm, rss;
						process_mem_usage(vm, rss);
						log_stream << "\t" << (int)vm << "\t" << (int)rss;
						
						
					} else {
						log_stream << "\t-\t-";
					}
					
					log_stream << endl;
					flush(log_stream);
				//}
			}

		} // end while(1) for read chunks
		// this thread has finished his task
		
		
		
		
		
		//if ((read_id >= limit_input_reads) && (limit_input_reads != -1)) {
		//	break;
		//}
		
		
		// CLEAN UP PRIVATE STUFF
		
		delete countOfOverlappingKmers;
		delete lastreadseen;
		
		delete vector_of_offsets;
		delete alignment_column;
		delete aminoacid_occurence;
		delete aminoacid_scores;
			
			
		for (int i = 0 ; i <max_protein_length ; ++i ) {
			delete [] cluster_aminoacid_counts[i];
		}
		delete [] cluster_aminoacid_counts;
	
			
		
	} // end of parallelization #########################
	
	
	if ((int) reads_processed != (int) sequence_count) {
		cerr << "error: reads_processed != sequence_count" << endl;
		cerr << reads_processed << endl;
		cerr << sequence_count << endl;
		exit(1);
	}
	
	if (reads_processed % 10000 != 0 && reads_processed > 0) {
		//cerr << "count: " << read_id << endl;
				
		time(&end);
		#pragma omp flush(stat_real_cluster_count)
		log_stream << reads_processed << "\t" << stat_real_cluster_count ;
		
		log_stream << "\t" << difftime(end, begin);
		
		log_stream << endl;
				
	}
	//double vm, rss;
	process_mem_usage(vm, rss);
	log_stream << "-\t-\t-\t" << (int)vm << "\t" << (int)rss << endl;
	
	
	
	delete fasta_parser;
	//cout << "------------- "<< " X1 " << *(inputSequences->at(0).second) << endl;

	
	
	// --------------------------------------------------------------------------------------------------
	// write consensus sequences to output file
	
	if (this->outputfile.length() == 0) {
		
		string extension = getFileExtension(inputfile);
		string fileNoExt =  getFileNameWithoutExtension(inputfile);
		
		this->outputfile = fileNoExt;
		
		if (extension.compare("faa")== 0 || extension.compare("fa")== 0 || extension.compare("fas")== 0 || extension.compare("fasta")== 0 || extension.compare("txt")== 0) {
			
			
			this->outputfile.append(".clustgun");
		
		
			this->outputfile.append(".");
			this->outputfile.append(extension);
		} else {
			this->outputfile.append(".");
			this->outputfile.append(extension);
			this->outputfile.append(".clustgun");
		}
		
		cout << this->outputfile << endl;
		//exit(0);
		
		
		//this->outputfile = string(inputfile)+string(".consensus");
	}
	log_stream << " done. write results into output file \"" << outputfile << "\"... " << endl;

	if (list_all_members) {
		log_stream << " write member lists into list file \"" << listfile << "\"... " << endl;
	}

	
	//if ((*cluster_consensus_sequences)[1] != NULL) {
	//	cout << "check2: " << *(*cluster_consensus_sequences)[1] << endl;
	//} else {
	//	cout << "is NULL" << endl;
	//	
	//}

	
	
	ofstream output_stream(outputfile.c_str());
	if (! output_stream.is_open())
	{
		cerr << "error writing file " << outputfile << endl;
		exit(1);
	}
	
	ofstream list_stream;
	if (list_all_members) {
		list_stream.open(listfile.c_str());
		if (! list_stream.is_open())
		{
			cerr << "error writing file " << listfile << endl;
			exit(1);
		}

	}
		
	int totalclustercount = last_cluster+1;
	//cout << "maxclustercount: " << maxclustercount << endl;
	//cout << "------------- "<< " X2 " << *(inputSequences->at(0).second) << endl;
	
	
	
	vector<short> * vector_of_offsets;
	vector<char> * alignment_column;
	vector<bool> * aminoacid_occurence;
	vector<int> * aminoacid_scores;
	
	vector_of_offsets = new vector<short>();
	vector_of_offsets->resize(max_protein_length);
	
	alignment_column = new vector<char>();
	alignment_column->resize(1000000);
	
	aminoacid_occurence = new vector<bool>(aminoacid_count, false);
	
	aminoacid_scores = new vector<int>(aminoacid_count, 0);
	
	
	int** cluster_aminoacid_counts = new int* [max_protein_length];
	
	for (int i = 0 ; i <max_protein_length ; ++i ) {
		
		cluster_aminoacid_counts[i]=new int[aminoacid_count];
		int * blabla = cluster_aminoacid_counts[i];
		for (int j = 0; j < aminoacid_count ; ++j  ) {
			blabla[j]=0;
		}
	}
	
	cluster_member_list::iterator mylist_it;
	
	for (int current_cluster = 0 ; current_cluster < totalclustercount; ++current_cluster){
	
		
		if ((*cluster_member_lists)[current_cluster] != NULL) {
			
			// update consensus first:
			
			
			
			
			
			
			mystring * consensus = (*cluster_consensus_sequences)[current_cluster];
			
			
			
			cluster_member_list*& mymemberlist = (*cluster_member_lists)[current_cluster];
			
			
			bool consensus_is_read = false;
			if (consensus == NULL) {
				consensus_is_read = true;
				// cluster of size 1 !
				
				// pointer to read
				consensus = &(*inputSequences)[mymemberlist->start->data_array1[0]].second ;
				
			}
				
			#ifdef DEBUG
			if (consensus == NULL) {
				cerr << "consensus == NULL" << endl;
				exit(1);
				
			}
			#endif
			
			
			// update consensus a last time
			if ( not consensus_is_read) {
				consensus->delete_data();
				delete consensus;
				
				consensus = computeConsensus(mymemberlist,
											  inputSequences, 
											  alignment_column,
											  aminoacid_occurence,
											  aminoacid_scores,
											  score_matrix,
											  cluster_aminoacid_counts);
				
				
				(*cluster_consensus_sequences)[current_cluster]= consensus;
				
			}
			
			
			
			
			output_stream << ">" << prefixname << current_cluster << " length=" << consensus->length() << " size=" << (*cluster_member_lists)[current_cluster]->getLength();
			
			if (true) {
				output_stream << " cov=";
				
				cluster_member_list*& mymemberlist = (*cluster_member_lists)[current_cluster];
				//mymemberlist->resetIterator();
				//bool start_loop = true;
				int tot_len = 0;
				
				for (mylist_it = mymemberlist->begin(); mylist_it != mymemberlist->end(); ++mylist_it) {
				//while (mymemberlist->nextElement()) {
					int seq_id = mylist_it.getFirst();
					//short offset = mymemberlist->getSecond();
					//cout << seq_id << endl;
					//if (start_loop) {
					//	start_loop = false;
					//} else {
					//	output_stream << " " ; 
					//}
					
					tot_len += (*inputSequences)[seq_id].second.length();
					
				}
				
				
				double avgcov = (double)tot_len / (double) consensus->length();
				int avgcon_print = (int)(avgcov * 10); // restrict to one digit after decimal point
				avgcov = (double) avgcon_print / (double) 10;
				output_stream << avgcov ;
			}
			
			
			output_stream << endl; // end of fasta description line
			output_stream << consensus->c_str() << endl; // protein sequence
			
			
			if (list_all_members) {
				list_stream << prefixname << current_cluster << " ";
				
				cluster_member_list*& mymemberlist = (*cluster_member_lists)[current_cluster];
				//mymemberlist->resetIterator();
				bool start_loop = true;
				//while (mymemberlist->nextElement()) {
				
				for (mylist_it = mymemberlist->begin(); mylist_it != mymemberlist->end(); ++mylist_it) {	
					int seq_id = mylist_it.getFirst();
					short offset = mylist_it.getSecond();
					//cout << seq_id << endl;
					if (start_loop) {
						start_loop = false;
					} else {
						list_stream << " " ;
					}
					
					list_stream << "(" << offset << ")" << (*inputSequences)[seq_id].first.c_str();
				}
				list_stream << endl;
				
			}
			
			//if (current_cluster == 98240) {
			//if(true){
			
			//if (mymemberlist->getLength() > 10000) {
			//if (current_cluster == 0) {
			if (false) {
				cluster_member_list*& mymemberlist = (*cluster_member_lists)[current_cluster];
				//mymemberlist->resetIterator();
				//while(mymemberlist->nextElement()) {
				for (mylist_it = mymemberlist->begin(); mylist_it != mymemberlist->end(); ++mylist_it) {
					mystring * tempseq = new mystring(((*inputSequences)[mylist_it.getFirst()].second).c_str());
					cout << string(mylist_it.getSecond()+20, '_') << tempseq->c_str() << endl;
					delete tempseq;
				}
				cout << consensus->c_str() << endl;
				//exit(0);
			}
			
			
			if (not consensus_is_read) {
				consensus->delete_data();
				delete consensus;
				//cluster_consensus_sequences->set_writer_lock();
				(*cluster_consensus_sequences)[current_cluster] = NULL;
				//cluster_consensus_sequences->unset_writer_lock();
			}
			//cout << "> cluster" << current_cluster << endl;
			//cout << *consensus << endl;
		}
		
	}
	
	
	//cout << "------------- "<< " X3 " << *(inputSequences->at(0).second) << endl;
	output_stream.close();
	if (list_all_members) {
		list_stream.close();
	}
	
	
	log_stream << "file written..." << endl;
	
	log_stream << "clean memory..." << endl;
	
	// --- CLEAN UP STUFF --- (only to make valgrind happy)
	
	
	
	
	//cout << cluster_consensus_sequences->name << endl;
	cluster_consensus_sequences->deleteContentPointers();
	//cout << "huhu" << endl;
	delete cluster_consensus_sequences;
	
		
	delete cluster_locks;
		
	//cout << cluster_member_lists->name << endl;
	cluster_member_lists->deleteContentPointers();

	delete cluster_member_lists;
	
	delete vector_of_offsets;
	
	delete alignment_column;
	//cout << "------------- "<< " X4 " << *(inputSequences->at(0).second) << endl;
	//cout << inputSequences->name << endl;
	//inputSequences->deleteContentPairOfPointers();
	delete inputSequences;
	delete inputSequencesData;
	
	for (int i = 0 ; i <max_protein_length ; ++i ) {
			
		delete [] cluster_aminoacid_counts[i];
	}
	delete [] cluster_aminoacid_counts;
		
	#ifndef USE_HashTableSimple
	HashTable::iterator cluster_kmer_hash_it;
	for (cluster_kmer_hash_it = cluster_kmer_hash->begin(); cluster_kmer_hash_it != cluster_kmer_hash->end(); ++cluster_kmer_hash_it) {
		delete cluster_kmer_hash_it->second;
		cluster_kmer_hash_it->second = NULL;
	}
	#endif
	
	delete cluster_kmer_hash;
		
	delete [] score_matrix;
	
	
	
	delete zcat_command;
	delete aminoacid_occurence;
	delete aminoacid_scores;
	log_stream << "cleaning done." << endl;	
	
	date = exec("date");
	log_stream << date << endl;
	log_stream << "end" << endl;
	
	time(&end);
	log_stream << "total time in seconds: " << difftime(end, begin) << endl;
	log_stream << "total time in hours: " << difftime(end, begin)/3600 << endl;
	#ifdef TIME
	
	log_stream << "total time for kmer-iteration: (seconds:nanoseconds) " << kmersearch_total.tv_sec << ":" << kmersearch_total.tv_nsec << endl;
	log_stream << "total time for validation: (seconds:nanoseconds)     " << validation_total.tv_sec << ":" << validation_total.tv_nsec << endl;
	log_stream << "total time for overlap: (seconds:nanoseconds)        " << overlap_total.tv_sec << ":" << overlap_total.tv_nsec << endl;

	int totalclocks=clock();
	log_stream << "total time for everything:     (clock ticks) " << totalclocks << " (seconds) " << (float) totalclocks/CLOCKS_PER_SEC << endl;
	
	//log_stream << "clock: " <<  clock() << " " << clock()/CLOCKS_PER_SEC<< endl;

	#endif
}

void usage(boost::program_options::options_description& options) {
	cerr << endl << "clustgun (git version: " << GIT_REF << " date: " << GIT_DATE << ")" << endl;
	cerr << endl;
	cerr << "Usage: clustgun [--option] input.faa" << endl;
	cerr << endl;
	cerr << options << endl;
	cerr << endl;
}








int main(int argc, const char * argv[])	{
	
	
	
	

	
	// ----------------------------------
	// define program options
	namespace po = boost::program_options;
	
	
	
	string name_help = string("fasta description prefix name (default \"");
			name_help.append(prefixname);
			name_help.append("\")");
	
	
	
	po::options_description options_visible("Options");
	options_visible.add_options()
		
	//	("sort",										"sort input sequences by length (slow!)")
		
		("kmernum",				po::value< int >(&cluster_kmer_overlap_threshold)->default_value(1),	"minimum number of k-mers required")
		("min_overlap_length",	po::value< int >(&min_overlap_length)->default_value(13),				"")
		("kmerlength",			po::value< int >(&kmerlength)->default_value(5) ,						"recommended: 5-7")
		("avgScoreThreshold",	po::value< double >(&avgScoreThreshold)->default_value(3),					"")
		//("windowLength",		po::value< int >(&windowLength)->default_value(10),						"length of sliding window")
		//("windowScoreThreshold",po::value< int >(&windowScoreThreshold)->default_value(0),				"min avg score in sliding window")
		("blosum",				po::value< string >(&blosum_file)->default_value("BLOSUM62"),			"")
		("threads",				po::value< int >(&thread_count_requested)->default_value(0),						"thread count, default(0) all threads")
	
	
		("output",				po::value< string >(),													"output file")
		("name",				po::value< string >(&prefixname)->default_value("cluster"),				"fasta description prefix name")
		("help",																						"display this information");
	
	po::options_description options_hidden("Hidden");
	options_hidden.add_options()
	("input-file", po::value< vector<string> >(), 	"input file");
	
	
	po::options_description all("Allowed options");
	all.add(options_visible).add(options_hidden);
	
	po::positional_options_description p;
	p.add("input-file", -1);
	
	po::variables_map vm;
	try {

		po::store(po::command_line_parser(argc, argv).options(all).positional(p).run(), vm);
		po::notify(vm);
	} catch ( const boost::program_options::error& e ) {
        cerr << "Error: " << e.what() << std::endl;
		exit(1);
    }

	
	// ----------------------------------
	// evaluate program options
	
	string logfile;
	if (vm.count("output")) {
		logfile=vm["output"].as< string >();
	} 
	
	
	if (vm.count("help")) {
		usage(options_visible);
		exit(0);
	}
	
	string input_file;
	string listfile;
	
	if (vm.count("input-file")) {
	
		
		vector< string > input_files_vec =  vm["input-file"].as< vector<string> > ();
		
		if (input_files_vec.size() > 1) {
			cerr << "error: multiple input files. currently only one input file is supported..." << endl;
			exit(1);
		}
		
		input_file = input_files_vec[0];
		
		// multiple input files:
		//cout << "Input files are: " << endl;
		//for (int i = 0 ; i < input_files_vec.size(); ++i) {
		//	cout << input_files_vec[i] << endl;
		//}
		
		if (not vm.count("output") ) {
			//logfile = input_file;
			logfile = getFileNameWithoutExtension(input_file);
			listfile = string(logfile);
			logfile.append(".clustgun.log");
			listfile.append(".clustgun.list");
		} else {
			logfile = getFileNameWithoutExtension(logfile);
			listfile = string(logfile);
			logfile.append(".log");
			listfile.append(".list");
		}
		
		

		
	} else {
		usage(options_visible);
		cerr << endl;
		cerr << "error: input file is missing..." << endl;
		exit(1);
	}

	
	ofstream ofs(logfile.c_str());
	TeeDevice my_tee(cerr, ofs); 
	log_stream.open(my_tee);
	
	
	log_stream << "clustgun (git version: " << GIT_REF << " date: " << GIT_DATE << ")" << endl;
	
	log_stream << "clustgun started with these arguments:" << endl;
	for (int i=0 ; i < argc-1 ; ++i) {
		log_stream << argv[i] << " ";
	}
	log_stream << argv[argc -1] << endl;
	
	
	// print parameters:
	std::map<std::string, boost::program_options::variable_value>::iterator vm_it;
	for ( vm_it=vm.begin() ; vm_it != vm.end(); vm_it++ ) {
		//log_stream << vm_it->second.value().type().name() << endl;
		
		// will not show array like --input
		if (vm_it->second.value().type() == typeid(string) ) {
			if ((*vm_it).second.as<string>().empty()) {
				if (vm.count((*vm_it).first) ) {
					log_stream << (*vm_it).first << " => " << "true" << endl;
				} else {
					log_stream << (*vm_it).first << " => " << "false" << endl;
				}
			} else {
				log_stream << (*vm_it).first << " => " << (*vm_it).second.as<string>() << endl;
			}
		} else if (vm_it->second.value().type() == typeid(int) ) {
			log_stream << (*vm_it).first << " => " << (*vm_it).second.as<int>() << endl;
		} else if (vm_it->second.value().type() == typeid(double) ) {
			log_stream << (*vm_it).first << " => " << (*vm_it).second.as<double>() << endl;
		} else if (vm_it->second.value().type() == typeid(vector<string>))  {
			log_stream << (*vm_it).first << " => " << (*vm_it).second.as<vector<string > >()[0] << endl;
		} else   {
			log_stream << "unknown parameter type" << endl;
		}
		//log_stream << (*vm_it).first << " => " << endl;
		//log_stream << (*vm_it).second->as< string >() << endl;
		//log_stream << (*vm_it).first << " => " << (*vm_it).second << endl;
	}
	log_stream << endl;
	
	
	Clustgun * my_pc = new Clustgun();
	
	my_pc->listfile = string(listfile);
	
	if (vm.count("output")) {
		my_pc->outputfile=vm["output"].as< string >();
	}
	
	//if (vm.count("blosum")) {
	//	blosum_file=vm["output"].as< string >();
	//}

	
	if (min_overlap_length < 1) {
		cerr << "error: parameter does make no sense" << endl;
		exit(1);
	}
	

	
	if (avgScoreThreshold < 0) {
		cerr << "error: parameter does make no sense" << endl;
		exit(1);
	}
	
	
	//if (vm.count("list")) {
	my_pc->list_all_members = true;
	//}
	
	if (vm.count("sort")) {
		my_pc->sort_input_seq = true;
	}
	
	
	

	if (cluster_kmer_overlap_threshold < 1) {
		cerr << "error: parameter does make no sense" << endl;
		exit(1);
	}
	
	// ----------------------------------
	// run
	
	my_pc->cluster(input_file);
	
	delete my_pc;
	log_stream.flush();
    log_stream.close();
		
	return 0;
}

