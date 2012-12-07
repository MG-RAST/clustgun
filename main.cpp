//
//  main.cpp
//  clustgun
//
//  Created by Wolfgang Gerlach on 5/9/12.
//  Copyright (c) 2012 University of Chicago. All rights reserved.
//


#include "main.h"


#ifdef TIME
timespec diff(timespec start, timespec end)
{
	timespec temp;
	if ((end.tv_nsec-start.tv_nsec)<0) {
		temp.tv_sec = end.tv_sec-start.tv_sec-1;
		temp.tv_nsec = 1000000000+end.tv_nsec-start.tv_nsec;
	} else {
		temp.tv_sec = end.tv_sec-start.tv_sec;
		temp.tv_nsec = end.tv_nsec-start.tv_nsec;
	}
	return temp;
}

timespec add_time(timespec time1, timespec time2) {
	timespec result;
	result.tv_sec = time1.tv_sec + time2.tv_sec ;
    result.tv_nsec = time1.tv_nsec + time2.tv_nsec ;
    if (result.tv_nsec >= 1000000000L) {		/* Carry? */
        result.tv_sec++ ;  result.tv_nsec = result.tv_nsec - 1000000000L ;
    }
	
    return (result) ;
	
}


#endif


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





int getOffsets(const char * sequence, int seq_len, vector<short> * vector_of_offsets, sparse_hash_map<int, kmer_appearance_list * , hash<int>, eqint> &cluster_kmer_hash, int max_cluster) {
	
	sparse_hash_map<int, kmer_appearance_list * , hash<int>, eqint>::iterator cluster_kmer_hash_it;
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
	
	while (mykmerit->nextKmer()) {
		int code = mykmerit->code;
		short readpos = mykmerit->kmer_start_pos;
		
		
		cluster_kmer_hash_it = cluster_kmer_hash.find(code);
		
		if (cluster_kmer_hash_it != cluster_kmer_hash.end()) {
			kmer_appearance_list * mylist = cluster_kmer_hash_it->second;
			//cout << "---- "<< endl;
			mylist->resetIterator();
			while (mylist->nextElement()) {
				int clusterhit = mylist->getFirst();
							
				if (clusterhit == max_cluster) {
					number_of_overlapping_kmers_seen++;
					//cout << "match: " << string_int_2_kmer(code) << endl;
					short clusterpos = mylist->getSecond();
					
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
			
			
		} 
		
		
		
	} // end while ( kmer iteration )
	delete mykmerit;
	return number_of_overlapping_kmers_seen;
}

bool sortbyPairSecond( const pair<int, short>& i, const pair<int, short>& j ) {
    return j.second > i.second;
}

//bool sortbyPairSecondLength( const pair<string * , string * >& i, const pair<string * , string * >& j ) {
//    return j.second->length() > i.second->length();
//}


	
string * computeConsensus(cluster_member_list * mymemberlist,
						  HashedArrayTree<pair<mystring, mystring > > * inputSequences,
						  vector<char> * alignment_column,
						  vector<bool> * aminoacid_occurence,
						  vector<int> * aminoacid_scores,
						  short * score_matrix) {
		
	vector<pair<int, short> > sorted_members;
	
	mymemberlist->resetIterator();
	while (mymemberlist->nextElement()){
		#ifdef DEBUG
		try {
		#endif
			sorted_members.push_back(pair<int, short>(mymemberlist->getFirst(), mymemberlist->getSecond()));
		#ifdef DEBUG
		} catch (bad_alloc& ba) {
			cerr << "error: (compteConsensus, push_back) bad_alloc caught: " << ba.what() << endl;
			exit(1);
		}	
		#endif
	}
	
	//for (int member = 0; member < sorted_members.size(); ++member) {
	//	cout << "member: " << member << " offset(org): " << sorted_members[member].second << endl;
	//}
	
	// sort according offset of the members
	sort(sorted_members.begin(),sorted_members.end(),sortbyPairSecond);
	
	//for (int member = 0; member < sorted_members.size(); ++member) {
	//	cout << "member: " << member << " offset(org): " << sorted_members[member].second << endl;
	//}
	
	#ifdef DEBUG
	if (sorted_members.size() == 0) {
		cerr << "error: sorted_members.length() == 0" << endl;
		exit(1);
	}
	#endif
	
	int base_offset = sorted_members[0].second;
	//cout << "base_offset: " << base_offset << endl;
	string * new_consensus_seq;
	#ifdef DEBUG
	try {
	#endif
		new_consensus_seq = new string();
	#ifdef DEBUG
	} catch (bad_alloc& ba) {
		cerr << "error: (computeConsensus) bad_alloc caught: " << ba.what() << endl;
		exit(1);
	}
	#endif
	// update offsets (leftmost member has now offset zero):
	for (int member = 0; member < (int) sorted_members.size(); ++member) {
		//cout << "member: " << member << " offset(org): " << sorted_members[member].second << endl;
		
		#ifdef DEBUG
		int member_offset = sorted_members.at(member).second - base_offset;
		#else
		int member_offset = sorted_members[member].second - base_offset;
		#endif
		
		#ifdef DEBUG
		if (member_offset < 0) {
			cerr <<  "member_offset < 0 : " << member_offset << endl;
			exit(1);
		}
		#endif
		
		sorted_members[member].second = member_offset;
		
		//cout << sorted_members[member].second[] << endl;
		//cout << string(member_offset, '_') << *((*inputSequences)[sorted_members[member].first].second)  << endl;
	}
	// do update also in list !!!
	mymemberlist->resetIterator();
	while (mymemberlist->nextElement()){
		mymemberlist->getSecond() -= base_offset;
	}
	
	
	int current_ali_position = 0;
	int first_member = 0;
	//int last_member;
	
	// compute new consensus
				
	
	while (true) { // iterate columns
		int last_alignment_char = -1;
		bool saw_aa = false;
		for (int member = first_member; member < sorted_members.size(); ++member) {
			if (sorted_members[member].second > current_ali_position) {
				break; // because it is sorted, there can be no further members...
			}
			
			#ifdef DEBUG
			int member_length;
			try {
			 member_length = (inputSequences->at(sorted_members[member].first).second).length();
			} catch (out_of_range& oor) {
				cerr << "Out of Range error:(some_length) " << oor.what() << endl;
				exit(1);
			}
			//if (sorted_members[member].second + some_length - 1 < current_ali_position ){
			#else
			int member_length = ((*inputSequences)[sorted_members[member].first].second).length();
			#endif
			
			if (last_alignment_char >= 100 ) {
				break; // 100 aa's should be enough for consensus
			}
				
			
			
			if (sorted_members[member].second + member_length - 1 < current_ali_position ){
			//if (sorted_members[member].second + ((*inputSequences)[sorted_members[member].first].second)->length() - 1 < current_ali_position ){	
				if (!saw_aa) {
					first_member = member + 1;
				}
				
			} else {
				saw_aa = true;
				//cout << member << endl;
				
				#ifdef DEBUG
				char aa;
				try {
				aa = (inputSequences->at(sorted_members[member].first).second).at(current_ali_position-sorted_members[member].second);
				} catch (out_of_range& oor)	{
					cerr << "Out of Range error:(huhuhuh) " << oor.what() << endl;
					exit(1);
				}

				#else
				char aa = ((*inputSequences)[sorted_members[member].first].second)[current_ali_position-sorted_members[member].second];
				#endif			
				//cout << aa << endl;
				
				last_alignment_char++;
				#ifdef DEBUG
				try {
				alignment_column->at(last_alignment_char) = aa;
				} catch (out_of_range& oor)	{
					cerr << "Out of Range error:(alignment_column->at(last_alignment_char) = aa;) " << oor.what() << endl;
					exit(1);
				}
				#else
				(*alignment_column)[last_alignment_char] = aa;
				#endif
			}
			//cout << string(sorted_members[member].second, '_') << *((*inputSequences)[sorted_members[member].first].second)  << endl;
		}
		//exit(0);
		
		//find consensus aa:
		if (last_alignment_char < 0 ) {
			break;
		}
		
		triplet<bool, char, int> major_aa =  majority_vote<char>(alignment_column, last_alignment_char+1);
		
		#ifdef DEBUG
		char test_char = major_aa.second;
		
		
		if ((int)major_aa.second == 0) {
			cerr << "B    (int)major_aa.second == 0" << endl;
			exit(1);
			
		}
		if ((int)major_aa.second == 0) {
			cerr << "BB    (int)major_aa.second == 0" << endl;
			exit(1);
			
		}
		#endif
			
		if (! major_aa.first) {
			// majority vote did not find a majority, now I use blosum scores to find aminoacid that yields highest score
		
			
			//cerr << "no majority on aa"  << endl;
			//exit(1);
			
			// init:
			for (aminoacid aa = 0; aa < aminoacid_count; ++aa) {
				(*aminoacid_occurence)[aa] = false;
				(*aminoacid_scores)[aa] = 0;
			}
			
			// check which aa are in column
			for (int i = 0; i <= last_alignment_char; ++i) {
				#ifdef DEBUG
				aminoacid aa;
				try {
				aa = aminoacid_ASCII2int[alignment_column->at(i)];
				} catch (out_of_range& oor) {
					cerr << "Out of Range error:(aminoacid aa = aminoacid_ASCII2int[alignment_column->at(i)]) " << oor.what() << endl;
					exit(1);
				}
				
				if (aa != -1 ) {
					try {
						aminoacid_occurence->at(aa) = true;
					} catch (out_of_range& oor) {
						cerr << "Out of Range error:(aminoacid_occurence->at(aa) = true;) " << oor.what() << endl;
						exit(1);
					}
				}
				
					
				#else
				aminoacid aa = aminoacid_ASCII2int[(*alignment_column)[i]];
				if (aa != -1 ) {
					(*aminoacid_occurence)[aa] = true;
				}
				#endif 
				
			}
			#ifdef DEBUG
			if (major_aa.second != test_char ) {
				cerr << "A major_aa.second != test_char" << endl;
				exit(1);
			}
			#endif
			// for each aa count its overall score
			aminoacid max_aa;
			int max_aa_score = INT_MIN;
			for (aminoacid aa = 0; aa < aminoacid_count; ++aa) {
				#ifdef DEBUG
				bool occbool;
				try {
					occbool = aminoacid_occurence->at(aa);
				} catch (out_of_range& oor) {
					cerr << "Out of Range error:(occbool = aminoacid_occurence->at(aa)) " << oor.what() << endl;
					exit(1);
				}	
				if (occbool) {
				#else
				if ((*aminoacid_occurence)[aa]) {
				#endif	
					// walk through cloumn and compute score
					for (int i = 0; i <= last_alignment_char; ++i) {
						#ifdef DEBUG
						char aa_i;
						try {
						aa_i = alignment_column->at(i);
						} catch (out_of_range& oor) {
							cerr << "Out of Range error:(aa_i = alignment_column->at(i);) " << oor.what() << endl;
							exit(1);
						}		
						#else
						char aa_i = (*alignment_column)[i];
						#endif
						
						#ifdef DEBUG
						if (aminoacid_int2ASCII[aa]+256*aa_i >= 256*256) {
							cerr << "error: aminoacid_int2ASCII[aa]+256*aa_i >= 256*256" << endl;
						}
						#endif
						
						#ifdef DEBUG
						try {
						aminoacid_scores->at(aa) += score_matrix[aminoacid_int2ASCII[aa]+256*aa_i];
						} catch (out_of_range& oor) {
							cerr << "Out of Range error:(aminoacid_scores->at(aa) +=) " << oor.what() << endl;
							exit(1);
						}			
						#else
						(*aminoacid_scores)[aa] += score_matrix[aminoacid_int2ASCII[aa]+256*aa_i];
						#endif
						
					}
					//cout << aminoacid_int2ASCII[aa] << ": "<< (*aminoacid_scores)[aa] << endl;
					
					#ifdef DEBUG
					try {
					if (aminoacid_scores->at(aa) > max_aa_score) {
						max_aa_score = aminoacid_scores->at(aa);
						max_aa = aa;
					}
					} catch (out_of_range& oor) {
						cerr << "Out of Range error:(aminoacid_scores->at(aa) > max_aa_score) " << oor.what() << endl;
						exit(1);
					}	
					#else
					if ((*aminoacid_scores)[aa] > max_aa_score) {
						max_aa_score = (*aminoacid_scores)[aa];
						max_aa = aa;
					}
					#endif
					
				}
			}
#ifdef DEBUG			
			if (max_aa_score == INT_MIN) {
				cerr << "max_aa_score == INT_MIN"<< endl;
				
				exit(1);
			}
			
			if ((int)major_aa.second == 0) {
				cerr << "B4    (int)major_aa.second == 0" << endl;
				exit(1);
				
			}
			if (major_aa.second != test_char ) {
				cerr << "major_aa.second != test_char" << endl;
				exit(1);
			}
				
			if (max_aa < 0 ) {
				cerr << "max_aa < 0" << endl;
				exit(1);
			}
#endif			
			//cout << "AA that won: " << aminoacid_int2ASCII[max_aa] << endl;
			major_aa.second = aminoacid_int2ASCII[max_aa];
#ifdef DEBUG			
			if ((int)major_aa.second == 0) {
				cerr << "A    (int)major_aa.second == 0" << endl;
				exit(1);
			}
#endif				
			//exit(0);
		}
		
		
		
#ifdef DEBUG	
		if ((int)major_aa.second == 0) {
			cerr << "D    (int)major_aa.second == 0" << endl;
			cerr << (int) test_char << endl;
			if ( major_aa.first) {
				cerr << "yes" << endl;
			} else {
				cerr << "no" << endl;
			}
			for (int i = 0; i < (last_alignment_char+1); i++) {
				cerr << i << ": " << (*alignment_column)[i] << endl;
			}
			exit(1);
			
		}
#endif	
		
		new_consensus_seq->append(1, major_aa.second);
#ifdef DEBUG			
		if ((int)major_aa.second == 0) {
			cerr << "C    (int)major_aa.second == 0" << endl;
			
			exit(1);
			
		}
#endif			
		if (first_member >= sorted_members.size()) {
			
			break;
		}
		
		//cout << "--" << endl;
		current_ali_position++;
	}
	
#ifdef DEBUG	
	for(int i = 0 ; i < new_consensus_seq->length(); ++i) {
		if ( ((int) (*new_consensus_seq)[i]) == 0 ) {
			cerr << "error: (*new_consensus_seq)[i] == 0"<< endl;
			cerr << *new_consensus_seq << endl;
			cerr << i << endl;
			
			exit(1);
		}
		
	}
#endif		
	return new_consensus_seq;
}

void addConsensusSequence(const char * sequence, int seq_len, int cluster, sparse_hash_map<int, kmer_appearance_list * , hash<int>, eqint>& cluster_kmer_hash) {

	
	sparse_hash_map<int, kmer_appearance_list * , hash<int>, eqint>::iterator cluster_kmer_hash_it;
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
	mykmer_insertion_it->reset(0);

	while (mykmer_insertion_it->nextKmer()) {
		int code = mykmer_insertion_it->code;
		int pos = mykmer_insertion_it->kmer_start_pos;
		
		cluster_kmer_hash_it = cluster_kmer_hash.find(code);
		kmer_appearance_list * mylist;
		if (cluster_kmer_hash_it != cluster_kmer_hash.end()) {
			mylist = cluster_kmer_hash_it->second;
			
			
			//cout << "append to old list: " << string_int_2_kmer(code) << endl;
			//exit(0);
			
		} else {
			#ifdef DEBUG
			try {
#endif
				mylist = new kmer_appearance_list();
				#ifdef DEBUG
			} catch (bad_alloc& ba) {
				cerr << "error: (new kmer_appearance_list) bad_alloc caught: " << ba.what() << endl;
				exit(1);
			}		
#endif
			cluster_kmer_hash[code]=mylist;
			//cout << "add new list" << endl;
		}
		mylist->append(cluster,pos);
		//if (code == 3566)
		//	cout << "insertX: "<< clusterhit << " " << pos << endl;
		//if (code == 1579170)
		//	cout << "insert: "<< clusterhit << " " << pos << endl;
	}// edn while

	delete mykmer_insertion_it;
}

// removes kmers from index
void removeConsensusSequences(string * consensus_old, int cluster, sparse_hash_map<int, kmer_appearance_list * , hash<int>, eqint>& cluster_kmer_hash) {

	sparse_hash_map<int, kmer_appearance_list * , hash<int>, eqint>::iterator cluster_kmer_hash_it;
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
		
		cluster_kmer_hash_it = cluster_kmer_hash.find(code);
		kmer_appearance_list * mylist;
		if (cluster_kmer_hash_it != cluster_kmer_hash.end()) {
			
			mylist = cluster_kmer_hash_it->second;
			mylist->resetIterator();
			#ifdef DEBUG
			bool found_kmer = false;
			#endif
			while (mylist->nextElement()) {
				//if (cluster == 606 && code == 2644139) {cout << "mylist->currentArrayPosition: " << mylist->currentArrayPosition << " mylist->lastArrayPosition: " << mylist->lastArrayPosition << " mylist->getFirst(): "<< mylist->getFirst() << endl; }
				if (mylist->getFirst() == cluster) {
					//if (cluster == 606 && code == 2644139) {mylist->print();}
					#ifdef DEBUG
					found_kmer = true;
					#endif
					mylist->eraseElement();
					//if (cluster == 606 && code == 2644139) {cout << "after erase element" << endl; mylist->print();}
					//if (mylist->getLength() > 1) {
					//exit(0);
					//}
					break;
				}
				
			}
			#ifdef DEBUG
			if (! found_kmer) {
				cerr << "error: did not find kmer I wanted to delete! cluster: " << cluster << " kmer: " << code << " " << string_int_2_kmer(code) << endl;
				cerr << *consensus_old << endl;
				mylist->print();
				exit(1);
			}
			#endif
			//cout << "append to old list: " << string_int_2_kmer(code) << endl;
			//exit(0);
			
		} else {
			cerr << "error: kmer for deletion not found!" << endl;
			exit(1);
			//cout << "add new list" << endl;
		}
		//mylist->append(last_cluster,pos);
	} // end while
	
	delete mykmer_deletion_it;


}


bool computeSequenceOverlap(int offset, string * a, string * b, short * score_matrix, int majority_vote_count) {
	int i = max(0, offset);
	int j = i - offset;
	
	//cout << "i: " << i << endl;
	//cout << "j: " << j << endl;
	
	if (false) {
		//int shift = 0;
		//if (j<0) {
		//shift = -j;
		//}
		cout << string(j, ' ') << *a << endl;
		cout << string(i, ' ') << *b << endl;
	}
	
	//if (j<0) {
	//	exit(0);	
	//}
	
	
	
	int a_len = a->length();
	int b_len = b->length();
	
	int overlap_length = min((a_len - i), (b_len - j));
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
	
	
	int start_i = i;
	int end_i = i+overlap_length;
	for (int i = start_i; i< end_i ; ++i) {
		tot_len++;
		j = i - offset;
		
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


void sortInputSequences(HashedArrayTree<pair<string *, string * > > * inputSequences) {
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
				pair<string * , string * > pp = (*inputSequences)[j];
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
		aminoacid_ASCII2int[aa]=i;
		aa = tolower ( aa );
		aminoacid_ASCII2int[aa]=i;
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


	
		
	
	//HashedArrayTree<string * > * myHAT = new HashedArrayTree<string * >(3, new string("empty"));
	//cout << "gotA: " << myHAT->hash->at(0) << endl;
	
	//myHAT->push_back(new string("test1"));
	
	//myHAT->push_back(new string("test2"));
	//myHAT->push_back(new string("test3"));
	//myHAT->push_back(new string("test4"));
	//myHAT->push_back(new string("test5"));
//	myHAT->push_back(new string("test6"));
//	myHAT->push_back(new string("test7"));
//	myHAT->push_back(new string("test8"));
//	myHAT->push_back(new string("test9"));
//	
	//cout << "gotA: " << *(*myHAT)[0] << endl;
	//cout << "gotA: " << *(myHAT->at(100)) << endl;
	//cout << "gotA: " << myHAT->hash->at(0) << endl;
	//cout << "got: " << myHAT->hash->at(1) << endl;
	
	//cout << "gotB: " << myHAT->hash->at(0)->size() << endl;
	//cout << "gotB: " << *(myHAT->hash->at(0)->at(0)) << endl;
	//cout << "got: " << myHAT->hash->at(0)->at(1) << endl;
	
	//exit(0);
	
	sparse_hash_map<int, int, hash<int>, eqint> kmer_count_hash;
	sparse_hash_map<int, int, hash<int>, eqint>::iterator it;
	
	FASTA_Parser * fasta_parser;
	
	int total_substrings=0;
	
	
	int total_read_count = 0;
	string descr;
	string * fasta_sequence;

	// ----------------------------------------
	// read sequences into memory
	
	HashedArrayTreeString * inputSequencesData = new HashedArrayTreeString(23); // 2^23 = 8MB chunks
	
	
	HashedArrayTree<pair<mystring, mystring > > * inputSequences;
	#ifdef DEBUG
	try {
	#endif
		inputSequences = new HashedArrayTree<pair<mystring, mystring > >(20); // 20 for 2^20=1MB chunks
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
		
		
		
		if (descr.length() > 0 && fasta_sequence->length() > min_overlap_length) { // won't make sense to keep sequences short as the min_overlap_length
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
			
			delete fasta_sequence;
			
			//cout << "pushed: " << inputSequences->size() << endl;	
			//cout << *sequence << endl;
			//exit(0);
							
			   
			if (total_read_count % 1000000 == 0 ) {
				cerr << "total_read_count: " << total_read_count << endl;
			}
						   
			
		} else {
			ignore_short_seq++;
			//cerr << "warning: sequence was not accepted..." << endl;
			delete fasta_sequence;
		}
		
		
		//if (total_read_count >= limit_input_reads && limit_input_reads != -1) {
		//	cerr << "WARNING: numer of reads for debugging purposes limited !!!!!!" << endl;
		//	break;
		//}
	}
	
	if (ignore_short_seq > 0 ) {
		log_stream << "warning: ignored " << ignore_short_seq << " sequences, because they were shorter than the mininmal overlap length." << endl;
	}
	
	cerr << "total_read_count: " << total_read_count << endl;
	

	
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
	
	
	// -------------------------------------------------------------------------
	
	
		
	sparse_hash_map<int, kmer_appearance_list * , hash<int>, eqint> cluster_kmer_hash;
	sparse_hash_map<int, kmer_appearance_list * , hash<int>, eqint>::iterator cluster_kmer_hash_it;
	
	
	// --------- CLUSTERS ----------- // cluster objects would have been nice, but I am afraid of memory inefficiency
	HashedArrayTree<short> * countOfOverlappingKmers;
	#ifdef DEBUG
	try {
	#endif
		countOfOverlappingKmers = new HashedArrayTree<short>(20, 0);
	#ifdef DEBUG
	} catch (bad_alloc& ba) {
		cerr << "error: (countOfOverlappingKmers) bad_alloc caught: " << ba.what() << endl;
		exit(1);
	}	
	#endif
	
	#ifdef DEBUG
	countOfOverlappingKmers->name = string("countOfOverlappingKmers");
	#endif
	HashedArrayTree<int> * lastreadseen;
	#ifdef DEBUG
	try {
	#endif
		lastreadseen = new HashedArrayTree<int>(20, -1); // this avoids initialization
	#ifdef DEBUG
	} catch (bad_alloc& ba) {
		cerr << "error: (lastreadseen) bad_alloc caught: " << ba.what() << endl;
		exit(1);
	}	
	#endif
	
	#ifdef DEBUG
	lastreadseen->name = string("lastreadseen");
	#endif
	HashedArrayTree<string * > * cluster_consensus_sequences;
	#ifdef DEBUG
	try {
	#endif
		cluster_consensus_sequences = new HashedArrayTree<string * >(20, NULL);
	#ifdef DEBUG
	} catch (bad_alloc& ba) {
		cerr << "error: (cluster_consensus_sequences) bad_alloc caught: " << ba.what() << endl;
		exit(1);
	}
	#endif
	
	#ifdef DEBUG
	cluster_consensus_sequences->name = string("cluster_consensus_sequences");
	#endif	
	//HashedArrayTree<short> * cluster_consensus_offsets = new HashedArrayTree<short>(20, 0); // offset relative to seed read
	
	HashedArrayTree<cluster_member_list * > * cluster_member_lists;
	#ifdef DEBUG	
	try {
	#endif		
		cluster_member_lists = new HashedArrayTree<cluster_member_list * >(20, NULL);
	#ifdef DEBUG		
	} catch (bad_alloc& ba) {
		cerr << "error: (cluster_member_lists) bad_alloc caught: " << ba.what() << endl;
		exit(1);
	}
	#endif	
	#ifdef DEBUG
	cluster_member_lists->name = string("cluster_member_lists");
	#endif
	
	vector<short> * vector_of_offsets = new vector<short>();
	vector_of_offsets->resize(max_protein_length);
	
	vector<char> * alignment_column = new vector<char>();
	alignment_column->resize(1000000);
	
	vector<bool> * aminoacid_occurence = new vector<bool>(aminoacid_count, false);
	vector<int> * aminoacid_scores = new vector<int>(aminoacid_count, 0);
	
	int last_cluster=-1;

	
	
	
	countOfOverlappingKmers->reserve(hat_increase_steps);
	lastreadseen->reserve(hat_increase_steps);
	cluster_consensus_sequences->reserve(hat_increase_steps);
	cluster_member_lists->reserve(hat_increase_steps);
	
	// -------------------------------------------------------------------------
	
	//iterate through all reads to create clusters and match against existing clusters
	log_stream << "#reads\t#clusters";
	
	log_stream << "\tseconds";
	
	log_stream << "\tVM[kb]\tRSS[kb]";
	
	log_stream << endl;


	
	process_mem_usage(vm, rss);
	log_stream << "0\t0\t0\t" << (int)vm << "\t" << (int)rss << endl;	
	
	
	size_t sequence_count = inputSequences->size();
	
	if (total_read_count != sequence_count) {
		cerr << "total_read_count != sequence_count" << endl;
		exit(1);
	}
	
	
	int stat_real_cluster_count=0;
	
    const int top_n_clusters = 3;
    int max_cluster_id_array[top_n_clusters]; 
    int max_cluster_kmer_count_array[top_n_clusters];
    short max_cluster_offsets[top_n_clusters];

	mystring * sequence = 0;
	
	for (size_t read_id = 0 ; read_id < sequence_count; read_id++){
	
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
		clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &kmersearch_start);
		#endif
		while (mykmerit->nextKmer()) {
			int code = mykmerit->code;
			//int pos = mykmerit->kmer_start_pos;
			
			//cout << "here code: " << code<< endl;
			cluster_kmer_hash_it = cluster_kmer_hash.find(code);
			
			// check if k-mer is indexed:
			if (cluster_kmer_hash_it != cluster_kmer_hash.end()) {
				kmer_appearance_list * mylist = cluster_kmer_hash_it->second;
				//cout << "---- "<< endl;
				mylist->resetIterator();
				// check all occurrences of that k-mer
				while (mylist->nextElement()) {
					int clusterhit = mylist->getFirst();
					
					
					#ifdef DEBUG
					if (clusterhit > last_cluster) {
						cerr << "clusterhit > last_cluster" << endl;
						cerr << "clusterhit: " << clusterhit << endl;
						cerr << "last_cluster: " << last_cluster << endl;
						exit(1);
					}
					
					if (clusterhit > lastreadseen->capacity()) {
						cerr << "clusterhit > lastreadseen->capacity()" << endl;
						cerr << "clusterhit: " << clusterhit << endl;
						cerr << "lastreadseen->capacity(): " << lastreadseen->capacity() << endl;
						exit(1);
					}
					#endif
					
					
					
					//cerr << "clusterhit: " << clusterhit << endl;
					if ((*lastreadseen)[clusterhit] == read_id) {
						#ifdef DEBUG
						if (clusterhit > countOfOverlappingKmers->capacity()) {
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
						if (clusterhit >= lastreadseen->capacity() ){
							cerr << "clusterhit >= lastreadseen->capacity()" << endl;
							cerr << clusterhit << " " << lastreadseen->capacity() << endl;
							exit(1);
						}
						
						if (clusterhit >= countOfOverlappingKmers->capacity() ){
							cerr << "clusterhit >= countOfOverlappingKmers->capacity()" << endl;
							cerr << clusterhit << " " << countOfOverlappingKmers->capacity() << endl;
							exit(1);
						}
						
						#endif
						
						(*lastreadseen)[clusterhit] = read_id;
						(*countOfOverlappingKmers)[clusterhit]=1;
						//cout << "not seen before" << endl;
					}
				} // end while
				
			} 				
			
			
		} // end while ( kmer iteration )
		#ifdef TIME
		clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &kmersearch_end);
		kmersearch_total = add_time( overlap_total, diff(kmersearch_start, kmersearch_end));
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
		clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &validation_start);
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
				
				
				if (max_cluster_id_array[cluster_it] == -1) {
					break;
				}
				
				// store overlaps in vector_of_offsets
				int number_of_overlapping_kmers_seen = getOffsets(sequence->c_str(),
																  sequence->length(),
																  vector_of_offsets,	// <-- output
																  cluster_kmer_hash,
																  max_cluster_id_array[cluster_it]);
				
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
				
				
				
				if (overlap_found) {
		
					//check overlap:
					int cluster = max_cluster_id_array[cluster_it];
					
					#ifdef DEBUG
					string * consensus_seq;
					try {
					consensus_seq = cluster_consensus_sequences->at(cluster);
					} catch (out_of_range& oor) {
						cerr << "Out of Range error:(consensus_seq = cluster_consensus_sequences->at(cluster)) " << oor.what() << endl;
						exit(1);
					}
					#else
					string * consensus_seq = (*cluster_consensus_sequences)[cluster];
					#endif
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
						try {
						consensus_seq = new string(inputSequences->at(cluster_read).second.c_str());
						} catch (out_of_range& oor) {
							cerr << "Out of Range error:(consensus_seq = inputSequences->at(cluster_read).second) " << oor.what() << endl;
							exit(1);
						}
						#else
						consensus_seq = new string((*inputSequences)[cluster_read].second.c_str());
						#endif
					}
					
					string help_sequence = string(sequence->c_str());
					overlap_found = computeSequenceOverlap(majority_offset, consensus_seq, &help_sequence, score_matrix, major.third);
					
				
					//exit(1);
					//computeSequenceOverlap(int offset, string * a, string * b, short * score_matrix)
				}
				
				if (overlap_found) {
					validated_overlaps++;
					//if (read_id == 414) {
					//	cout<< "::: " << cluster_it  << " " << validated_overlaps << endl;
					//}
					if (cluster_it != (validated_overlaps-1)) {
						max_cluster_id_array[validated_overlaps-1] = max_cluster_id_array[cluster_it];
						max_cluster_kmer_count_array[validated_overlaps-1] = max_cluster_kmer_count_array[cluster_it];
					}
					max_cluster_offsets[validated_overlaps-1] = majority_offset;
					
				} else {
					max_cluster_id_array[cluster_it] = -1;
					max_cluster_kmer_count_array[cluster_it] = 0;
				}
				
			}
			
			
			
				
			
			
			
			
		} // end if (max_cluster_id_array[0] >= 0)
		
		#ifdef TIME
		clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &validation_end);
		validation_total = add_time( validation_total, diff(validation_start, validation_end));
		//cerr << validation_start << " " << validation_end << endl;
		#endif
		
		
		
		// ************ process overlaps ************
		//cout << "validated_overlaps"  << endl;
		#ifdef TIME
		clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &overlap_start);
		#endif
		while (true) {
			if (validated_overlaps == 0) {
				
				// did not found matches, insert as new cluster
				stat_real_cluster_count++;
				last_cluster++;
				
				// make sure the HAT-arrays have enough memory allocated
				// this has the advantage I can use []-operator instead of at()-function
				
				if ((last_cluster % hat_increase_steps) == 0) {
					int reserve_max = last_cluster+hat_increase_steps;
					
					countOfOverlappingKmers->reserve(reserve_max);
					lastreadseen->reserve(reserve_max);
					cluster_consensus_sequences->reserve(reserve_max);
					cluster_member_lists->reserve(reserve_max);
					
				}
				
				
				// insert seed as first member:
				
				#ifdef DEBUG
				
				cluster_member_list * templist;
				try {
					templist = cluster_member_lists->at(last_cluster);
				} catch (out_of_range& oor) {
					cerr << "Out of Range error:(templist = cluster_member_lists->at(last_cluster)) " << oor.what() << endl;
					exit(1);
				}
				
				cluster_member_list*& mymemberlist = cluster_member_lists->at(last_cluster);
				#else
				cluster_member_list*& mymemberlist = (*cluster_member_lists)[last_cluster]; // alias to a pointer
				#endif			
				//cout << mymemberlist << endl;
				assert( (mymemberlist == NULL) );
				#ifdef DEBUG					
				try {
				#endif
					mymemberlist = new cluster_member_list();
				#ifdef DEBUG
				} catch (bad_alloc& ba) {
					cerr << "error: (new cluster_member_list) bad_alloc caught: " << ba.what() << endl;
					exit(1);
				}
				#endif				
				mymemberlist->append(read_id, 0);
				
				
				
				addConsensusSequence(sequence->c_str(), sequence->length(), last_cluster, cluster_kmer_hash);
				
				break;
				
			} else if (validated_overlaps == 1) { // read has match to exactly one cluster
				
				int clusterhit = max_cluster_id_array[0];
				
				
				#ifdef DEBUG
				
				cluster_member_list * templist;
				try {
					templist = cluster_member_lists->at(clusterhit);
				} catch (out_of_range& oor) {
					cerr << "Out of Range error:(templist = cluster_member_lists->at(clusterhit)) " << oor.what() << endl;
					exit(1);
				}
				
				cluster_member_list*& mymemberlist = cluster_member_lists->at(clusterhit);
				#else
				cluster_member_list*& mymemberlist = (*cluster_member_lists)[clusterhit]; // alias to a pointer
				#endif	
								
				if (mymemberlist == NULL){
					#ifdef DEBUG
					try {
					#endif
						mymemberlist = new cluster_member_list();
					#ifdef DEBUG
					} catch (bad_alloc& ba) {
						cerr << "error: (new cluster_member_list) bad_alloc caught: " << ba.what() << endl;
						exit(1);
					}
					#endif
				}
				
				mymemberlist->append(read_id, max_cluster_offsets[0]);
				
				
				// find old consensus sequence
				#ifdef DEBUG
				string * consensus_old;
				try {
				consensus_old = cluster_consensus_sequences->at(clusterhit);
				} catch (out_of_range& oor) {
					cerr << "Out of Range error:(consensus_old = cluster_consensus_sequences->at(clusterhit)) " << oor.what() << endl;
					exit(1);
				}
				#else
				string * consensus_old = (*cluster_consensus_sequences)[clusterhit];
				#endif				
				if (consensus_old == NULL) {
					// consensus seuquence does not exist yet, seed member had been used
					int seed_id = mymemberlist->start->data_array1[0];
					#ifdef DEBUG
					try {
					consensus_old = inputSequences->at(seed_id).second;
					} catch (out_of_range& oor) {
						cerr << "Out of Range error:(consensus_old = inputSequences->at(seed_id).second) " << oor.what() << endl;
						exit(1);
					}
					#else
					consensus_old = new string((*inputSequences)[seed_id].second.c_str());
					#endif
				}
				
				
				
				// check if new member is a perfect full-length overlap
				bool extends_consensus = true;
				
				#ifdef DEBUG
				int member_len;
				try {
				member_len = inputSequences->at(read_id).second->length();
				} catch (out_of_range& oor) {
					cerr << "Out of Range error:(member_len = inputSequences->at(read_id).second->length()) " << oor.what() << endl;
					exit(1);
				}
				#else
				int member_len = (*inputSequences)[read_id].second.length();
				#endif					
				//bool perfect_overlap = false;
				if (max_cluster_offsets[0] >= 0 && max_cluster_offsets[0]+member_len <= consensus_old->length() ) { // check for full-length overlap
					extends_consensus = false;
				} 
				
				
				int membercount = mymemberlist->getLength();
				//cout << "membercount: " << membercount << endl;
				
				
				if (clusterNeedsUpdate(membercount) || extends_consensus) { // mymemberlist->getLength() >= 100
				//if (true) {	
					
					
					
					// ============ REMOVE OLD CONSENSUS K-MERS =============
					
					
					
					
					// remove consensus sequence kmers for cluster "clusterhit"
					removeConsensusSequences(consensus_old, clusterhit, cluster_kmer_hash);
					
					
					// delete old consensus, but not seed !
					
					#ifdef DEBUG
					string * tempstr;
					try {
					tempstr = cluster_consensus_sequences->at(clusterhit);
					} catch (out_of_range& oor) {
						cerr << "Out of Range error:(tempstr = cluster_consensus_sequences->at(clusterhit)) " << oor.what() << endl;
						exit(1);
					}
					if (tempstr != NULL) {
						delete consensus_old;
						cluster_consensus_sequences->at(clusterhit) = NULL;
					}
					#else					
					if ((*cluster_consensus_sequences)[clusterhit] != NULL) {
						delete consensus_old;
						(*cluster_consensus_sequences)[clusterhit] = NULL;
					}
					#endif					
					//consensus_old = NULL;					
					
					
					// ============= SORT, UPDATE OFFSETS AND COMPUTE CONSENSUS ==============
					// copy members from list to array (temporarly)
					
					string * new_consensus_seq = computeConsensus(mymemberlist,  
																  inputSequences, 
																  alignment_column,
																  aminoacid_occurence,
																  aminoacid_scores,
																  score_matrix);
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
					//cout << "write " << new_consensus_seq << " at " << clusterhit << endl;
					//cluster_consensus_sequences->
					
					//cout <<"con: " << *new_consensus_seq << endl;
					
					
					
					
					// insert kmers of new_consensus_seq into hash
					
					addConsensusSequence(new_consensus_seq->c_str(), new_consensus_seq->length(), clusterhit, cluster_kmer_hash);
					
					
					
				}
				break;
				
			} else { //if (validated_overlaps >= 2) {
				
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
					
				// Cluster 1
				int cluster_1 = max_cluster_id_array[previous_master_cluster_index];
				cluster_member_list * cluster_1_memberlist = (*cluster_member_lists)[cluster_1];
				string * consensus_1 = (*cluster_consensus_sequences)[cluster_1];
				
				if (consensus_1 == NULL) {
					int cluster_read = cluster_1_memberlist->start->data_array1[0];
//					#ifdef DEBUG
//					consensus_1 = inputSequences->at(cluster_read).second;
//					#else
//					consensus_1 = (*inputSequences)[cluster_read].second;
//					#endif
					consensus_1 = new string((*inputSequences)[cluster_read].second.c_str());
				}
				
				// Cluster 2
				int cluster_2 = max_cluster_id_array[second_cluster_index];
				cluster_member_list * cluster_2_memberlist = (*cluster_member_lists)[cluster_2];
				string * consensus_2 = (*cluster_consensus_sequences)[cluster_2];
				
				
				if (consensus_2 == NULL) {
					int cluster_read = cluster_2_memberlist->start->data_array1[0];
//					#ifdef DEBUG
//					consensus_2 = inputSequences->at(cluster_read).second;
//					#else
//					consensus_2 = (*inputSequences)[cluster_read].second;
//					#endif
					consensus_2 = new string((*inputSequences)[cluster_read].second.c_str());
				}
				
				//cout << cluster_1 << endl;
				//cout << cluster_2 << endl;
				
				//cout << max_cluster_offsets[previous_master_cluster_index] << endl;
				//cout << max_cluster_offsets[second_cluster_index] << endl;
				
				int offset_diff = max_cluster_offsets[previous_master_cluster_index] - max_cluster_offsets[second_cluster_index];
				//cout << "seq: " << *sequence << endl;
				//cout << "o: " << offset_diff << endl;
				
				
				
				bool has_overlap = computeSequenceOverlap(offset_diff, consensus_1, consensus_2, score_matrix, 0);
				
				//cout << (has_overlap?string("true"):string("false")) << endl;
				
				if (! has_overlap) {
					// in this case we fall back to the simple procedure of adding one read to one cluster
					validated_overlaps = 1;
					
				} else {
				
					
					
										
					
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
						
						master_cluster_member_list = (*cluster_member_lists)[cluster_1];
						slave_cluster_member_list = (*cluster_member_lists)[cluster_2];
						
						
					} else {
						
						master_cluster_index =  second_cluster_index;
						slave_cluster_index = previous_master_cluster_index;
						
						master_cluster_member_list = (*cluster_member_lists)[cluster_2];
						slave_cluster_member_list = (*cluster_member_lists)[cluster_1];
						
					}
					
					master_cluster = max_cluster_id_array[master_cluster_index];
					slave_cluster = max_cluster_id_array[slave_cluster_index];
					
					//master_consensus_sequence = consensus_2;
					//slave_consensus_sequence = consensus_1;
					
					slave_offset = max_cluster_offsets[master_cluster_index] - max_cluster_offsets[slave_cluster_index];
					newread_offset = max_cluster_offsets[master_cluster_index];
					
					//cout << "C" << endl;
					
					// update and copy slave members to master
					
					slave_cluster_member_list->resetIterator();
					while (slave_cluster_member_list->nextElement()) {
						master_cluster_member_list->append(slave_cluster_member_list->getFirst(), slave_cluster_member_list->getSecond() + slave_offset);
					}
					
					// append current read as member:
					master_cluster_member_list->append(read_id, newread_offset);
					
					// delete stuff
					delete slave_cluster_member_list;
					(*cluster_member_lists)[slave_cluster] = NULL;
					
					removeConsensusSequences(consensus_1, cluster_1, cluster_kmer_hash);
					removeConsensusSequences(consensus_2, cluster_2, cluster_kmer_hash);
					
					
					
					
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
					
					// since both are copies, I can delete now both:
					delete consensus_1;
					delete consensus_2;
					(*cluster_consensus_sequences)[ cluster_1] = NULL;
					(*cluster_consensus_sequences)[ cluster_2] = NULL;
					
					// visualize cluster:
					//master_cluster_member_list->resetIterator();
					//while(master_cluster_member_list->nextElement()) {
					//	string * tempseq = ((*inputSequences)[master_cluster_member_list->getFirst()].second);
					//	cout << string(master_cluster_member_list->getSecond()+20, '_') << *tempseq << endl;
					//}
					
					
					// create new consensus:
					
					
					string * new_consensus_seq = computeConsensus(master_cluster_member_list,  
																  inputSequences, 
																  alignment_column,
																  aminoacid_occurence,
																  aminoacid_scores,
																  score_matrix);
					//cout << "D3" << endl;
					assert( (*cluster_consensus_sequences)[master_cluster]==NULL );
					
					#ifdef DEBUG
					try {
					cluster_consensus_sequences->at(master_cluster) = new_consensus_seq;
					} catch (out_of_range& oor) {
						cerr << "Out of Range error:(cluster_consensus_sequences->at(master_cluster) = new_consensus_seq) " << oor.what() << endl;
						exit(1);
					}	
					#else
					(*cluster_consensus_sequences)[master_cluster] = new_consensus_seq;
					#endif
					//cout << "write2 " << new_consensus_seq << " at " << master_cluster << endl;
					//cout <<"merged con: " << *new_consensus_seq << endl;
					//exit(0);
					//cout << "E" << endl;
					
					
					// insert kmers of new_consensus_seq into hash
					
					addConsensusSequence(new_consensus_seq->c_str(), new_consensus_seq->length(), master_cluster, cluster_kmer_hash);
					
					previous_master_cluster_index = master_cluster_index;
					stat_real_cluster_count--;
					
					break;				
				}
				
				
				
				
				
				
				//cout << "leave merge process " << endl;
				//exit(0);
			} 		
			
			
		  
		} // end while true	
		
		#ifdef TIME
		//overlap_end=clock();
		//overlap_total = overlap_end - overlap_start;
		clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &overlap_end);
		overlap_total = add_time( overlap_total, diff(overlap_start, overlap_end));
		#endif
		
		if ((read_id +1) % 10000 == 0 && read_id > 0) {
			//cerr << "count: " << read_id << endl;
				
			time(&end);
					
			
			log_stream << read_id+1 << "\t" << stat_real_cluster_count ;
			
			log_stream << "\t" << difftime(end, begin);
			
			
			if ((read_id +1) % 500000 == 0 && read_id > 0) {
				double vm, rss;
				process_mem_usage(vm, rss);
				log_stream << "\t" << (int)vm << "\t" << (int)rss;

				 
			} else {
				log_stream << "\t-\t-";
			}
			
			log_stream << endl;
			flush(log_stream);
		}
		//if ((read_id >= limit_input_reads) && (limit_input_reads != -1)) {
		//	break;
		//}
	}
	
	if (sequence_count % 10000 != 0 && sequence_count > 0) {
		//cerr << "count: " << read_id << endl;
				
		time(&end);
		
		log_stream << sequence_count << "\t" << stat_real_cluster_count ;
		
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
	
	
	for (int current_cluster = 0 ; current_cluster < totalclustercount; ++current_cluster){
		if ((*cluster_member_lists)[current_cluster] != NULL) {
			
			// update consensus first:
			
			
			
			
			
			
			string * consensus = (*cluster_consensus_sequences)[current_cluster];
			cluster_member_list*& mymemberlist = (*cluster_member_lists)[current_cluster];
			
			bool consensus_is_read = false;
			if (consensus == NULL) {
				consensus_is_read = true;
				// cluster of size 1 !
				
				
				consensus = new string((*inputSequences)[mymemberlist->start->data_array1[0]].second.c_str()) ;
			}
				
			#ifdef DEBUG
			if (consensus == NULL) {
				cerr << "consensus == NULL" << endl;
				exit(1);
				
			}
			#endif
			
			
			//if (not clusterNeedsUpdate(mymemberlist->getLength()) ) {
			if ( not consensus_is_read) {	
				delete consensus;
				
				consensus = computeConsensus(mymemberlist,  
											  inputSequences, 
											  alignment_column,
											  aminoacid_occurence,
											  aminoacid_scores,
											  score_matrix);
				
				
				(*cluster_consensus_sequences)[current_cluster]= consensus;
				
			}
			
			
			
			
			output_stream << ">" << prefixname << current_cluster << " length=" << consensus->length() << " size=" << (*cluster_member_lists)[current_cluster]->getLength();
			
			if (true) {
				output_stream << " cov=";
				
				cluster_member_list*& mymemberlist = (*cluster_member_lists)[current_cluster];
				mymemberlist->resetIterator();
				//bool start_loop = true;
				int tot_len = 0;
				while (mymemberlist->nextElement()) {
					int seq_id = mymemberlist->getFirst();
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
			output_stream << *consensus << endl; // protein sequence
			
			
			if (list_all_members) {
				list_stream << prefixname << current_cluster << " ";
				
				cluster_member_list*& mymemberlist = (*cluster_member_lists)[current_cluster];
				mymemberlist->resetIterator();
				bool start_loop = true;
				while (mymemberlist->nextElement()) {
					int seq_id = mymemberlist->getFirst();
					short offset = mymemberlist->getSecond();
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
				mymemberlist->resetIterator();
				while(mymemberlist->nextElement()) {
					string * tempseq = new string(((*inputSequences)[mymemberlist->getFirst()].second).c_str());
					cout << string(mymemberlist->getSecond()+20, '_') << *tempseq << endl;
					delete tempseq;
				}
				cout << *consensus << endl;
				//exit(0);
			}
			
			delete consensus;
			if (not consensus_is_read) {
				(*cluster_consensus_sequences)[current_cluster] = NULL;
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
	
	//cout << countOfOverlappingKmers->name << endl;
	delete countOfOverlappingKmers;
	//cout << lastreadseen->name << endl;
	delete lastreadseen;
	
	//cout << cluster_consensus_sequences->name << endl;
	cluster_consensus_sequences->deleteContentPointers();
	//cout << "huhu" << endl;
	delete cluster_consensus_sequences;
	
	//cout << cluster_member_lists->name << endl;
	cluster_member_lists->deleteContentPointers();

	delete cluster_member_lists;
	
	delete vector_of_offsets;
	
	delete alignment_column;
	//cout << "------------- "<< " X4 " << *(inputSequences->at(0).second) << endl;
	//cout << inputSequences->name << endl;
	//inputSequences->deleteContentPairOfPointers();
	delete inputSequences;
	
	
	for (cluster_kmer_hash_it = cluster_kmer_hash.begin(); cluster_kmer_hash_it != cluster_kmer_hash.end(); ++cluster_kmer_hash_it) {
		delete cluster_kmer_hash_it->second;
		cluster_kmer_hash_it->second = NULL;
	}
	
	
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
	
		("output",				po::value< string >(),													"output file")
		("name",				po::value< string >(&prefixname)->default_value("cluster"),				"fasta description prefix name")
		("list",																						"list all members and their offsets of a cluster")
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
	
	log_stream << "clustgun git version: " << GIT_REF << endl;
	
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
	
	
	Clustgun my_pc = Clustgun();
	
	my_pc.listfile = string(listfile);
	
	if (vm.count("output")) {
		my_pc.outputfile=vm["output"].as< string >();
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
	
	
	if (vm.count("list")) {
		my_pc.list_all_members = true;
	}
	
	if (vm.count("sort")) {
		my_pc.sort_input_seq = true;
	}
	
	
	

	if (cluster_kmer_overlap_threshold < 1) {
		cerr << "error: parameter does make no sense" << endl;
		exit(1);
	}
	
	// ----------------------------------
	// run
	
	my_pc.cluster(input_file);
	
	log_stream.flush();
    log_stream.close();
		
	return 0;
}

