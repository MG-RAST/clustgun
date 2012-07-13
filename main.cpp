//
//  main.cpp
//  clustgun
//
//  Created by Wolfgang Gerlach on 5/9/12.
//  Copyright (c) 2012 University of Chicago. All rights reserved.
//


#include "main.h"

int str2int(string& text){
	
	int number = std::atoi( text.c_str() );
	
	return number;
	
}


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





short int * readBLOSUM(const char * blosum_file) {
	
	string line;
	
	int matrix_length = 256*256;
	short int * matrix = new short int [matrix_length];
	
	for (int i = 0 ; i < matrix_length; i++) {
		matrix[i] = SHRT_MAX;
	}
	
	char * row = 0 ;
	
	ifstream myfile (blosum_file);
	if (myfile.is_open())
	{
		while (! myfile.eof() )	{
			char * pch;
			getline (myfile,line, '\n');
			
			
			//cout << line << endl;
			
			
			
			if (line.length() == 0) {
				continue;
			}
			
			if (line[0] == '#') {
				continue;
			}
			
			
			char buf[1024];
			strncpy(buf, line.c_str(), sizeof(buf) - 1);
			buf[sizeof(buf) - 1] = '\0';
			
			if (line[0] == ' ') {
				row = new char [256];
				
				pch = strtok (buf," ");
				int i = 0;
				while (pch != NULL)
				{
					//cout << i  << ": " <<  pch[0] << endl;
					char c = pch[0];
					
					//if (c<0 || c > 256) {
					if (c<0) {
						cerr << "c invalid" << endl;
						exit(EXIT_FAILURE);
					}
					
					row[i]=c;
					//printf ("%s\n",pch);
					pch = strtok (NULL, " ");
					i++;
				}
				
				
				continue;
			}
			
			if (row == 0) {
				cerr << "row=0" << endl;
				exit(EXIT_FAILURE);
			}
			
			// normal rows:
			pch = strtok (buf," ");
			
			char c = pch[0];
			
			//if (c<0 || c > 256) {
			if (c<0) {
				cerr << "c invalid" << endl;
				exit(EXIT_FAILURE);
			}
			
			pch = strtok (NULL, " ");
			
			int i = 0;
			while (pch != NULL) {
				string score_str(pch);
				short int score_int = (short int) str2int(score_str);
				
				//cout << c << "," << row[i] << ": " << score_int << endl;
				matrix[c+256*row[i]]=score_int;
				matrix[c+32+256*row[i]]=score_int;
				matrix[c+256*(row[i]+32)]=score_int;
				matrix[c+32+256*(row[i]+32)]=score_int;
				pch = strtok (NULL, " ");
				i++;
			}
			
			
			
		}
		
		myfile.close();
		
		if (row != 0) delete [] row;
		
	}
	else {
		cerr << "Error: Unable to open file " << blosum_file << endl;
		exit(1);
	}
	
	return matrix;
}



int getOffsets(string * sequence, vector<short> * vector_of_offsets, sparse_hash_map<int, kmer_appearance_list * , hash<int>, eqint> &cluster_kmer_hash, int max_cluster) {
	
	sparse_hash_map<int, kmer_appearance_list * , hash<int>, eqint>::iterator cluster_kmer_hash_it;
	KmerIterator * mykmerit = new KmerIterator(sequence,0, kmerlength, aminoacid_int2ASCII, aminoacid_ASCII2int, aminoacid_count);
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
					
					(*vector_of_offsets)[number_of_overlapping_kmers_seen-1]=clusterpos-readpos;
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
						  HashedArrayTree<pair<string *, string * > > * inputSequences,
						  vector<char> * alignment_column,
						  vector<bool> * aminoacid_occurence,
						  vector<int> * aminoacid_scores,
						  short * score_matrix) {
		
	vector<pair<int, short> > sorted_members;
	
	mymemberlist->resetIterator();
	while (mymemberlist->nextElement()){
		sorted_members.push_back(pair<int, short>(mymemberlist->getFirst(), mymemberlist->getSecond()));
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
	string * new_consensus_seq = new string();

	// update offsets (leftmost member has now offset zero):
	for (int member = 0; member < sorted_members.size(); ++member) {
		//cout << "member: " << member << " offset(org): " << sorted_members[member].second << endl;
		int member_offset = sorted_members[member].second - base_offset;
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
			
			if (sorted_members[member].second + ((*inputSequences)[sorted_members[member].first].second)->length() - 1 < current_ali_position ){
				if (!saw_aa) {
					first_member = member + 1;
				}
				
			} else {
				saw_aa = true;
				//cout << member << endl;
				char aa = (*((*inputSequences)[sorted_members[member].first].second))[current_ali_position-sorted_members[member].second];
				
				//cout << aa << endl;
				
				last_alignment_char++;
				(*alignment_column)[last_alignment_char] = aa;
			}
			//cout << string(sorted_members[member].second, '_') << *((*inputSequences)[sorted_members[member].first].second)  << endl;
		}
		//exit(0);
		
		//find consensus aa:
		if (last_alignment_char < 0 ) {
			break;
		}
		
		pair<bool, char > major_aa =  majority_vote<char>(alignment_column, last_alignment_char+1);
		
		#ifdef DEBUG
		char test_char = major_aa.second;
		#endif
		
		if ((int)major_aa.second == 0) {
			cerr << "B    (int)major_aa.second == 0" << endl;
			exit(1);
			
		}
		if ((int)major_aa.second == 0) {
			cerr << "BB    (int)major_aa.second == 0" << endl;
			exit(1);
			
		}
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
				aminoacid aa = aminoacid_ASCII2int[(*alignment_column)[i]];
				(*aminoacid_occurence)[aa] = true;
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
				if ((*aminoacid_occurence)[aa]) {
					// walk through cloumn and compute score
					for (int i = 0; i <= last_alignment_char; ++i) {
						char aa_i = (*alignment_column)[i];
#ifdef DEBUG
						if (aminoacid_int2ASCII[aa]+256*aa_i > 256*256) {
							cerr << "error: aminoacid_int2ASCII[aa]+256*aa_i > 256*256" << endl;
						}
#endif
						
						(*aminoacid_scores)[aa] += score_matrix[aminoacid_int2ASCII[aa]+256*aa_i];
					}
					//cout << aminoacid_int2ASCII[aa] << ": "<< (*aminoacid_scores)[aa] << endl;
					if ((*aminoacid_scores)[aa] > max_aa_score) {
						max_aa_score = (*aminoacid_scores)[aa];
						max_aa = aa;
					}
					
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
				cout << i << ": " << (*alignment_column)[i] << endl;
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

void addConsensusSequence(string * new_consensus_seq, int cluster, sparse_hash_map<int, kmer_appearance_list * , hash<int>, eqint>& cluster_kmer_hash) {

	
	sparse_hash_map<int, kmer_appearance_list * , hash<int>, eqint>::iterator cluster_kmer_hash_it;
	
	KmerIterator * mykmer_insertion_it = new KmerIterator(new_consensus_seq, 0, kmerlength, aminoacid_int2ASCII, aminoacid_ASCII2int, aminoacid_count);
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
			mylist = new kmer_appearance_list();
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
	
	
	KmerIterator * mykmer_deletion_it = new KmerIterator(consensus_old, 0, kmerlength, aminoacid_int2ASCII, aminoacid_ASCII2int, aminoacid_count);
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


bool computeSequenceOverlap(int offset, string * a, string * b, short * score_matrix) {
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
	
	int min_required_length = (int) ((double)min(a_len, b_len)*(double)min_overlap_fraction);
	
	if (overlap_length < min_required_length) {
		//cout << "failed ==============================" << endl;
		return false;
	}
	
	
	
	
	//exit(0);
	int total_score = 0;
	int tot_len = 0;
	
	int aa_score;
	
	int currentWinLen = 0;
	int windowScore = 0;
	int index;
	
	int start_i = i;
	int end_i = i+overlap_length;
	for (int i = start_i; i< end_i ; ++i) {
		tot_len++;
		j = i - offset;
		
		//cout << "i " << (*a)[i] << " "  << (*b)[j] << " " << i << " " << j  << endl;
		
		
		
		
		index = (*a)[i]+256*(*b)[j];
		
		aa_score = score_matrix[index];
		
		total_score += aa_score;
		windowScore += aa_score;
		
		if (currentWinLen < windowLength) {
			//cout << "_windows score is : " << windowScore << endl;
			currentWinLen++;
		} else {
			
			//remove first aa score, of aa that has left window
			
			index = (*a)[i-windowLength]+256*(*b)[j-windowLength];
			
			aa_score = score_matrix[index];
			windowScore -= aa_score;
			
			
			//cout << "windows score is : " << windowScore << endl;
			if (windowScore < windowScoreThreshold) {
				//cerr << "windowScore is too small" << endl;
				return false;
				//exit(1);
				
			}
		}
				
				
		
		//cout << "score: " << score << endl;
		
		
		
	}
	
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


		
  
	std::cout << "Hello, World!\n";

	system("date");
	time_t begin, end; 
	time(&begin);
	
	
	// create mappings of amino acid ASCII to int and back.
	
	cout << "initialize alphabet mapping... \"" << aminoacid_int2ASCII << "\""<< endl;
	for (aminoacid i = 0 ; i<20; ++i) {
		char aa = aminoacid_int2ASCII[i];
		//cout << aa << " " << i << endl;
		aminoacid_ASCII2int[aa]=i;
		aa = tolower ( aa );
		aminoacid_ASCII2int[aa]=i;
		//cout << aa << " " << i << endl;
	}
  
	
	short * score_matrix = readBLOSUM(blosum_file);
	
	
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
	string * sequence;

	// ----------------------------------------
	// read sequences into memory
	HashedArrayTree<pair<string *, string * > > * inputSequences 
				= new HashedArrayTree<pair<string *, string * > >(20); // 20 for 2^20=1MB chunks
	
	#ifdef DEBUG
	inputSequences->name = string("inputSequences");
	#endif
	
	fasta_parser = new FASTA_Parser(file, true, zcat_command);
	//string my_sequence;
	while (fasta_parser->getNextDescriptionLine(descr)) {
		
		//cerr << "descr: " << descr<< endl;
		sequence = fasta_parser->getSequence();
		//my_sequence = *sequence;
		
		
		
		if (descr.length() > 0 && sequence->length() >= minimal_input_sequence_length) {
			//cout << "huhu: " << *sequence << endl;
			total_read_count++;
			inputSequences->push_back(pair<string *, string * >(new string(descr), sequence));
			//cout << "pushed: " << inputSequences->size() << endl;	
			//cout << *sequence << endl;
			//exit(0);
							
			   
				
						   
			
		} else {
			cerr << "warning: sequence was not accepted..." << endl;
			delete sequence;
		}
		
		if (total_read_count % 100000 == 0) cerr << "total_read_count: " << total_read_count << endl;
		//if (total_read_count >= limit_input_reads && limit_input_reads != -1) {
		//	cerr << "WARNING: numer of reads for debugging purposes limited !!!!!!" << endl;
		//	break;
		//}
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
	
	if (sort_input_seq) {
		sortInputSequences(inputSequences);
	}
	
	
	
	
	
	bool kmerabundance = false;
	
	if (kmerabundance) {
		vector< pair<int, int> > * array;
		vector< pair<int, int> >::iterator array_it;
		
		int sum;
		int totsum;
		
		
		array = new vector< pair<int, int> >;
		
		// print hash table / put in array
		for ( it=kmer_count_hash.begin() ; it != kmer_count_hash.end(); it++ ) {
			if ((*it).second > low_abundance_threshold) {
				// cout << string_int_2_kmer((*it).first) << " => " << (*it).second << endl;
				array->push_back(pair<int, int >((*it).first, (*it).second));
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
	HashedArrayTree<short> * countOfOverlappingKmers = new HashedArrayTree<short>(20, 0);
#ifdef DEBUG
	countOfOverlappingKmers->name = string("countOfOverlappingKmers");
#endif
	HashedArrayTree<int> * lastreadseen = new HashedArrayTree<int>(20, -1); // this avoids initialization
#ifdef DEBUG
	lastreadseen->name = string("lastreadseen");
#endif
	HashedArrayTree<string * > * cluster_consensus_sequences = new HashedArrayTree<string * >(20, NULL);
#ifdef DEBUG
	cluster_consensus_sequences->name = string("cluster_consensus_sequences");
#endif	
	//HashedArrayTree<short> * cluster_consensus_offsets = new HashedArrayTree<short>(20, 0); // offset relative to seed read
	
	
	
	HashedArrayTree<cluster_member_list * > * cluster_member_lists = new HashedArrayTree<cluster_member_list * >(20, NULL);

#ifdef DEBUG
	cluster_member_lists->name = string("cluster_member_lists");
#endif
	
	vector<short> * vector_of_offsets = new vector<short>();
	vector_of_offsets->reserve(max_protein_length);
	
	vector<char> * alignment_column = new vector<char>();
	alignment_column->reserve(max_protein_length);
	
	vector<bool> * aminoacid_occurence = new vector<bool>(aminoacid_count, false);
	vector<int> * aminoacid_scores = new vector<int>(aminoacid_count, 0);
	
	int last_cluster=-1;

	
	
	
	countOfOverlappingKmers->reserve(hat_increase_steps);
	lastreadseen->reserve(hat_increase_steps);
	cluster_consensus_sequences->reserve(hat_increase_steps);
	cluster_member_lists->reserve(hat_increase_steps);
	
	// -------------------------------------------------------------------------
	
	//iterate through all reads to create clusters and match against existing clusters
	cout << "read count vs cluster count" << endl;
	
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

	
	for (size_t read_id = 0 ; read_id < sequence_count; read_id++){
	
		//cout << "read_id: " << read_id<< endl;
	
		sequence = (*inputSequences)[read_id].second;
		
					
		KmerIterator * mykmerit = new KmerIterator(sequence, 0, kmerlength, aminoacid_int2ASCII, aminoacid_ASCII2int, aminoacid_count);
		
		//int max_cluster = -1;
		//int max_cluster_kmer_count = 0;
		
		
		
		for (int i = 0; i < top_n_clusters; ++i) {
			max_cluster_id_array[i]=-1;
			max_cluster_kmer_count_array[i]=0;
			
		}
		
		//cout << "------------- "<< read_id << " " << *(inputSequences->at(0).second) << endl;

		
		// check if read has matches to known clusters, iterate through read kmers
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
					//cerr << "clusterhit: " << clusterhit << endl;
					if ((*lastreadseen)[clusterhit] == read_id) {
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
						(*lastreadseen)[clusterhit] = read_id;
						(*countOfOverlappingKmers)[clusterhit]=1;
						//cout << "not seen before" << endl;
					}
				} // end while
				
			} 				
			
			
		} // end while ( kmer iteration )
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
		
		
		int validated_overlaps=0;
		
		// validate overlaps
		if (max_cluster_id_array[0] >= 0) { // number of k-mer-hits indicates a hit, needs to be verified.
			// (probably) found cluster for the read !
			
			//int seed_sequence_number =  (*cluster_seedread)[max_cluster];
			
			if (false) {
				cout << "yeah --------------------" << endl;
				cout << "read: " << read_id << endl;
				cout << "readseq: " << *sequence << endl;
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
				int number_of_overlapping_kmers_seen = getOffsets(sequence,
																  vector_of_offsets,	// <-- output
																  cluster_kmer_hash,
																  max_cluster_id_array[cluster_it]);
				
				// majority vote on offsets				
				pair<bool, short> major = majority_vote<short>(vector_of_offsets, number_of_overlapping_kmers_seen);
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
					string * consensus_seq = (*cluster_consensus_sequences)[cluster];
					if (consensus_seq == NULL) {
						cluster_member_list * cluster_memberlist = (*cluster_member_lists)[cluster];
						int cluster_read = cluster_memberlist->start->data_array1[0];
						consensus_seq = inputSequences->at(cluster_read).second;
						
					}
					
					
					overlap_found = computeSequenceOverlap(majority_offset, consensus_seq, sequence, score_matrix);
					
				
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
		
		
		
		
		
		// ************ process overlaps ************
		//cout << "validated_overlaps"  << endl;
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
				
				cluster_member_list*& mymemberlist = (*cluster_member_lists)[last_cluster]; // alias to a pointer
				
				//cout << mymemberlist << endl;
				assert( (mymemberlist == NULL) );
					
				
				mymemberlist = new cluster_member_list;
				mymemberlist->append(read_id, 0);
				
				
				
				addConsensusSequence(sequence, last_cluster, cluster_kmer_hash);
				
				break;
				
			} else if (validated_overlaps == 1) { // read has match to exactly one cluster
				
				int clusterhit = max_cluster_id_array[0];
				
				
				
				cluster_member_list*& mymemberlist = (*cluster_member_lists)[clusterhit]; // alias to a pointer
				
				if (mymemberlist == NULL){
					mymemberlist = new cluster_member_list();
				}
				
				mymemberlist->append(read_id, max_cluster_offsets[0]);
				
				
				// find old consensus sequence
				string * consensus_old = (*cluster_consensus_sequences)[clusterhit];
				
				if (consensus_old == NULL) {
					// consensus seuquence does not exist yet, seed member had been used
					int seed_id = mymemberlist->start->data_array1[0];
					consensus_old = (*inputSequences)[seed_id].second;
				}
				
				
				
				// check if new member is a perfect full-length overlap
				bool extends_consensus = true;
				if (true) {
					int member_len = (*inputSequences)[read_id].second->length();
					
					//bool perfect_overlap = false;
					if (max_cluster_offsets[0] >= 0 && max_cluster_offsets[0]+member_len <= consensus_old->length() ) { // check for full-length overlap
						extends_consensus = false;
						//string * newMember = (*inputSequences)[read_id].second;
						//perfect_overlap = true;
						//for (int i = 0 ; i < member_len; ++i ) {
					//	if ((*newMember)[i] != (*consensus_old)[max_cluster_offsets[0] + i]) {
					//			//cout << i << " "  << (*newMember)[i] << " " << (*consensus_old)[max_cluster_offsets[0] + i] << endl;
					//			perfect_overlap = false;
					//			break;
					//		}
					//		
					//	}
		
						
						//exit(0);
					} 
				}
				
				int membercount = mymemberlist->getLength();
				//cout << "membercount: " << membercount << endl;
				
				
				if (clusterNeedsUpdate(membercount) || extends_consensus) { // mymemberlist->getLength() >= 100
				//if (true) {	
					
					
					
					// ============ REMOVE OLD CONSENSUS K-MERS =============
					
					
					
					
					// remove consensus sequence kmers for cluster "clusterhit"
					removeConsensusSequences(consensus_old, clusterhit, cluster_kmer_hash);
					
					
					// delete old consensus, but not seed !
					if ((*cluster_consensus_sequences)[clusterhit] != NULL) {
						delete consensus_old;
						(*cluster_consensus_sequences)[clusterhit] = NULL;
					}
					//consensus_old = NULL;					
					
					
					// ============= SORT, UPDATE OFFSETS AND COMPUTE CONSENSUS ==============
					// copy members from list to array (temporarly)
					
					string * new_consensus_seq = computeConsensus(mymemberlist,  
																  inputSequences, 
																  alignment_column,
																  aminoacid_occurence,
																  aminoacid_scores,
																  score_matrix);
					
					(*cluster_consensus_sequences)[clusterhit] = new_consensus_seq;
					//cout << "write " << new_consensus_seq << " at " << clusterhit << endl;
					//cluster_consensus_sequences->
					
					//cout <<"con: " << *new_consensus_seq << endl;
					
					
					
					
					// insert kmers of new_consensus_seq into hash
					
					addConsensusSequence(new_consensus_seq, clusterhit, cluster_kmer_hash);
					
					
					
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
					consensus_1 = inputSequences->at(cluster_read).second;
				}
				
				// Cluster 2
				int cluster_2 = max_cluster_id_array[second_cluster_index];
				cluster_member_list * cluster_2_memberlist = (*cluster_member_lists)[cluster_2];
				string * consensus_2 = (*cluster_consensus_sequences)[cluster_2];
				
				
				if (consensus_2 == NULL) {
					int cluster_read = cluster_2_memberlist->start->data_array1[0];
					consensus_2 = inputSequences->at(cluster_read).second;
				}
				
				//cout << cluster_1 << endl;
				//cout << cluster_2 << endl;
				
				//cout << max_cluster_offsets[previous_master_cluster_index] << endl;
				//cout << max_cluster_offsets[second_cluster_index] << endl;
				
				int offset_diff = max_cluster_offsets[previous_master_cluster_index] - max_cluster_offsets[second_cluster_index];
				//cout << "seq: " << *sequence << endl;
				//cout << "o: " << offset_diff << endl;
				
				bool has_overlap = computeSequenceOverlap(offset_diff, consensus_1, consensus_2, score_matrix);
				
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
					
					if ( (*cluster_consensus_sequences)[ cluster_1] != NULL ) {
						delete (*cluster_consensus_sequences)[ cluster_1];
						(*cluster_consensus_sequences)[ cluster_1] = NULL;
						
					}
					if ( (*cluster_consensus_sequences)[ cluster_2] != NULL ) {
						delete (*cluster_consensus_sequences)[ cluster_2];
						(*cluster_consensus_sequences)[ cluster_2] = NULL;
						
					}
					
					
					
					
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
					
					(*cluster_consensus_sequences)[master_cluster] = new_consensus_seq;
					//cout << "write2 " << new_consensus_seq << " at " << master_cluster << endl;
					//cout <<"merged con: " << *new_consensus_seq << endl;
					//exit(0);
					//cout << "E" << endl;
					
					
					// insert kmers of new_consensus_seq into hash
					
					addConsensusSequence(new_consensus_seq, master_cluster, cluster_kmer_hash);
					
					previous_master_cluster_index = master_cluster_index;
					stat_real_cluster_count--;
					
					break;				
				}
				
				
				
				
				
				
				//cout << "leave merge process " << endl;
				//exit(0);
			} 		
			
			
		  
		} // end while true	
			
		
		if ((read_id +1) % 10000 == 0 && read_id > 0) {
			//cerr << "count: " << read_id << endl;
			#ifdef TIME		
			time(&end);
			#endif		
			
			cout << read_id+1 << "\t" << stat_real_cluster_count 
			#ifdef TIME
			<< "\t" << difftime(end, begin)
			#endif			
			<< endl;
		}
		//if ((read_id >= limit_input_reads) && (limit_input_reads != -1)) {
		//	break;
		//}
	}
	
	if (sequence_count % 10000 != 0 && sequence_count > 0) {
		//cerr << "count: " << read_id << endl;
		cout << sequence_count << "\t" << stat_real_cluster_count  << endl;
	}
	
	delete fasta_parser;
	//cout << "------------- "<< " X1 " << *(inputSequences->at(0).second) << endl;

	
	
	// --------------------------------------------------------------------------------------------------
	// write consensus sequences to output file
	string outputfile = string(inputfile)+string(".consensus");
	cerr << " done. write results into output file \"" << outputfile << "\"... " << endl;

	

	
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
				
				
				consensus = (*inputSequences)[mymemberlist->start->data_array1[0]].second ;
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
			
			
			
			
			output_stream << "> cluster" << current_cluster << " length=" << consensus->length() << " size=" << (*cluster_member_lists)[current_cluster]->getLength();
			
			if (list_all_members) {
				output_stream << " ";
				
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
						output_stream << " " ; 
					}
					
					output_stream << "(" << offset << ")" << *((*inputSequences)[seq_id].first);
				}
				
				
			} 
			
			output_stream << endl;
			
			output_stream << *consensus << endl;
			
			
			
			//if (current_cluster == 98240) {
			//if(true){
			
			//if (mymemberlist->getLength() > 10000) {
			//if (current_cluster == 0) {
			if (false) {
				cluster_member_list*& mymemberlist = (*cluster_member_lists)[current_cluster];
				mymemberlist->resetIterator();
				while(mymemberlist->nextElement()) {
					string * tempseq = ((*inputSequences)[mymemberlist->getFirst()].second);
					cout << string(mymemberlist->getSecond()+20, '_') << *tempseq << endl;
				}
				cout << *consensus << endl;
				//exit(0);
			}
			
			
			//if (! is_read) {
			//	delete consensus;
			//	(*cluster_consensus_sequences)[current_cluster] = NULL;
			//}
			//cout << "> cluster" << current_cluster << endl;
			//cout << *consensus << endl;
		}
		
	}
	//cout << "------------- "<< " X3 " << *(inputSequences->at(0).second) << endl;
	output_stream.close();
	
	

	
	
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
	inputSequences->deleteContentPairOfPointers();
	delete inputSequences;
	
	
	for (cluster_kmer_hash_it = cluster_kmer_hash.begin(); cluster_kmer_hash_it != cluster_kmer_hash.end(); ++cluster_kmer_hash_it) {
		
		delete cluster_kmer_hash_it->second;
	}
	
	
	delete [] score_matrix;
	
	
	
	delete zcat_command;
	delete aminoacid_occurence;
	delete aminoacid_scores;
		
	system("date");
	cerr << "end" << endl;
	
}

void usage(boost::program_options::options_description& options) {
	cerr << endl;
	cerr << "Usage: clustgun [--option] input.fas" << endl;
	cerr << endl;
	cerr << options << endl;
	cerr << endl;
}








int main(int argc, const char * argv[])	{
	
	
	namespace po = boost::program_options;
	
	po::options_description options_visible("Options");
	options_visible.add_options()
		("sort", 	"sort input sequences by length (slow!)")
		("list", 	"list all members and their offsets of a cluster")
		("help", 	"display this information");
	
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

	
	if (vm.count("help")) {
		usage(options_visible);
		exit(0);
	}
	
	string input_file;

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
		
		
		
	} else {
		usage(options_visible);
		cerr << endl;
		cerr << "error: input file is missing..." << endl;
		exit(1);
	}

	
	
	
	//if (argc != 2) {
	//	usage(options_visible);
	//	exit(1);
	//}
	
	
	
	
	Clustgun my_pc = Clustgun();
	
	if (vm.count("list")) {
		my_pc.list_all_members = true;
	}
	
	if (vm.count("sort")) {
		my_pc.sort_input_seq = true;
	}
	
	my_pc.cluster(input_file);
	
		
	return 0;
}

