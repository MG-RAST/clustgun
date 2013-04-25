

#include "kmer_iterator.hpp"




string string_int_2_kmer(int kmer_code, int kmerlength, const char * aminoacid_int2ASCII, int aminoacid_count) {
	
	
	string kmerstring(kmerlength, '.');
	
	
	for (int i = 0 ; i < kmerlength-1; ++i){
		//cout << "kmer_code: " << kmer_code << endl; 
		//cout << "aminoacid_count: " << aminoacid_count << endl;
		//cout << "decode: " << (kmer_code % aminoacid_count) << endl;
		kmerstring[i] = aminoacid_int2ASCII[kmer_code % aminoacid_count];
		kmer_code =  kmer_code / aminoacid_count;
	}
	//cout << "decode: " << (kmer_code % aminoacid_count) << endl;
	kmerstring[kmerlength-1] = aminoacid_int2ASCII[kmer_code % aminoacid_count];
	
	//cout << "got: " <<kmerstring << endl;
	
	return kmerstring;
}



KmerIterator::KmerIterator(const char* seq, int seqlen, int startpos, int kmerlength, const char* aminoacid_int2ASCII, aminoacid* aminoacid_ASCII2int, int aminoacid_count){
	
	this->sequence = seq;
	this->seqlen = seqlen;
	
	
	this->kmer_start_pos = startpos-1;
	this->kmerlength = kmerlength;
	this->aminoacid_int2ASCII = aminoacid_int2ASCII;
	this->aminoacid_ASCII2int = aminoacid_ASCII2int;
	this->aminoacid_count = aminoacid_count;
	
	
	this->code = 0;
	this->coded_length = 0;
	
	
	for (int i = 0 ; i <= kmerlength ; ++i) {
		powertable[i] = ipow( aminoacid_count , i);
	}
	
	#ifdef DEBUG
	if (strlen(seq) != (size_t)seqlen) {
		cerr <<"error: strlen(seq) != seqlen" << endl;
		cerr << strlen(seq) << endl;
		cerr << seqlen << endl;
		exit(1);
	}
	#endif
	
}

void KmerIterator::reset(int startpos){
	this->code = 0;
	this->coded_length = 0;
	this->kmer_start_pos = startpos-1;
}


bool KmerIterator::nextKmer(){
	
	kmer_start_pos++;
	
	if (coded_length == kmerlength) {
		//cout << "code: " << code << endl;
		//string kmer = string_int_2_kmer(code);
		//cout << "kmer: " << kmer << endl;
		
		int front = code % aminoacid_count;
		//cout << "front: "<< front << endl;
		//cout << "front: "<< aminoacid_int2ASCII[front] << endl;
		
		code -= front;
		code /= aminoacid_count;
		coded_length--;
	}
	
	while (coded_length < kmerlength) {
		if (kmer_start_pos > seqlen-kmerlength) return false;
		
		int new_character_position = kmer_start_pos+coded_length;
		
		#ifdef DEBUG
		char c = sequence[new_character_position]; 
		if ((int)c < 0 ) { // that could happen if char is signed by default... I should check that...
			cerr << "(int)c < 0" << endl;
			exit(1);
		}
		int aa = aminoacid_ASCII2int[c]; 
		#else
		//int aa = aminoacid_ASCII2int[(*sequence)[kmer_start_pos+kmerlength-1]];
		int aa = aminoacid_ASCII2int[(int) sequence[new_character_position]];
		#endif
		
		if (aa == -1) {
			//cout << flush;
			//for (int i = 0; i < sequence->length(); ++i) {
			//	cout << i << ": " << sequence->at(i) << " " << (int) (sequence->at(i)) << endl;
			//}
			if (sequence[new_character_position] != 'X') {
				cerr << "warning A:  amino acid not accepted \"" << sequence[new_character_position] << "\" code:" << (int) sequence[new_character_position] << " pos: " << new_character_position << endl;
				cerr << "seq: " << sequence << endl;
				
				//#ifdef DEBUG
				exit(1);
				//#endif
			}
			//std::exit(1);
			// and set kmer_start_pos XXXXXXXXXXXXX
			kmer_start_pos = new_character_position+1;
			coded_length = 0; // use this to skip character
			this->code =0;
			
			
		} else {
			code +=  aa * powertable[coded_length]; //(int) ipow( aminoacid_count , coded_length); // comment: could have replaced multiplication with a look-up, but mulit seems to single instruction now...
			coded_length++;
		}
		
			
			
			
		
	} // end while
	
		
	return true;
}

