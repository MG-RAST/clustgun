

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
	
	
	this->code =0;
	this->coded_length = 0;
	
}

void KmerIterator::reset(int startpos){
	this->code =0;
	this->coded_length = 0;
	this->kmer_start_pos = startpos-1;
}


bool KmerIterator::nextKmer(){
	
	kmer_start_pos++;
	
	
	
	while (true) {
		if (kmer_start_pos > seqlen-kmerlength) return false;
	// normal shift by one
		if (coded_length == kmerlength) { 
			//cout << "code: " << code << endl;
			//string kmer = string_int_2_kmer(code);
			//cout << "kmer: " << kmer << endl;
			
			int front = code % aminoacid_count; 
			//cout << "front: "<< front << endl;
			//cout << "front: "<< aminoacid_int2ASCII[front] << endl;
			
			code -= front;
			code /= aminoacid_count;
			//cout << "code: " << code << endl;
			
			//int aa = aminoacid_ASCII2int[(*sequence)[kmer_start_pos+kmerlength-1]];
			
	#ifdef DEBUG
			char c = sequence[kmer_start_pos+kmerlength-1]; 
			if ((int)c < 0 ) { // that could happen if char is signed by default... I should check that...
				cerr << "(int)c < 0" << endl;
				exit(1);
			}
			int aa = aminoacid_ASCII2int[c]; 
	#else
			//int aa = aminoacid_ASCII2int[(*sequence)[kmer_start_pos+kmerlength-1]];
			int aa = aminoacid_ASCII2int[sequence[kmer_start_pos+kmerlength-1]];
	#endif
			
			if (aa == -1) {
				//cout << flush;
				//for (int i = 0; i < sequence->length(); ++i) {
				//	cout << i << ": " << sequence->at(i) << " " << (int) (sequence->at(i)) << endl;
				//}
				if (sequence[kmer_start_pos+kmerlength-1] != 'X') {
					cerr << "warning A:  amino acid not accepted \"" << sequence[kmer_start_pos+kmerlength-1] << "\" pos: " << kmer_start_pos+kmerlength-1 << endl;
					cerr << "seq: " << sequence << endl;
				}
				//std::exit(1);
				// and set kmer_start_pos XXXXXXXXXXXXX
				kmer_start_pos = kmer_start_pos+kmerlength;
				coded_length = 0; // use this to skip  characters
				this->code =0;
				
				// continue while loop
				
			} else {
			
			//cout << "new: " << (*sequence)[kmer_start_pos+kmerlength-1] << endl;
			
			
				code +=  aa * (int) pow((double)aminoacid_count , (double)kmerlength-1); 
				return true;
			}
			//cout << *(this->sequence)<< endl;
			//cout << "code: " << code << endl;
			//string kmer = string_int_2_kmer(code, kmerlength, aminoacid_int2ASCII, aminoacid_count);
			//cout << "new kmer: " << kmer << endl;
			
			//exit(0);
			
			
		} else { //(coded_length < kmerlength) {
			
			//cout << "kmer_start_pos: " << kmer_start_pos << endl;
			//cout << "seqlen-kmerlength: " << seqlen-kmerlength << endl;
			//cout << "read: " << (*sequence)[kmer_start_pos+coded_length] << endl;
			
			
			
			
			
			int aa = aminoacid_ASCII2int[sequence[kmer_start_pos+coded_length]];
			
			//cout << "aa: " << aa << endl;
			if (aa == -1) {
				if (sequence[kmer_start_pos+coded_length] != 'X') {
					cerr << "warning B: amino acid not accepted \"" << sequence[kmer_start_pos+coded_length] << "\"" << endl;
					cerr << "seq: " << sequence << endl;
				}
				//exit(1);
				// and set kmer_start_pos XXXXXXXXXXXXX
				kmer_start_pos = kmer_start_pos+coded_length;
				coded_length = 0;
				this->code =0;
			}
			//cout << "add: " << aa * pow((double)aminoacid_count , (double)coded_length) << endl;
			
			//cout << "aa: " << aa << endl;
			// cout << "mulitply: " << pow((double)aminoacid_count , (double)coded_length) << endl;
			code+=  aa * (int) pow((double)aminoacid_count , (double)coded_length);
			// cout << "code: " << code << endl;
			coded_length++;
			
			if (coded_length == kmerlength) {
				return true;
			}
			
			
		}
	}
	//cout << "code produced: " << code << endl;
	//
	//cout << *(this->sequence)<< endl;
	//cout << "code: " << code << endl;
	//string kmer = string_int_2_kmer(code, kmerlength, aminoacid_int2ASCII, aminoacid_count);
	//cout << "new kmer: " << kmer << endl;
	cerr << "error: should not reached here" << endl;
	exit(1);
	
	return false;
}

