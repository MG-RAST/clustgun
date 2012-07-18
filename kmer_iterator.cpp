

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



KmerIterator::KmerIterator(string* seq, int startpos, int kmerlength, const char* aminoacid_int2ASCII, aminoacid* aminoacid_ASCII2int, int aminoacid_count){
	
	this->sequence = seq;
	
	
	
	
	this->kmer_start_pos = startpos-1;
	this->kmerlength = kmerlength;
	this->aminoacid_int2ASCII = aminoacid_int2ASCII;
	this->aminoacid_ASCII2int = aminoacid_ASCII2int;
	this->aminoacid_count = aminoacid_count;
	
	this->seqlen = sequence->length();
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
		int aa = aminoacid_ASCII2int[sequence->at(kmer_start_pos+kmerlength-1)]; // with boundary check on sequence
#else
		int aa = aminoacid_ASCII2int[(*sequence)[kmer_start_pos+kmerlength-1]];
#endif
		
		
		if (aa == -1) {
			cerr << "error A:  amino acid not accepted \"" << (*sequence)[kmer_start_pos+kmerlength-1] << "\" pos: " << kmer_start_pos+kmerlength-1 << endl;
			cerr << "seq: " << *sequence << endl;
			for (int i = 0; i < sequence->length(); ++i) {
				cout << i << ": " << sequence->at(i) << " " << (int) (sequence->at(i)) << endl;
			}
			std::exit(1);
			// and set kmer_start_pos XXXXXXXXXXXXX
			coded_length = 0;
		}
		
		//cout << "new: " << (*sequence)[kmer_start_pos+kmerlength-1] << endl;
		
		
		code +=  aa * (int) pow((double)aminoacid_count , (double)kmerlength-1); 
		
		//cout << *(this->sequence)<< endl;
		//cout << "code: " << code << endl;
		//string kmer = string_int_2_kmer(code, kmerlength, aminoacid_int2ASCII, aminoacid_count);
		//cout << "new kmer: " << kmer << endl;
		
		//exit(0);
		
		return true;
	}
	
	// build code
	while (coded_length < kmerlength) {
		
		//cout << "kmer_start_pos: " << kmer_start_pos << endl;
		//cout << "seqlen-kmerlength: " << seqlen-kmerlength << endl;
		//cout << "read: " << (*sequence)[kmer_start_pos+coded_length] << endl;
		
		
		
		
#ifdef DEBUG
		int aa = aminoacid_ASCII2int[sequence->at(kmer_start_pos+coded_length)];
#else
		int aa = aminoacid_ASCII2int[(*sequence)[kmer_start_pos+coded_length]];
#endif
		
		//cout << "aa: " << aa << endl;
		if (aa == -1) {
			cerr << "error B: amino acid not accepted \"" << (*sequence)[kmer_start_pos+coded_length] << "\"" << endl;
			exit(1);
			// and set kmer_start_pos XXXXXXXXXXXXX
			coded_length = 0;
		}
		//cout << "add: " << aa * pow((double)aminoacid_count , (double)coded_length) << endl;
		
		//cout << "aa: " << aa << endl;
		// cout << "mulitply: " << pow((double)aminoacid_count , (double)coded_length) << endl;
		code+=  aa * (int) pow((double)aminoacid_count , (double)coded_length);
		// cout << "code: " << code << endl;
		coded_length++;
		
	}
	//cout << "code produced: " << code << endl;
	//
	//cout << *(this->sequence)<< endl;
	//cout << "code: " << code << endl;
	//string kmer = string_int_2_kmer(code, kmerlength, aminoacid_int2ASCII, aminoacid_count);
	//cout << "new kmer: " << kmer << endl;

	
	return true;
}

