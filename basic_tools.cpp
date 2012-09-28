

#include "basic_tools.hpp"


double round( double value )
{
	return floor( value + 0.5 );
}

string stream2string(std::stringstream & stream) {
	std::string s;
	std::ostringstream os;
	os<<stream.rdbuf();
	s=os.str();
	return s;
}

bool FileExists(const char * filename) {
	struct stat stFileInfo;
	bool blnReturn;
	int intStat;
	
	// Attempt to get the file attributes
	intStat = stat(filename,&stFileInfo);
	if(intStat == 0) {
		// We were able to get the file attributes
		// so the file obviously exists.
		blnReturn = true;
	} else {
		// We were not able to get the file attributes.
		// This may mean that we don't have permission to
		// access the folder which contains this file. If you
		// need to do that level of checking, lookup the
		// return values of stat which will give you
		// more details on why stat failed.
		blnReturn = false;
	}
	
	return(blnReturn);
}


void checkFile(string * file) {
	if (file == 0) {
		return;
	}
	
	if (! FileExists((char *) file->c_str())) {
		cerr << "error reading file: " << *file << endl;
		exit(EXIT_FAILURE);
	}
	
}

string int2str (int n) {
	stringstream ss;
	ss << n;
	return ss.str();
}

int str2int(string& text){
	
	int number = std::atoi( text.c_str() );
	
	return number;
	
}

void reverseComplementDNA(string& sequence) {
	int len = sequence.length();
	
	// A-T
	// G-C
	
	// reverse:
	int half = len/2;
	for (int i=0; i<half; i++) {
		swap(sequence[i], sequence[len-1-i]);
	}
	
	// complement:
	char c;
	for (int i=0; i<len; i++) {
		c = sequence[i];
		toupper(c);
		switch (c) {
			case 'A':
				c='T';
				break;
			case 'C':
				c='G';
				break;
			case 'G':
				c='C';
				break;
			case 'T':
				c='A';
				break;
			case 'U':
				c='C';
				break;
			default:
				//cerr << "error at " << i << " with " << sequence << endl;
				//exit(1);
				c='N';
		}
		sequence[i]=c;
	}
	
	
}
