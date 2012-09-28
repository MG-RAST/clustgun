

#include "read_blosum.hpp"

short int getScore(char a, char b, short int * blosum_matrix) {
		
	return blosum_matrix[a + 256*b];
	
}

short int getScoreSave(char a, char b, short int * blosum_matrix) {
	if (a < 0) {
		cerr << "error: a < 0" << endl;
		exit(1);
	}
	if (b < 0) {
		cerr << "error: b < 0" << endl;
		exit(1);
	}
	
	return blosum_matrix[a + 256*b];
	
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
