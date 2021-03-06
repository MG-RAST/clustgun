



#include "fasta_parser.hpp"



string getFASTADescription (string & line) {
	size_t start = line.find_first_not_of(" \t", 1);
	
	if (start == string::npos)  {
		cerr << "error: no char found" << endl;
		exit(1);
	}
	
	size_t end = line.find_first_of(" \t\n", start);
	
	string identifier;
	
	if (end == string::npos)  {
		identifier = line.substr(start);
	} else {
		identifier = line.substr(start, end-start);
	}
	return identifier;
}



FASTA_Parser::~FASTA_Parser(){
	
	if (fp != 0) {
		pclose(fp);
		fp = 0;
	}
	delete rna_stream;
	delete fpstream;
}


FASTA_Parser::FASTA_Parser(string command, bool getOnlyIdentfier) {
	
	this->start(command, getOnlyIdentfier);
	
}

FASTA_Parser::FASTA_Parser(string file, bool getOnlyIdentfier, string * zcat_bin) {
	
	string zcat_command = *zcat_bin;
	zcat_command.append(" -f ");
	zcat_command.append(file);
	
	this->start(zcat_command, getOnlyIdentfier);
	
}

FASTA_Parser::FASTA_Parser(string file, bool getOnlyIdentfier, bool use_zcat) {
	
	#ifdef __APPLE__
		string * zcat_bin = new string("/usr/bin/gzcat"); 
	#else
		string * zcat_bin = new string("/bin/zcat"); 
	#endif

	string zcat_command = *zcat_bin;
	zcat_command.append(" -f ");
	zcat_command.append(file);
	
	this->start(zcat_command, getOnlyIdentfier);
	
	delete zcat_bin;
}

void FASTA_Parser::start(string command, bool getOnlyIdentfier) {
	
	this->getOnlyIdentfier=getOnlyIdentfier;
	
	expected_sequence_length=0;
	
	
	fp = popen(command.c_str(), "r");
	if (fp == NULL) {
		cerr << "error: " << command << endl;
		exit(EXIT_FAILURE);
	}
	
	
	namespace io = boost::iostreams;
	
	io::file_descriptor_flags flags= io::close_handle;
	
	
	
	fpstream = new io::stream_buffer<io::file_descriptor_source>(fileno(fp), flags);
	
	rna_stream = new std::istream(fpstream);
	//rna_stream = std::ifstream(fpstream);
	//rna_stream.open(fpstream, ifstream::in);
	
	string line;
	
	
	while (1) {
		
		std::getline(*rna_stream, line);
		if (rna_stream->eof()) {
			return;
		}
		
		if (line.length() == 0) {
			continue;
		}
		
		if (line[0] == '#' ) {
			continue;
		}
		break;
	}
	
	if (line.length() < 2) {
		cerr << "FASTA_Parser: error parsing" << endl;
	}
	if (line[0] != '>') {
		cerr << "FASTA_Parser: error parsing" << endl;
	}
	if (getOnlyIdentfier) {
		this->description_line = getFASTADescription(line);
	} else {
		this->description_line = line;
	}
	have_read_descr = false;
	
}


bool FASTA_Parser::getNextDescriptionLine(string& descr) {
	//cerr << "getNextDescriptionLine" << endl;
	//if ( *rna_stream == false ) {
	//std::istream & rna_stream_ref = *rna_stream;
	sequence = 0;
	
	if ( rna_stream->fail() ) {
		return false;
	}
	
	// jump over sequence to next description line
	if (have_read_descr == true) {
		bool found_descr=false;
		string line;
		while (not rna_stream->fail()) {
			
			std::getline(*rna_stream, line);
			//cerr << "x:" << line << endl;
			
			if (line.length() == 0) {
				continue;
			}
			
			if (line[0] == '#' ) {
				continue;
			}
			
			
			if (line[0] == '>' ) {
				found_descr=true;
				//exit(1);
				//this->description_line = getFASTADescription(line);
				if (getOnlyIdentfier) {
					this->description_line = getFASTADescription(line);
				} else {
					this->description_line = line;
				}
				break;
			}
			
		}
		if (! found_descr) {
			return false;
		}
		
	}
	
	descr = this->description_line;
	have_read_descr = true;
	
	return true;
}




string * FASTA_Parser::getSequence() {
	string line;
	
	
	if (sequence != 0) {
		cerr << "error: sequence has already been read via getSequence()" << endl;
		exit(1);
	}
	
	sequence = new string();
	
	if (expected_sequence_length > 0) {
		sequence->reserve(expected_sequence_length);
	}
	
	//cerr << "getSequence" << endl;
	
	//cout <<"line: "<< line << endl;
	while (true) {
		std::getline(*rna_stream, line);
		if (rna_stream->fail()) {
			break;
		}
		//cout << "fail: " << rna_stream->fail() << endl;
		
		//cerr << "y:" << line << endl;
		if (line.length() == 0) {
			continue;
		}
		
		if (line[0] == '#' ) {
			continue;
		}
		
		if (line[0] == '>' ) {
			
			//this->description_line = getFASTADescription(line);
			if (getOnlyIdentfier) {
				this->description_line = getFASTADescription(line);
			} else {
				this->description_line = line;
			}
			have_read_descr = false;
			sequence->reserve(sequence->length());
			return sequence;
		} else {
			int lastpos = line.length()-1;
			if (line[lastpos] == '\n') {
				line.erase(lastpos); // remove trailing newline
			}
			//sequence->append("@");
			//cout << "append" << endl;
			sequence->append(line);
		}
	}
	
	//file_closed = true;
	if (fp != 0) {
		pclose(fp);
		fp = 0;
	}
	sequence->reserve(sequence->length());
	return sequence;
	
}
