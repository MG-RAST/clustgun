
#ifndef fasta_parser_h
#define fasta_parser_h


#include <iostream>
#include <fstream>

#include <boost/iostreams/device/file_descriptor.hpp>
#include <boost/iostreams/stream.hpp>


using namespace std;


string getFASTADescription (string & line);

class FASTA_Parser {

private:
	string * sequence;
	
public:
	
	FILE *fp;
	std::istream * rna_stream;
	boost::iostreams::stream_buffer<boost::iostreams::file_descriptor_source> * fpstream;
	
	string description_line;
	bool have_read_descr;
	int expected_sequence_length;
	bool getOnlyIdentfier;
	//bool file_closed;
	
	FASTA_Parser(string command, bool getOnlyIdentfier);
	FASTA_Parser(string file, bool getOnlyIdentfier, string * zcat);
	FASTA_Parser(string file, bool getOnlyIdentfier, bool use_zcat);
	
	void start(string command, bool getOnlyIdentfier);
	
	~FASTA_Parser();
	string * getSequence();
	bool getNextDescriptionLine(string& descr);
	
	
};

#endif
