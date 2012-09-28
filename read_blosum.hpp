#ifndef read_blosum_h
#define read_blosum_h


#include <cstdlib>
#include <fstream>
#include <string>
#include <iostream>

#include "basic_tools.hpp"

using namespace std;

short int * readBLOSUM(const char * blosum_file);
short int getScore(char a, char b, short int * blosum_matrix);


#endif
