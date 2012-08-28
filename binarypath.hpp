#ifndef binary_path_h
#define binary_path_h


#include <cstdlib>
#include <cstring>
#include <iostream>

//OSX
#ifdef __APPLE__
#include <mach-o/dyld.h>

#else
//Linux
#include <stdio.h>
#include <unistd.h>
#endif


using namespace std;



string getbinarypath();


#endif
