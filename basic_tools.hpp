#ifndef basic_tools_h
#define basic_tools_h


#include <string>
#include <sstream>
#include <iostream>
#include <sys/stat.h>
#include <algorithm>
#include <map>
#include <unistd.h>
#include<fstream>
#include <cmath>
#include <time.h>


using namespace std;





std::string getFileExtension(const std::string& FileName);
std::string getFileNameWithoutExtension(const std::string& FileName);

void process_mem_usage(double& vm_usage, double& resident_set);


timespec diff(timespec start, timespec end);
timespec add_time(timespec time1, timespec time2);

int msleep(unsigned long milisec);
int nsleep(unsigned long nanosec);
std::string exec(const char* cmd);
double round( double value );
string stream2string(std::stringstream & stream);
void checkFile(string * file);
bool FileExists(const char * filename);
int str2int(string& text);
string int2str (int n);
void reverseComplementDNA(string & sequence);
int * get_nucleobase_ascii2num();

//template <class T1, class T2, class T3> class triplet;

inline int ipow(int base, int exp)
{
    int result = 1;
    while (exp)
    {
        if (exp & 1)
            result *= base;
        exp >>= 1;
        base *= base;
    }
	
    return result;
}

template <class T1, class T2, class T3>
class triplet {
public:
	T1 first;
	T2 second;
	T3 third;
	
	triplet(T1 a, T2 b, T3 c) {
		first = a;
		second = b;
		third = c;
	}
	
	triplet() {}
	
};



template <class T>
void DeleteObject (T * myobject) {
	
	
    if (myobject == 0) return;
	
	delete myobject;
	
	return;
}


struct delete_object
{
public:
	template <typename T>
	void operator()(T *ptr){ DeleteObject(ptr);}
};


template <class T>
void DeleteContainerWithPointers (T * container) {
	
	
    if (container == 0) return;
	for_each( container->begin(), container->end(), delete_object() );
	delete container;
	
	return;
}




// Functor for deleting pointers in map.
template<class A, class B>
struct DeleteMapFntor
{
    // Overloaded () operator.
    // This will be called by for_each() function.
    bool operator()(pair<A,B> x) const
    {
        // Assuming the second item of map is to be
        // deleted. Change as you wish.
        DeleteObject( x.second );
        return true;
    }
};

template <class A, class B>
void DeleteMapWithPointers (map< A, B> * my_map) {
	
	
    if (my_map == 0) return;
	for_each( my_map->begin(), my_map->end(), DeleteMapFntor<A, B>() );
	delete my_map;
	
	return;
}

template <class A, class B>
void DeleteMapWithPointers (multimap< A, B> * my_map) {
	
	
    if (my_map == 0) return;
	for_each( my_map->begin(), my_map->end(), DeleteMapFntor<A, B>() );
	delete my_map;
	
	return;
}


template <class A, class B>
pair< typename map< A , B >::iterator , bool > mapInsert (map< A, B> * mymap, A a , B b) {
	pair< typename map< A , B >::iterator , bool > ret;

	ret = mymap->insert(pair<A,B>(a, b));
		
	return ret;
}


//template <class A, class B>
//typename multimap<A, B>::iterator  mmapInsert (multimap< A, B> * mymap, A a , B b) {
//	multimap< A , B >::iterator ret;
	
//	ret = mymap->insert(pair<A,B>(a, b));
	
//	return ret;
//};

#endif
