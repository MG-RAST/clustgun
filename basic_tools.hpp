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

using namespace std;


void process_mem_usage(double& vm_usage, double& resident_set);

std::string exec(const char* cmd);
double round( double value );
string stream2string(std::stringstream & stream);
void checkFile(string * file);
bool FileExists(const char * filename);
int str2int(string& text);
string int2str (int n);
void reverseComplementDNA(string & sequence);


//template <class T1, class T2, class T3> class triplet;

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
};


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
};




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
};

template <class A, class B>
void DeleteMapWithPointers (multimap< A, B> * my_map) {
	
	
    if (my_map == 0) return;
	for_each( my_map->begin(), my_map->end(), DeleteMapFntor<A, B>() );
	delete my_map;
	
	return;
};


template <class A, class B>
pair< typename map< A , B >::iterator , bool > mapInsert (map< A, B> * mymap, A a , B b) {
	pair< typename map< A , B >::iterator , bool > ret;

	ret = mymap->insert(pair<A,B>(a, b));
		
	return ret;
};


//template <class A, class B>
//typename multimap<A, B>::iterator  mmapInsert (multimap< A, B> * mymap, A a , B b) {
//	multimap< A , B >::iterator ret;
	
//	ret = mymap->insert(pair<A,B>(a, b));
	
//	return ret;
//};

#endif
