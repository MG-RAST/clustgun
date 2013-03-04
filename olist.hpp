//
//  olist.hpp
//  clustgun
//
//  Created by Wolfgang Gerlach on 5/9/12.
//  Copyright (c) 2012 University of Chicago. All rights reserved.
//

#ifndef clustgun_olist_hpp
#define clustgun_olist_hpp

//#include <iterator>

#include <boost/iterator/iterator_facade.hpp> // header-only lib

#include <boost/dynamic_bitset.hpp>

#include "basic_tools_omp.hpp" // for the lock


using namespace std;



// unrolled linked list

//const int lengthofmyolistarrays = 10;
const int lengthofmyolistarrays = 32;

template <class T1, class T2>
class OList_iterator;


template <class T1, class T2>
class OListArray {
public:
	OListArray *nextArray;
	OListArray *prevArray;
	// array of pairs would have been nice, but does not seem to be memory efficient
	T1 data_array1[lengthofmyolistarrays];
	T2 data_array2[lengthofmyolistarrays];
	
};





template <class T1, class T2>
class OList {
public:
	
	typedef OList_iterator<T1, T2> iterator;
	
	
	
	OListArray<T1, T2> *start;
	OListArray<T1, T2> *lastArray;
	
	int lastArrayPosition; // last element in last array
	int length;

	#ifdef DEBUG
	bool test_write_access;
	#endif
	
	OList<T1, T2>(){
		start=NULL;
		lastArray=NULL;
		lastArrayPosition=-1;
		length = 0;
		
		//currentListArray = NULL;
		//currentArrayPosition = -1;
#ifdef DEBUG
		test_write_access=false;
#endif
		
	};
	
	~OList<T1, T2>(){
		
		OListArray<T1, T2> * currentListArray = start;
		while (currentListArray != NULL) {
			OListArray<T1, T2> *delArray = currentListArray;
			currentListArray = currentListArray->nextArray;
			delete delArray;
		}
	}
	
	
	
	T1& getFirst_ext(OListArray<T1, T2> * currentListArray, int currentArrayPosition){
		
		#ifdef DEBUG
		if (currentListArray == NULL ) {
			cerr << "error: (getFirst) this->currentListArray == NULL" << endl;
			exit(1);
		}
		
		if (currentArrayPosition < 0) {
			cerr << "error: (getFirst) this->currentArrayPosition < 0" << endl;
			exit(1);
		}
		#endif
		
		return currentListArray->data_array1[currentArrayPosition];
	}
	
	T2& getSecond_ext(OListArray<T1, T2> * currentListArray, int currentArrayPosition){
		
		#ifdef DEBUG
		if (currentListArray == NULL ) {
			cerr << "error: (getSecond) currentListArray == NULL" << endl;
			exit(1);
		}
		
		if (currentArrayPosition < 0) {
			cerr << "error: (getSecond) currentArrayPosition < 0" << endl;
			exit(1);
		}
		#endif
		
		return currentListArray->data_array2[currentArrayPosition];
	}
	
	
	
	void erase(iterator & it) {
		OListArray<T1, T2> * currentListArray = it.position.first;
		int currentArrayPosition = it.position.second;
		eraseElement_ext(currentListArray, currentArrayPosition);
	}
	
	
	
	void eraseElement_ext(OListArray<T1, T2> * currentListArray, int currentArrayPosition) {
		
		// overwrite entry with last entrythis->
		if (currentListArray != this->lastArray || currentArrayPosition != this->lastArrayPosition) { // not necessary if we want do delete the last element
			
			//cout << "this->currentArrayPosition: " << this->currentArrayPosition << endl;
			
			//for (int i  = 0; i <= 5 ; ++i) {
			//	cout << i << " __" << this->currentListArray->data_array1[i] << endl;
			//	cout << i << " __" << this->currentListArray->data_array2[i] << endl;
			//}
			
			//cout << "1) this->currentListArray->data_array1[this->currentArrayPosition]: " << this->currentListArray->data_array1[this->currentArrayPosition]<< endl;
			
			
			currentListArray->data_array1[currentArrayPosition]=this->lastArray->data_array1[this->lastArrayPosition];
			//cout << "2) this->currentListArray->data_array1[this->currentArrayPosition]: " << this->currentListArray->data_array1[this->currentArrayPosition]<< endl;
			currentListArray->data_array2[currentArrayPosition]=this->lastArray->data_array2[this->lastArrayPosition];
		}
		
		// delete last element, to avoid duplicate
		this->pop_back();
		
	}
	
	
	// delete last element
	void pop_back(){
		if (this->lastArrayPosition>0) {
			this->lastArrayPosition--;
		} else {
			// have to delete last array completly:
			OListArray<T1, T2> * secondLastArray = this->lastArray->prevArray;
			
			delete this->lastArray;
			this->lastArray = secondLastArray; // might also be NULL, ok
			//cout << "this->lastArray: " << this->end<< endl;
			if (this->lastArray != NULL) {
				this->lastArray->nextArray = NULL;
				this->lastArrayPosition = lengthofmyolistarrays-1;
			} else {
				this->start = NULL; // list contains no arrays
				this->lastArrayPosition = -1;
			}
			
		}
		this->length--;
		
	}
	
	
	int getLength() {
		return this->length;
	}
	
	
	void print() {
		//cout << "this->currentArrayPosition: " << this->currentArrayPosition << endl;
		if (start == NULL) {
			cout << "empty" << endl;
			return;
		}
		
		OListArray<T1, T2> *currentListArray_tmp = start;
		int length = 0;
		while (currentListArray_tmp != lastArray) {
			for (int i =0 ; i< lengthofmyolistarrays; ++i) {
				//cout << (length + i) << ": " <<  currentListArray_tmp->data_array1[i] << endl;
				cout << (length + i) << ": (" <<  currentListArray_tmp->data_array1[i] << "," << currentListArray_tmp->data_array2[i] << ")"<< endl;
				
			}
			cout << "--" << endl;
			length += lengthofmyolistarrays;
			currentListArray_tmp = currentListArray_tmp->nextArray;
		}
		for (int i =0 ; i<= lastArrayPosition; ++i) {
			cout << (length + i) << ": (" <<  currentListArray_tmp->data_array1[i] << "," << currentListArray_tmp->data_array2[i] << ")"<< endl;
			
		}
		length += (lastArrayPosition+1);
		return;
	}
	
	// for iterators to use.. (move code into iterator ???)
	pair<OListArray<T1, T2> *, int> nextElement_extIt(OListArray<T1, T2> * currentListArray, int currentArrayPosition) {
		//cerr << "nextElement_extIt" << endl;
		if (currentListArray == NULL) return make_pair((OListArray<T1, T2> *) NULL, -1); //<OListArray<T1, T2> *, int>
		
		if (currentListArray->nextArray != NULL) {
			
			if (currentArrayPosition < lengthofmyolistarrays-1) {
				// increase in same array
				currentArrayPosition++;
			} else {
				//go in next array
				currentArrayPosition = 0;
				currentListArray=currentListArray->nextArray;
			}
		} else {
			// last array
			if (currentArrayPosition < lastArrayPosition) {
				currentArrayPosition++;
			} else {
				return make_pair(currentListArray, lastArrayPosition+1); // this points to the element after the last element, for end()-method //<OListArray<T1, T2> *, int
			}
		}
		
		return make_pair(currentListArray, currentArrayPosition); //<OListArray<T1, T2> *, int>
	}
	
	pair<OListArray<T1, T2> *, int>  previousElement_extIt(OListArray<T1, T2> * currentListArray, int currentArrayPosition) { 
		if (currentListArray == NULL) return false;
		
		if (currentArrayPosition > 0) {
			currentArrayPosition--;
			
		} else {
			
			if (currentListArray->prevArray == NULL) {
				return false;
			}
			
			currentListArray = currentListArray->prevArray;
			currentArrayPosition = lengthofmyolistarrays - 1;
		}
		
		
		return make_pair(currentListArray, currentArrayPosition); //<OListArray<T1, T2> *, int>
	} 
	
	void append(T1 a, T2 b){
		//cout << "append" << endl;
	
#ifdef DEBUG
#pragma omp critical(test_write_access)
		{
			if (test_write_access) {
				cerr << "error: olist test_write_access" << endl;
				exit(1);
			}
			
			test_write_access=true;
		}
		
#endif
		
		if (this->start == NULL) {
			OListArray<T1, T2> * temparray;
			#ifdef DEBUG
			try {
			#endif
				temparray = new OListArray<T1, T2>(); 
			#ifdef DEBUG
			} catch (bad_alloc& ba) {
				cerr << "error: (olist) bad_alloc caught: " << ba.what() << endl;
				exit(1);
			}
			#endif
			this->start = temparray;
			this->lastArray = this->start;
			this->lastArrayPosition = 0;
			
			//cout << "appendcreatefirst" << endl;
		} else if (lastArrayPosition < lengthofmyolistarrays-1) {
			this->lastArrayPosition++;
			//cout << "appendext" << endl;
		} else { // last position in array
			#ifdef DEBUG
			if (this->lastArray->nextArray != NULL) {
				cerr << "error: (append) this->lastArray->nextArray != NULL" << endl;
				exit(1);	
			}
			#endif
			
			OListArray<T1, T2> * newarray;
			#ifdef DEBUG
			try {
			#endif
				newarray = new OListArray<T1, T2>();
			#ifdef DEBUG
			} catch (bad_alloc& ba) {
				cerr << "error: (olist) bad_alloc caught: " << ba.what() << endl;
				exit(1);
			}
			#endif
			
			// link old end and new end:
			this->lastArray->nextArray = newarray;
			newarray->prevArray = this->lastArray;
			
			this->lastArray = newarray;
			this->lastArrayPosition = 0;
			
			//cout << "appendnew-----" << endl;
		}
		
		this->lastArray->data_array1[lastArrayPosition]=a;
		this->lastArray->data_array2[lastArrayPosition]=b;
		
		this->length++;
		
#ifdef DEBUG
#pragma omp critical(test_write_access)
		{
			if (not test_write_access) {
				cerr << "error: olist NOT test_write_access" << endl;
				exit(1);
			}
			
			test_write_access=false;
		}
		
#endif
		
		return;
		
	}
	
	
	
	iterator begin() {
		if (start != NULL) {
			
			return iterator( this, make_pair(start, (int) 0) ); // <OListArray<T1, T2> *, int>
		}
		
		//return iterator( this, make_pair(lastArray, lastArrayPosition+1) );
		return end();
	}
	
	iterator end() {
		return iterator( this, make_pair(lastArray, lastArrayPosition +1) ); // has to point to element after last element...
	}
	
	
	
};



// OList with reader-writer lock
template <class T1, class T2>
class OList_rwlock : public ReaderWriterLock, public OList<T1,T2> {
	public:
	OList_rwlock(){};
	OList_rwlock(int thread_count): ReaderWriterLock(thread_count) {};
	
};


// this iterator has the problem that is cannot increment without a pointer to the main data object
template <class T1, class T2>
class OList_iterator : public boost::iterator_facade<	OList_iterator<T1, T2>,
														pair< OListArray<T1, T2> *, int>,
														boost::forward_traversal_tag
> {
	public:
		
		typedef pair< OListArray<T1, T2> *, int> OListPosPair;
		typedef OList_iterator<T1, T2> _base;

		OList_iterator(){}
	
//		OList_iterator()
//			:	olist(0)
//		{
//			position = make_pair((OListArray<T1, T2> *) NULL, 0) ; //<OListArray<T1, T2> *, int >
//		}
//	
//		OList_iterator(OList<T1, T2> * list)
//			:	olist(list)
//		{
//			position = make_pair((OListArray<T1, T2> *) NULL, 0) ; //< OListArray<T1, T2> *, int >
//		}
//	
		OList_iterator(OList<T1, T2> * list, OListPosPair position)
			:	olist(list),
				position(position)
		{}
			
	

	
		T1& getFirst(){
			return olist->getFirst_ext(position.first, position.second);
		}
			
		T2& getSecond(){
			return olist->getSecond_ext(position.first, position.second);
		}
	

	
	private:
	
		
	
		OList<T1, T2> * olist;
		OListPosPair position;
		
	
		friend class boost::iterator_core_access;
		friend class OList<T1, T2>; // Olist erase function needs access to position
	
		void increment() {
			//m_node = m_node->next();
			position = olist->nextElement_extIt(position.first, position.second);
		}
	
		bool equal(_base const& other) const{
			return  ( this->position == other.position );
			//return  ( (this->position.first == other.position.first) && (this->position.second == other.position.second) );
			//return this->m_node == other.m_node;
		}
	
		
		
		// returns a pair
		OListPosPair& dereference() const {
			
			return make_pair< T1, T2>(getFirst_ext(position.first, position.second), getSecond_ext(position.first, position.second));
			//return *m_node;
		}
};


#endif
