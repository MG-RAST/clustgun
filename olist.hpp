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

using namespace std;

// unrolled linked list

const int lengthofmyolistarrays = 10;


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
	OListArray<T1, T2> *start;
	OListArray<T1, T2> *end;
	
	int lastArrayPosition; // last element in last array
	int length;
	
	// for iteration purposes; using real c++ iterator would be nice:
	OListArray<T1, T2> *currentListArray;
	int currentArrayPosition;
	
	OList<T1, T2>(){
		start=NULL;
		end=NULL;
		lastArrayPosition=-1;
		length = 0;
		
		currentListArray = NULL;
		currentArrayPosition = -1;
	};
	
	~OList<T1, T2>(){
		
		currentListArray = start;
		while (currentListArray != NULL) {
			OListArray<T1, T2> *delArray = currentListArray;
			currentListArray = currentListArray->nextArray;
			delete delArray;
		}
		
	}
	
	T1& getFirst(){
		return this->currentListArray->data_array1[this->currentArrayPosition];
	}
	
	T2& getSecond(){
		return this->currentListArray->data_array2[this->currentArrayPosition];
	}
	
	
	void eraseElement() {
		
		// overwrite entry with last entrythis->
		if (this->currentListArray != this->end || this->currentArrayPosition != lastArrayPosition) { // not necessary if we want do delete the last element
			//cout << "1) this->currentListArray->data_array1[this->currentArrayPosition]: " << this->currentListArray->data_array1[this->currentArrayPosition]<< endl;
			this->currentListArray->data_array1[this->currentArrayPosition]=this->end->data_array1[this->lastArrayPosition];
			//cout << "2) this->currentListArray->data_array1[this->currentArrayPosition]: " << this->currentListArray->data_array1[this->currentArrayPosition]<< endl;
			this->currentListArray->data_array2[this->currentArrayPosition]=this->end->data_array2[this->lastArrayPosition];
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
			OListArray<T1, T2> * secondLastArray = this->end->prevArray;
			
			delete this->end;
			this->end = secondLastArray; // might also be NULL, ok
			//cout << "this->end: " << this->end<< endl;
			if (this->end != NULL) {
				this->end->nextArray = NULL;
				this->lastArrayPosition = lengthofmyolistarrays-1;
			} else {
				this->start = NULL; // list contains no arrays
				this->lastArrayPosition = -1;
			}
			
		}
		this->length--;
		
	}
	
	void resetIterator(){
		currentListArray = start;
		currentArrayPosition = -1;
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
		while (currentListArray_tmp != end) {
			for (int i =0 ; i< lengthofmyolistarrays; ++i) {
				cout << (length + i) << ": " <<  currentListArray_tmp->data_array1[i] << endl;
				
			}
			length += lengthofmyolistarrays;
			currentListArray_tmp = currentListArray_tmp->nextArray;
		}
		for (int i =0 ; i<= lastArrayPosition; ++i) {
			cout << (length + i) << ": (" <<  currentListArray_tmp->data_array1[i] << "," << currentListArray_tmp->data_array2[i] << ")"<< endl;
			
		}
		length += (lastArrayPosition+1);
		return;
	}
	
	// iterate through list and array elements
	bool nextElement() {
		if (currentListArray == NULL) return false;
		
		if (currentListArray->nextArray != NULL) {
			
			if (currentArrayPosition < lengthofmyolistarrays-1) {
				this->currentArrayPosition++;
			} else {
				//cout << "--" << endl;
				this->currentArrayPosition = 0;
				currentListArray=currentListArray->nextArray;
			}	
		} else {
			if (currentArrayPosition < lastArrayPosition) {
				this->currentArrayPosition++;
			} else {
				return false;
			}
		}
		
		
		return true;
	}
	
	bool previousElement() {
		if (this->currentListArray == NULL) return false;
		
		if (this->currentArrayPosition > 0) {
			this->currentArrayPosition--;
			
		} else {
			
			if (this->currentListArray->prevArray == NULL) {
				return false;
			}
			
			this->currentListArray = this->currentListArray->prevArray;
			this->currentArrayPosition = lengthofmyolistarrays - 1;
		}
		
		
		return true;
	} 
	
	void append(T1 a, T2 b){
		//cout << "append" << endl;
		
		if (this->start == NULL) {
			try {
				this->start = new OListArray<T1, T2>(); 
			} catch (bad_alloc& ba) {
				cerr << "error: (olist) bad_alloc caught: " << ba.what() << endl;
				exit(1);
			}
			this->end = this->start;
			this->lastArrayPosition = 0;
			//cout << "appendcreatefirst" << endl;
		} else if (lastArrayPosition < lengthofmyolistarrays-1) {
			this->lastArrayPosition++;
			//cout << "appendext" << endl;
		} else {
			if (this->end->nextArray != NULL) {
				cerr << "error: (append) this->end->nextArray != NULL" << endl;
				exit(1);	
			}
			try {
				this->end->nextArray = new OListArray<T1, T2>();
			} catch (bad_alloc& ba) {
				cerr << "error: (olist) bad_alloc caught: " << ba.what() << endl;
				exit(1);
			}
			this->end->nextArray->prevArray = this->end;
			this->end = this->end->nextArray;
			this->lastArrayPosition = 0;
			
			//cout << "appendnew-----" << endl;
		}
		
		this->end->data_array1[lastArrayPosition]=a;
		this->end->data_array2[lastArrayPosition]=b;
		
		this->length++;
		
		return;
		
	}
	
	
		
};


#endif
