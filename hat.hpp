

#ifndef clustgun_hat_h
#define clustgun_hat_h

#include <omp.h>

#include "basic_tools_omp.hpp"



class ReaderWriterLock;


template <class T>
class HashedArrayTree { //  powers of two 2^20=1048576
protected:
	vector<T * > *hash; // the inner arrays will be fixed size, the outer vector flexibel
	size_t arrayChunkSize; // choose something like 1MB
	size_t two_power;
	size_t mod_mask;
	bool initialize_elements;
	T initialization_value;
	bool size_is_capacity;
	
	size_t currentArray;
	size_t nextFreePos;
	
public:
	#ifdef DEBUG
	string name;
	#endif

	HashedArrayTree<T> (bool size_is_capacity, size_t two_power) {
		this->initialize_elements = false;
		
		init(size_is_capacity, two_power);
	}
	
	HashedArrayTree<T> (bool size_is_capacity, size_t two_power, T initialization_value) {
		
		this->initialize_elements = true;
		this->initialization_value = initialization_value;
		
		init(size_is_capacity, two_power);
		
	}
	
	void init(bool size_is_capacity, size_t two_power) {
		
		this->size_is_capacity = size_is_capacity;
		this->two_power = two_power;
		this->arrayChunkSize = (size_t) pow((double)2, (double)two_power);
		this->mod_mask=this->arrayChunkSize -1 ; // e.g. converts 1000 into 111
#ifdef DEBUG
		try {
#endif
			hash = new vector<T * >();
#ifdef DEBUG
		} catch (bad_alloc& ba) {
			cerr << "error: (hash HAT) bad_alloc caught: " << ba.what() << endl;
			exit(1);
		}
#endif
		currentArray = 0;
		nextFreePos = 0;
		
		
	}
	
	~HashedArrayTree<T>() {
		//cout << "XX destructor" << this->name << endl;
		//cout << hash->size() << endl;
		for (int i = 0; i < (int) hash->size(); ++i){
			//cout << "currentArray:" << currentArray << endl;
			//cout << "delete: " << i << " " << (*hash)[i] << endl;
			//cout << "sizeof: " << i << " " << sizeof((*hash)[i]) << endl;
			delete [] (*hash)[i];
			//(*hash)[i]=NULL;
		}
		delete hash;
	}
	
	void deleteContentPointers() {
		//cout << "XX deleteContentPointers" << endl;
		//vector<T>::iterator vec_it;
		//typename vector<T * >::iterator vec_it;
		
		//cout << this->name << endl;
		if (!initialize_elements) {
			cerr << "error: deleting elements" << endl;
			return;
		}
		
		size_t effective_size;
		if (size_is_capacity) {
			effective_size = this->capacity();
		} else {
			effective_size = this->size();
		}
		
		
		for (size_t i = 0 ; i < effective_size; ++i ){
			
			#ifdef DEBUG
			try {
			if ( this->at(i) != this->initialization_value) {
				delete this->at(i);
			}
			} catch (out_of_range& oor) {
				cerr << "Out of Range error:(HAT: this->at(i) != this->initialization_value) " << oor.what() << endl;
				exit(1);
			}		
			#else
			if ( (*this)[i] != this->initialization_value) {
				delete (*this)[i];
			}
			#endif
			
		}

		
	}

	
	void deleteContentPairOfPointers() {
		//cout << "XX deleteContentPairOfPointers" << endl;
		//vector<T>::iterator vec_it;
		//typename vector<T * >::iterator vec_it;
		//int count = 0;
		//cout << "-------------2 " << *(this->at(0).second) << endl;
		
		for (int i = 0; i < currentArray ; ++i){
			for (int j = 0 ; j < arrayChunkSize; ++j) {
				//count++;
				//cout << "want deleteA: " << i << "," << j << endl;
				delete (*hash)[i][j].first;
				delete (*hash)[i][j].second;
			}
		}
		//i = currentArray-1 (last array)
		for (int j = 0 ; j < nextFreePos; ++j) {
			//count++;
			//cout << "want deleteB: " << currentArray << "," << j << endl;
			//cout << hash->size() << endl;
			//cout << *((*hash)[currentArray][j].second) << endl;
			
			delete (*hash)[currentArray][j].first;
			delete (*hash)[currentArray][j].second;
		}
		//cout << "Number of pairs deleted: " << count << endl;
		
	}
	
	
	size_t size() {
		if (size_is_capacity) {
			return capacity();
		}
		return ((currentArray)*arrayChunkSize + nextFreePos);
	}
	
	size_t capacity() {
		
		if (hash->size()*arrayChunkSize != hash->size() << two_power) { // TODO
			cerr << "hash->size()*arrayChunkSize != hash->size() << two_power" << endl;
			exit(1);
		}
		
		return hash->size() << two_power; // should be same as multiplying with "2^two_power"
		//return hash->size()*arrayChunkSize;
	}
	
	void reserve (size_t new_min_size) { // size, not position !!
		//cerr << "t: " << omp_get_thread_num() << " reserve " << new_min_size << endl; 
		size_t array = (new_min_size-1) >> two_power;
		
		while (capacity() < new_min_size) {
			T * small_vec;
			#ifdef DEBUG
			try {
			#endif
				small_vec = new T [arrayChunkSize];
			#ifdef DEBUG
			} catch (bad_alloc& ba) {
				cerr << "error: (HAT reserve) bad_alloc caught: " << ba.what() << endl;
				exit(1);
			}
			#endif
			if (initialize_elements) {
				for (size_t i =0; i < arrayChunkSize; ++i)  {
					small_vec[i]=this->initialization_value;
				}
			}
	//cout << "reserve push A " << small_vec << endl;
			#ifdef DEBUG
			try {
			#endif
				hash->push_back(small_vec);
			#ifdef DEBUG	
			} catch (bad_alloc& ba) {
				cerr << "error: (HAT reserve, push_back) bad_alloc caught: " << ba.what() << endl;
				exit(1);
			}
			#endif
		}
		
		if (array >= hash->size() ) { // hash->size() needs to be bigger than array // TODO in debug
			cerr << "error: (reserve HAT) array+ >= hash->size() " << endl;
			cerr << "array: " << array << endl;
			cerr << "hash->size(): " << hash->size() << endl;
			cerr << "new_min_size: " << new_min_size << endl;
			cerr << "pos, new_min_size-1: " << new_min_size-1 << endl;
			cerr << "capacity(): " << capacity() << endl;
			
			exit(1);
		}
		
	}
	
	
	T& at(size_t x){
		// boundary check 

		size_t array = x >> two_power;         // shifts the insignificant bits away, e.g. two_power is 20 for one 1MB
		size_t arraypos = x & this->mod_mask;  // I get the least significant bits for local array
		
		if (not size_is_capacity) {
		// check size !
		//	reserve(x);
			if (array > currentArray) {
				cerr << "array > currentArray" << endl;
				exit(1);
			}
			if (array == currentArray && arraypos >= nextFreePos) {
				cerr << "arraypos >= nextFreePos" << endl;
				cerr << "array: " << array << endl;
				cerr << "arraypos: " << arraypos << endl;
				cerr << "x: " << x << endl;
				cerr << "size: " << this->size() << endl;
				exit(1);
			}
			
			
		} else {
		// check capacity!
			if (array >= hash->size()) { // DO NOT CONFUSE WITH this->size() !!!!
				cerr << "array > hash->size()" << endl;
				cerr << "array: " << array << endl;
				cerr << "hash->size(): " << hash->size() << endl;
				cerr << "arraypos: " << arraypos << endl;
				cerr << "x: " << x << endl;
				cerr << "t: " << omp_get_thread_num() << endl;
				exit(1);
			}
			
			if (arraypos >= arrayChunkSize) {
				cerr << "exit: arraypos >= arrayChunkSize" << endl;
				exit(1);
			}
		}
		
		//}
		
		
		T* yyy;
		try {
			yyy =hash->at(array);
		} catch (out_of_range& oor) {
			cerr << "Out of Range error:(HAT: yyy =hash->at(array); " << oor.what() << endl;
			cerr << "arraypos >= nextFreePos" << endl;
			cerr << "array: " << array << endl;
			cerr << "arraypos: " << arraypos << endl;
			cerr << "x: " << x << endl;
			cerr << "size: " << this->size() << endl;
			cerr << "hash->size(): " << hash->size() << endl;
			exit(1);
		}
		
		T& xxx = yyy[arraypos ];
		return xxx;
		
	
		
		//return (*hash)[x >> two_power][ x & this->mod_mask];
		
	}
	
	
	inline T& operator[] (size_t x) {
		
		#ifdef DEBUG
		return this->at(x);
		#endif
		
		return (*hash)[x >> two_power][ x & this->mod_mask]; 
	}
	
#ifdef DEBUG 
	int push_back(T element) {
#else
	void push_back(T element) {
#endif
		//cerr << "lastPos: " << nextFreePos << endl;
		
		T * small_vec;
		
		if (currentArray == 0 && nextFreePos == 0) {
			#ifdef DEBUG
			try {
			#endif
				small_vec = new T [arrayChunkSize];
				//cout << "A create array of size " << 	arrayChunkSize	 << endl;
			#ifdef DEBUG	
			} catch (bad_alloc& ba) {
				cerr << "error: (HAT) bad_alloc caught: " << ba.what() << endl;
				exit(1);
			}
			#endif
			if (this->initialize_elements) {
				for (size_t i =0; i < arrayChunkSize; ++i)  {
					small_vec[i]=this->initialization_value;
				} 
			}
			//cout << "new vector" << endl;
			#ifdef DEBUG
			try {
			#endif
				hash->push_back(small_vec);
			#ifdef DEBUG	
			} catch (bad_alloc& ba) {
				cerr << "error: (HAT push_back,push_back) bad_alloc caught: " << ba.what() << endl;
				exit(1);
			}	
			#endif
			//cout << hash->size() << endl;
			//cout << "push was successful" << endl;
		// choose to create new or use old vector:
		} else if (nextFreePos >= arrayChunkSize) {
			// new vector
			//cout << "1lastPos: " << nextFreePos << endl;
			 
			if (hash->size()-1 <= (currentArray+1) ) { // it might have been reserved earlier
				#ifdef DEBUG
				try {
				#endif
					small_vec = new T [arrayChunkSize];
	//cout << "B create array of size " << 	arrayChunkSize	 << endl;				
				#ifdef DEBUG	
				} catch (bad_alloc& ba) {
					cerr << "error: (HAT) bad_alloc caught: " << ba.what() << endl;
					exit(1);
				}
				#endif
				
				if (this->initialize_elements) {
					for (size_t i =0; i < arrayChunkSize; ++i)  {
						small_vec[i]=this->initialization_value;
					} 
				} 
				
				//cout << "new vector!" << endl;
				#ifdef DEBUG
				try {
				#endif
					hash->push_back(small_vec);
				#ifdef DEBUG
				} catch (bad_alloc& ba) {
					cerr << "error: (HAT push_back, push_back) bad_alloc caught: " << ba.what() << endl;
					exit(1);
				}
				#endif
			}
			
			//cout << "hash->size(): " << hash->size() << endl;
			currentArray++;
			nextFreePos=0;
			
		} else {
			// old vector
			//cout << "old vector: " << currentArray << endl;
			small_vec = (*hash)[currentArray];
			//cout << "done" << endl;
		}
		
		
		//cout << "write: " << nextFreePos << endl;
		small_vec[nextFreePos]=element;
		//*this)[nextFreePos] = element;
		//cout << "small_vec->size(): " << small_vec->size() << endl;
		nextFreePos++;
		#ifdef DEBUG
		return nextFreePos-1;
		#endif
	}	
	

	
};

	
class HashedArrayTreeString : public HashedArrayTree<char> {

	
	public: 
		HashedArrayTreeString(size_t two_power) :  HashedArrayTree<char> (false, two_power) {
			reserve(1);
		};	
	
		char * addSequence(const char * seq, int seqlen) {
			
			if (seqlen+1 >= (int) arrayChunkSize) {
				cerr << "error: seqlen+1 >= arrayChunkSize" << endl;
				cerr << " input sequence is bigger than internal arrays, please increase default array size" << endl;
				exit(1);
			}
			#ifdef DEBUG
			if ((int)strlen(seq) != seqlen) {
				cerr << "error: strlen(seq) != seqlen:  "<< strlen(seq) << " " << seqlen << endl;
				exit(1);
			}
			#endif
			if (nextFreePos+seqlen+1 >= arrayChunkSize) {
				// go in next chunk, leave rest of current array empty...
				//cout << endl << seq <<endl;
				//cout << "nextFreePos: " << nextFreePos << endl;
				//cout << "arrayChunkSize: " << arrayChunkSize << endl;
				//cout << "seqlen: " << seqlen << endl;
				//cout << "this->size(): " << this->size() << endl;
				
				reserve(hash->size()*arrayChunkSize + 1);
				currentArray++;
				nextFreePos = 0;
			}
			
			char * small_vec = (*hash)[currentArray];
			
			char * buffer_start = &small_vec[nextFreePos];
			
			//memcpy(buffer_start, seq, seqlen)		
			strcpy(buffer_start, seq);
			
			nextFreePos += seqlen+1;
			
			return buffer_start;
			
			
		}
	
		char * addData(const char * seq, int datalen) {
			
			if (datalen >= (int) arrayChunkSize) {
				cerr << "error: seqlen+1 >= arrayChunkSize" << endl;
				cerr << " input sequence is bigger than internal arrays, please increase default array size" << endl;
				exit(1);
			}
	
			if (nextFreePos+datalen >= arrayChunkSize) {
				// go in next chunk, leave rest of current array empty...
				//cout << endl << seq <<endl;
				//cout << "nextFreePos: " << nextFreePos << endl;
				//cout << "arrayChunkSize: " << arrayChunkSize << endl;
				//cout << "seqlen: " << seqlen << endl;
				//cout << "this->size(): " << this->size() << endl;
				
				reserve(hash->size()*arrayChunkSize + 1);
				currentArray++;
				nextFreePos = 0;
			}
			
			char * small_vec = (*hash)[currentArray];
			
			char * buffer_start = &small_vec[nextFreePos];
			
			memcpy(buffer_start, seq, datalen); // datalen would include terminal \0	
			// old: strcpy(buffer_start, seq);
			
			nextFreePos += datalen;
			
			return buffer_start;
			
			
		}
	
};

template <class T>
class HashedArrayTree_rwlock : public HashedArrayTree<T>, public ReaderWriterLock {
	public:
	
	HashedArrayTree_rwlock<T>(bool size_is_capacity, int thread_count): HashedArrayTree<T>(), ReaderWriterLock(thread_count) {};
	HashedArrayTree_rwlock<T>(bool size_is_capacity, int thread_count, size_t two_power): HashedArrayTree<T>(size_is_capacity, two_power), ReaderWriterLock(thread_count) {};
	HashedArrayTree_rwlock<T>(bool size_is_capacity, int thread_count, size_t two_power, T initialization_value): HashedArrayTree<T>(size_is_capacity, two_power, initialization_value), ReaderWriterLock(thread_count) {};
	
};

#endif
