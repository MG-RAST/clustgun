

#ifndef clustgun_hat_h
#define clustgun_hat_h



template <class T>
class HashedArrayTree { //  powers of two 2^20=1048576
private:
	vector<T * > *hash; // the inner vector will be fixed size, the outer flexibel
	size_t arrayChunkSize; // choose something like 1MB
	size_t two_power;
	size_t mod_mask;
	bool initialize_elements;
	T initialization_value;
	
	
	size_t currentArray;
	size_t nextFreePos;
	
public:
	#ifdef DEBUG
	string name;
	#endif
	
	~HashedArrayTree<T>() {
		//cout << "XX destructor" << this->name << endl;
		//cout << hash->size() << endl;
		for (int i = 0; i<hash->size(); ++i){
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
		
		
		for (int i = 0; i<hash->size()-1; ++i){
			//cout << "i: " << i << endl;
			for (int j = 0 ; j < arrayChunkSize; ++j) {
				//if (j > 1040000) {
				//	cout << i << " - " << j << endl;
				//	cout << (*hash)[i][j] << endl;
				//}
				if ((*hash)[i][j] != NULL ) {
					//cout << i << " - " << j << endl;
					//cout << (*hash)[i][j] << endl;
					delete (*hash)[i][j];
				}
			}
			
		}
		
		for (int j = 0 ; j < nextFreePos; ++j) {
			//cout << "j: " << j << endl;
			//if (j % 10000 == 0) {
			//	cout << (hash->size()-1) << " -- " << j << endl;
			//	cout << (*hash)[hash->size()-1][j] << endl;
			//}
			if ((*hash)[hash->size()-1][j] != NULL ) {
				
				
				delete (*hash)[hash->size()-1][j];
			}
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
		return ((currentArray)*arrayChunkSize + nextFreePos);
	}
	
	
	void reserve (size_t x) {
		
		size_t requested_array = (x >> two_power);
//cout << "reserve" << endl;
		while (hash->size() < requested_array) {
			T * small_vec;
			try {
				small_vec = new T [arrayChunkSize];
			} catch (bad_alloc& ba) {
				cerr << "error: (HAT reserve) bad_alloc caught: " << ba.what() << endl;
				exit(1);
			}
			if (initialize_elements) {
				for (size_t i =0; i < arrayChunkSize; ++i)  {
					small_vec[i]=this->initialization_value;
				}
			}
			//cout << "reserve push A " << small_vec << endl;
			try {
				hash->push_back(small_vec);
			} catch (bad_alloc& ba) {
				cerr << "error: (HAT reserve, push_back) bad_alloc caught: " << ba.what() << endl;
				exit(1);
			}	
		}
		
		if (requested_array*two_power < x) {
			T * small_vec;
			try {
				small_vec = new T [arrayChunkSize];
			} catch (bad_alloc& ba) {
				cerr << "error: (HAT reserve) bad_alloc caught: " << ba.what() << endl;
				exit(1);
			}	
			if (initialize_elements) {
				for (size_t i =0; i < arrayChunkSize; ++i)  {
					small_vec[i]=this->initialization_value;
				}
			}
			//cout << "reserve push B" << small_vec  << endl;
			try {
				hash->push_back(small_vec);
			} catch (bad_alloc& ba) {
				cerr << "error: (HAT reserve, push_back) bad_alloc caught: " << ba.what() << endl;
				exit(1);
			}	
		}
		//cout << "reserve len: " << hash->size() << endl;
		
	}
	
	
	inline T& at(size_t x){
		// boundary check and increase memory if needed
		if (initialize_elements) {
			reserve(x);
		} else {
			if ((x >> two_power) > currentArray) {
				cerr << "(x >> two_power) > currentArray" << endl;
				
				exit(1);
			}
			if ((x >> two_power) == (currentArray) && (x & this->mod_mask) >= nextFreePos) {
				cerr << "(x & this->mod_mask) >= nextFreePos" << endl;
				exit(1);
			}
		}
		
		return hash->at(x >> two_power)[ x & this->mod_mask];
	}
	
	
	inline T& operator[] (size_t x) {
		
		
		
		//size_t array = x / arrayChunkSize;
		//size_t arraypos = arrayChunkSize % x;
		
		//size_t array = x >> two_power;
		//size_t arraypos = x & this->mod_mask;
		
		
		//cout << "x: " << x<< endl;
		//cout << "this->mod_mask: " << this->mod_mask<< endl;
		//cout << "array: " << array<< endl;
		//cout << "  arraypos: " <<arraypos << endl;
		
		
		
		return (*hash)[x >> two_power][ x & this->mod_mask]; 
	}
	
	void push_back(T element) {
		//cout << "lastPos: " << nextFreePos << endl;
		
		T * small_vec;
		
		if (currentArray == 0 && nextFreePos == 0) {
			try {
				small_vec = new T [arrayChunkSize];
			} catch (bad_alloc& ba) {
				cerr << "error: (HAT) bad_alloc caught: " << ba.what() << endl;
				exit(1);
			}
			if (this->initialize_elements) {
				for (size_t i =0; i < arrayChunkSize; ++i)  {
					small_vec[i]=this->initialization_value;
				} 
			}
			//cout << "new vector" << endl;
			try {
				hash->push_back(small_vec);
			} catch (bad_alloc& ba) {
				cerr << "error: (HAT push_back,push_back) bad_alloc caught: " << ba.what() << endl;
				exit(1);
			}	
			//cout << hash->size() << endl;
			//cout << "push was successful" << endl;
		// choose to create new or use old vector:
		} else if (nextFreePos >= arrayChunkSize) {
			// new vector
			//cout << "1lastPos: " << nextFreePos << endl;
			 
			if (hash->size()-1 <= (currentArray+1) ) { // it might have been reserved earlier
				try {
					small_vec = new T [arrayChunkSize];
				} catch (bad_alloc& ba) {
					cerr << "error: (HAT) bad_alloc caught: " << ba.what() << endl;
					exit(1);
				}
				
				if (this->initialize_elements) {
					for (size_t i =0; i < arrayChunkSize; ++i)  {
						small_vec[i]=this->initialization_value;
					} 
				} 
				
				//cout << "new vector!" << endl;
				try {
					hash->push_back(small_vec);
				} catch (bad_alloc& ba) {
					cerr << "error: (HAT push_back, push_back) bad_alloc caught: " << ba.what() << endl;
					exit(1);
				}
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
		
	}	
	
	HashedArrayTree<T> (size_t two_power) {
		this->initialize_elements = false;
		
		this->two_power = two_power; 
		this->arrayChunkSize = (size_t) pow((double)2, (double)two_power); 
		this->mod_mask=this->arrayChunkSize -1 ; // e.g. converts 1000 into 111
		try {
			hash = new vector< T * >();
		} catch (bad_alloc& ba) {
			cerr << "error: (HAT hash) bad_alloc caught: " << ba.what() << endl;
			exit(1);
		}
		
		currentArray = 0;
		nextFreePos = 0;
		//this->reserve(arrayChunkSize);
	};
	
	HashedArrayTree<T> (size_t two_power, T initialization_value) {
		
		this->initialization_value = initialization_value;
		this->initialize_elements = true;
		
		this->two_power = two_power; 
		this->arrayChunkSize = (size_t) pow((double)2, (double)two_power); 
		this->mod_mask=this->arrayChunkSize -1 ; // e.g. converts 1000 into 111
		try {
			hash = new vector<T * >();
		} catch (bad_alloc& ba) {
			cerr << "error: (hash HAT) bad_alloc caught: " << ba.what() << endl;
			exit(1);
		}
		
		currentArray = 0;
		nextFreePos = 0;
		//this->reserve(arrayChunkSize);
		
	};
	
	
};



#endif
