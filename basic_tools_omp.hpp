#ifndef basic_tools_omp_h
#define basic_tools_omp_h

#include <omp.h>
#include <bitset>

#include "basic_tools.hpp"

#include "main.h"


const int max_thread_count=100; // warning: change value also in main.hpp
//extern int max_thread_count;

#ifdef DEBUG_DEADLOCK
extern int thread_status[max_thread_count];
#endif
extern int reads_per_thread[max_thread_count];

#ifdef DEBUG_DEADLOCK
extern int thread_count;
#endif

#ifdef DEBUG_DEADLOCK
void omp_set_lock_with_timeout(omp_lock_t * lock);
#endif

class ReaderWriterLock {
	
	
private:
	
	//boost::dynamic_bitset<> flush_needed;
	#ifdef USE_FLUSH_NEEDED
	bitset<CHAR_BIT> flush_needed;
	#endif
	
	#ifdef DEBUG_DEADLOCK
	//boost::dynamic_bitset<> reader_indicator;
	//bitset<CHAR_BIT> reader_indicator; // TODO need more flexible and efficient solution than CHAR_BIT!?
	vector<bool> reader_indicator;
	#endif
	
	
	omp_lock_t writer_lock;
	short reader_counter;
	
	#ifdef DEBUG
	int writer_thread;
	#endif
public:
	
	ReaderWriterLock(){}; // dummy constructor
	
	ReaderWriterLock(int threadcount){
		omp_init_lock(&writer_lock);
		
		this->reader_counter = 0; 
		#ifdef USE_FLUSH_NEEDED
		//flush_needed = boost::dynamic_bitset<>(threadcount);
		flush_needed = bitset<CHAR_BIT>();
		#endif
		
		
		#ifdef DEBUG_DEADLOCK
		#pragma omp critical(reader_indicator)
		{
		//reader_indicator  = boost::dynamic_bitset<>(threadcount, 0);
			//reader_indicator=bitset<CHAR_BIT>();
			reader_indicator=vector<bool>(CHAR_BIT, false);  // should be faster than bitset
		}
		#endif

	}
	
	
	void flush_if_needed() {
		#ifdef USE_FLUSH_NEEDED
		int threadid = omp_get_thread_num();
		#pragma omp critical(flushbit)
		{
			if (flush_needed[threadid]) {
			
				flush_needed[threadid]=0; // this should not be called concurrently.
			
			}
		}
		#else
		#pragma omp flush
		#endif
	}
	
	
	
	
	void set_reader_lock() {
		//#pragma omp atomic
		//reader_lock++;
		
		#ifdef DEBUG_DEADLOCK
		timespec time_start;
		timespec time_end;
		clock_gettime(CLOCK_REALTIME, &time_start);  // was CLOCK_PROCESS_CPUTIME_ID
		#endif
		
		// if reader can get write lock, continue...
		while ( ! omp_test_lock(&writer_lock)) {
			
			#ifdef DEBUG_DEADLOCK
			clock_gettime(CLOCK_REALTIME, &time_end);
			
			
			
			if ( diff(time_start, time_end).tv_sec > 50 ) {
				#pragma omp critical(cerr)
				{
					cerr << "thread: " << (int) omp_get_thread_num() << " warning: I wait for 5 sec now..! reader_counter:"<< reader_counter << endl;
					
											
					for (int tt = 0; tt < thread_count ; tt++) {
						cerr << thread_status[tt] << ",";
						
					}
					cerr << endl;
					for (int tt = 0; tt < thread_count ; tt++) {
						cerr << reads_per_thread[tt] << ",";
						
					}
					cerr << endl;
				
					for (int tt = 0; tt < thread_count ; tt++) {
						cerr << reader_indicator[tt] << ",";
						
					}
					cerr << endl;
				
					
				}
			}
		
			
			if ( diff(time_start, time_end).tv_sec > 5 ) {
				#pragma omp critical(cerr)
				{
					cerr << "thread: " << (int) omp_get_thread_num() << " timeout in wait for read access! " << endl;
						
												
					for (int tt = 0; tt < 10 ; tt++) {
						cerr << thread_status[tt] << ",";
						
					}
					cerr << endl;
					for (int tt = 0; tt < 10 ; tt++) {
						cerr << reads_per_thread[tt] << ",";
						
					}
					cerr << endl;
				
					for (int tt = 0; tt < 8 ; tt++) {
						cerr << reader_indicator[tt] << ",";
						
					}
					
					cerr << endl;
						
						
				}
				exit(1);
			}
			
			
			#endif
			//msleep(1); // TODO more often, or use individual locks !? !?
			//nsleep(100000); // no sleep seems to perform better (on small test set)
			//sleep(1);
		}
		#pragma omp atomic
		reader_counter++;
		
		#ifdef DEBUG_DEADLOCK_NO
		#pragma omp critical(reader_indicator)
		{
		reader_indicator[omp_get_thread_num()]=1;
		}
		#endif
		omp_unset_lock(&writer_lock);
		
		
		flush_if_needed();
		
		//#pragma omp flush // TODO remove !?
	}
	
	void unset_reader_lock() {
		#pragma omp atomic
		reader_counter--;
		
		#ifdef DEBUG_DEADLOCK
		if (reader_counter < 0) {
			cerr << "error: reader_counter < 0" << endl;
			exit(1);
		}
		#endif
		
		#ifdef DEBUG_DEADLOCK_NO
		#pragma omp critical(reader_indicator)
		{
			reader_indicator[omp_get_thread_num()]=0;
		}
		#endif
	}
	
	int test_writer_lock() {
		return omp_test_lock(&writer_lock);
	}
	
	void set_writer_lock() {
		//cerr<< "omp_set_lock: in HAT " <<  omp_get_thread_num() << endl;
		
		#ifdef DEBUG_DEADLOCK
		omp_set_lock_with_timeout(&writer_lock);
		#else
		omp_set_lock(&writer_lock);
		#endif
		
		#ifdef DEBUG
		#pragma omp critical(writer_thread)
		{
		writer_thread = omp_get_thread_num();
			//#pragma omp flush
		}
		#endif
		//cerr<< "ok,omp_set_lock: in HAT " <<  omp_get_thread_num() << endl;
		//check if I need to flush
		
		flush_if_needed();
		
		
		// wait for readers
		#ifdef DEBUG_DEADLOCK
		timespec time_start;
		timespec time_end;
		clock_gettime(CLOCK_REALTIME, &time_start);
		#endif
		// wait for other process to finish their tasks...
		
		
		//int blabla = omp_get_num_threads();
		#pragma omp flush(reader_counter)
		while (reader_counter > 0) {
			
			
			//if (reader_counter >= blabla) {
			//	cerr << "error: reader_counter= " << reader_counter << endl;
			//	exit(1);
			//}
			
			#ifdef DEBUG_DEADLOCK
			clock_gettime(CLOCK_REALTIME, &time_end);
			if ( diff(time_start, time_end).tv_sec > 30 ){
				cerr << "thread: " << (int) omp_get_thread_num() << " timeout in hash table, waiting for reader! " << endl;
				
				exit(1);
			}
			#endif
			//if (cluster_kmer_hash_reader == 1) {
			//	break;
			//}
			//cerr << "thread: " << (int) omp_get_thread_num() << " wait for other readers... "<< cluster_kmer_hash_reader << endl;
			//msleep(1);
			//nsleep(100000);
			#pragma omp flush(reader_counter)
		}
	#ifdef DEBUG
		clock_gettime(CLOCK_REALTIME, &time_end);
		int waited = diff(time_start, time_end).tv_sec;
		
		if (waited > 5) {
			cerr << "waited for read lock: " << waited << endl;
		}
		
	#endif
		
		
		return;
	}
	
	void unset_writer_lock() {
		#ifdef USE_FLUSH_NEEDED
		setAllFlushBits();
		#endif
		#pragma omp flush
#ifdef DEBUG
#pragma omp critical(writer_thread)
		{
			if (writer_thread != omp_get_thread_num()) {
				cerr << "writer_thread: " << writer_thread << endl;
				cerr << "omp_get_thread_num(): " << omp_get_thread_num() << endl;
			}
		}
#endif
		omp_unset_lock(&writer_lock);
		return;
	}
	
	// notify other threads that they need to flush
	#ifdef USE_FLUSH_NEEDED
	void setAllFlushBits() {
		#pragma omp critical(flushbit)
		{
			flush_needed.set();
			flush_needed[omp_get_thread_num()] = 0; //this thread does not need to flush now!
		}
	}
	#endif
};


class ScopedLock
{
	private:
		omp_lock_t * lock;
	
	public:
		//ScopedLock(): lock(NULL){};
	
		ScopedLock(omp_lock_t * lock){
			this->lock = lock;
			#ifdef DEBUG_DEADLOCK
			omp_set_lock_with_timeout(lock);
			#else
			omp_set_lock(lock);
			#endif
		}
	
	~ScopedLock() {
		if (lock != NULL) {
			omp_unset_lock(lock);
		}
	}
};

#ifdef DEBUG_DEADLOCK
void omp_set_lock_with_timeout(omp_lock_t * lock) {
	timespec time_start;
	timespec time_end;
	clock_gettime(CLOCK_REALTIME, &time_start);

	
	while (not omp_test_lock(lock)) {
		nsleep(100000);
		clock_gettime(CLOCK_REALTIME, &time_end);
		
		if ( diff(time_start, time_end).tv_sec > 5 ) {
			#pragma omp critical(cerr)
			{
				cerr << "thread: " << (int) omp_get_thread_num() << " timeout in getting lock! " << endl;
			
								
				for (int tt = 0; tt < thread_count ; tt++) {
					cerr << thread_status[tt] << ",";
					
				}
				cerr << endl;
				for (int tt = 0; tt < thread_count ; tt++) {
					cerr << reads_per_thread[tt] << ",";
					
				}
				cerr << endl;
			
				
				
			} // end pragma
			exit(1);
				
		} // end if
	} // end while
	
}
#endif

#endif

