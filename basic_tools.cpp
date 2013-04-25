

#include "basic_tools.hpp"





timespec diff(timespec start, timespec end)
{
	timespec temp;
	if ((end.tv_nsec-start.tv_nsec)<0) {
		temp.tv_sec = end.tv_sec-start.tv_sec-1;
		temp.tv_nsec = 1000000000+end.tv_nsec-start.tv_nsec;
	} else {
		temp.tv_sec = end.tv_sec-start.tv_sec;
		temp.tv_nsec = end.tv_nsec-start.tv_nsec;
	}
	return temp;
}

timespec add_time(timespec time1, timespec time2) {
	timespec result;
	result.tv_sec = time1.tv_sec + time2.tv_sec ;
    result.tv_nsec = time1.tv_nsec + time2.tv_nsec ;
    if (result.tv_nsec >= 1000000000L) {		/* Carry? */
        result.tv_sec++ ;  result.tv_nsec = result.tv_nsec - 1000000000L ;
    }
	
    return (result) ;
	
}


int msleep(unsigned long milisec)
{
    struct timespec req={0};
    time_t sec=(int)(milisec/1000);
    milisec=milisec-(sec*1000);
    req.tv_sec=sec;
    req.tv_nsec=milisec*1000000L;
    while(nanosleep(&req,&req)==-1)
		continue;
    return 1;
}

int nsleep(unsigned long nanosec) // 1.000.000 milisec = 1 nanosec !!!
{
    struct timespec req={0};
    //time_t sec=(int)(milisec/1000);
    //milisec=milisec-(sec*1000);
    req.tv_sec=0;
    req.tv_nsec=nanosec;
    while(nanosleep(&req,&req)==-1)
		continue;
    return 1;
}


int * get_nucleobase_ascii2num() {
	
	int * nucleobase_ascii2num;
	nucleobase_ascii2num = new int [256];
	
	for (int i=0; i< 256; i++) {nucleobase_ascii2num[i] = -1;}
	nucleobase_ascii2num['A']=0;
	nucleobase_ascii2num['T']=1;
	nucleobase_ascii2num['G']=2;
	nucleobase_ascii2num['C']=3;
	nucleobase_ascii2num['a']=0;
	nucleobase_ascii2num['t']=1;
	nucleobase_ascii2num['g']=2;
	nucleobase_ascii2num['c']=3;
	nucleobase_ascii2num['U']=1;
	nucleobase_ascii2num['u']=1;
	
	return nucleobase_ascii2num;
}

std::string getFileExtension(const std::string& FileName) {
	//cout <<"Z:" <<FileName << endl;
	size_t pos  = FileName.find_last_of(".");
    if(pos != std::string::npos) {
        return FileName.substr(pos+1);
	}
    return "";
}


std::string getFileNameWithoutExtension(const std::string& FileName) {
	size_t start  = FileName.find_last_of("/");
	
    if(start == std::string::npos) {
        start = 0;
	} else {
		start++;
	}
	
	size_t end = FileName.find_last_of(".");
	if(end != std::string::npos) {
		if (end < start) {
			end = std::string::npos;
		}
	}
	
	if(end == std::string::npos) {
		return FileName.substr(start);
	} 
	
    return FileName.substr(start, end-start);
}


void process_mem_usage(double& vm_usage, double& resident_set)
{
	// source: http://stackoverflow.com/questions/669438/how-to-get-memory-usage-at-run-time-in-c
	// usage: 
	// double vm, rss;
	// process_mem_usage(vm, rss);
	// cout << "VM: " << vm << "; RSS: " << rss << endl;
#ifdef __APPLE__
	vm_usage     = 0;
	resident_set = 0;
	return;
#endif
	
	using std::ios_base;
	using std::ifstream;
	using std::string;
	
	vm_usage     = 0.0;
	resident_set = 0.0;
	
	// 'file' stat seems to give the most reliable results
	//
	ifstream stat_stream("/proc/self/stat",ios_base::in);
	
	if (stat_stream.good() == false) {
		cerr << "warning: /proc/self/stat not found, not Linux ?" << endl;
		
		return;
		
	}
	
	// dummy vars for leading entries in stat that we don't care about
	//
	string pid, comm, state, ppid, pgrp, session, tty_nr;
	string tpgid, flags, minflt, cminflt, majflt, cmajflt;
	string utime, stime, cutime, cstime, priority, nice;
	string O, itrealvalue, starttime;
	
	// the two fields we want
	//
	unsigned long vsize;
	long rss;
	
	stat_stream >> pid >> comm >> state >> ppid >> pgrp >> session >> tty_nr
	>> tpgid >> flags >> minflt >> cminflt >> majflt >> cmajflt
	>> utime >> stime >> cutime >> cstime >> priority >> nice
	>> O >> itrealvalue >> starttime >> vsize >> rss; // don't care about the rest
	
	stat_stream.close();
	
	long page_size_kb = sysconf(_SC_PAGE_SIZE) / 1024; // in case x86-64 is configured to use 2MB pages
	vm_usage     = vsize / 1024.0;
	resident_set = rss * page_size_kb;
}


std::string exec(const char* cmd) {
    FILE* pipe = popen(cmd, "r");
    if (!pipe) return "ERROR";
    char buffer[128];
    std::string result = "";
    while(!feof(pipe)) {
        if(fgets(buffer, 128, pipe) != NULL)
			result += buffer;
    }
    pclose(pipe);
    return result;
}


double round( double value )
{
	return floor( value + 0.5 );
}

string stream2string(std::stringstream & stream) {
	std::string s;
	std::ostringstream os;
	os<<stream.rdbuf();
	s=os.str();
	return s;
}

bool FileExists(const char * filename) {
	struct stat stFileInfo;
	bool blnReturn;
	int intStat;
	
	// Attempt to get the file attributes
	intStat = stat(filename,&stFileInfo);
	if(intStat == 0) {
		// We were able to get the file attributes
		// so the file obviously exists.
		blnReturn = true;
	} else {
		// We were not able to get the file attributes.
		// This may mean that we don't have permission to
		// access the folder which contains this file. If you
		// need to do that level of checking, lookup the
		// return values of stat which will give you
		// more details on why stat failed.
		blnReturn = false;
	}
	
	return(blnReturn);
}


void checkFile(string * file) {
	if (file == 0) {
		return;
	}
	
	if (! FileExists((char *) file->c_str())) {
		cerr << "error: cannot read file: " << *file << endl;
		exit(EXIT_FAILURE);
	}
	
}

string int2str (int n) {
	stringstream ss;
	ss << n;
	return ss.str();
}

int str2int(string& text){
	
	int number = std::atoi( text.c_str() );
	
	return number;
	
}

void reverseComplementDNA(string& sequence) {
	int len = sequence.length();
	
	// A-T
	// G-C
	
	// reverse:
	int half = len/2;
	for (int i=0; i<half; i++) {
		swap(sequence[i], sequence[len-1-i]);
	}
	
	// complement:
	char c;
	for (int i=0; i<len; i++) {
		c = sequence[i];
		c = toupper(c);
		switch (c) {
			case 'A':
				c='T';
				break;
			case 'C':
				c='G';
				break;
			case 'G':
				c='C';
				break;
			case 'T':
				c='A';
				break;
			case 'U':
				c='C';
				break;
			default:
				//cerr << "error at " << i << " with " << sequence << endl;
				//exit(1);
				c='N';
		}
		sequence[i]=c;
	}
	
	
}
