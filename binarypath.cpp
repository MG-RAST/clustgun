
#include "binarypath.hpp"



string getbinarypath() {

	string str;
	
	#ifdef __APPLE__
	char path[1024];
	char resolved_name[1024];
	
	char * real;
	
	uint32_t size = sizeof(path);
	if (_NSGetExecutablePath(path, &size) == 0) {
		//cout << "executable path is " << path << endl;
		
		real = realpath(path, resolved_name);
		if (real == NULL ) {
			cerr << "error converting to real filename" << endl;
			exit(1);
			
		}
		//cout << "resolved path is " << real << endl;
	}
	else {
		cerr << "getbinarypath(): buffer too small; need size\n" << size << endl;
		exit(1);
	}
	
	str = string(real);
	
	#else

	// ###Linux###
	
	char buffer[BUFSIZ];
	size_t foo = readlink("/proc/self/exe", buffer, BUFSIZ);
	printf("%s\n", buffer);

	
	str = string(buffer);
	#endif
	
	
	size_t found = str.find_last_of('/');
	if (found == string::npos) {
		cerr << "error: binary path does not contain \"/\" character: " << str << endl;
		exit(1);
	}
	str = str.substr(0,found);
	
	//cout << "go: " << str << endl;
	return str;
}

