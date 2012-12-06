
# which OS
UNAME := $(shell uname)

CXX=g++

# test for >= 4.3
GCCVERSIONGTEQ43 := $(shell expr `gcc -dumpversion | cut -f1,2 -d.` \> 4.3)

ifeq "$(GCCVERSIONGTEQ43)" "1"
    override CFLAGS += -std=c++0x
endif

ifdef PREFIX
override INCLUDES += -I${PREFIX}/local/include
override LFLAGS += -L${PREFIX}/local/lib
override LDFLAGS += -Wl,-R ${PREFIX}/local/lib
endif

override CFLAGS += -O3 -m64
# override directive allows me to append values to variables from makefile arguments without overwriting
override LDFLAGS +=



override LIBS += -lboost_iostreams -lboost_system -lboost_program_options
#override LIBS += -lboost_iostreams-mt -lboost_system-mt -lboost_program_options-mt


SOURCES=main.cpp kmer_iterator.cpp fasta_parser.cpp binarypath.cpp read_blosum.cpp basic_tools.cpp
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=clustgun

all:  $(SOURCES) $(EXECUTABLE) 

clean :
	rm -f $(EXECUTABLE) $(OBJECTS) git_ref.h

debug: CFLAGS += -DDEBUG
debug: all
	
$(EXECUTABLE): git_ref.h $(OBJECTS) 
	$(CXX) -o $(EXECUTABLE)  $(OBJECTS)  $(LIBS) git_ref.h $(LFLAGS) $(LDFLAGS)
	rm -f git_ref.h

.cpp.o:
	$(CXX) -c $(CFLAGS) $(INCLUDES) $< -o $@  

git_ref.h:
	echo "#define GIT_REF \"`git rev-parse --verify --short HEAD 2>/dev/null || echo unknown`\"" > git_ref.h
	echo "#define GIT_DATE \"`git show -s --format="%ci" HEAD 2>/dev/null || echo unknown`\"" >> git_ref.h



#g++ -o clustgun main.cpp kmer_iterator.cpp fasta_parser.cpp binarypath.cpp -lboost_iostreams-mt -lboost_system-mt  -lboost_program_options-mt -O3 -DDEBUG -DTIME
