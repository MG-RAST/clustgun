
# which OS
UNAME := $(shell uname)

CXX=g++

# gwich gcc version
GCC_MAJOR    := $(shell $(CXX) -dumpversion 2>&1 | cut -d'.' -f1)
GCC_MINOR    := $(shell $(CXX) -dumpversion 2>&1 | cut -d'.' -f2)

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


ifeq ($(UNAME), Linux)
# Linux
override LIBS += -lboost_iostreams -lboost_system -lboost_program_options
else
# OSX, e.g. Darwin
override LIBS += -lboost_iostreams-mt -lboost_system-mt -lboost_program_options-mt
endif

SOURCES=main.cpp kmer_iterator.cpp fasta_parser.cpp binarypath.cpp
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=clustgun



all: $(SOURCES) $(EXECUTABLE)

clean :
	rm -f $(EXECUTABLE) $(OBJECTS)

debug: CFLAGS += -DDEBUG
debug: all
	
$(EXECUTABLE): $(OBJECTS) 
	$(CXX) -o $(EXECUTABLE)  $(OBJECTS)  $(LIBS) $(LFLAGS) $(LDFLAGS)

.cpp.o:
	$(CXX) -c $(CFLAGS) $(INCLUDES) $< -o $@  

#g++ -o clustgun main.cpp kmer_iterator.cpp fasta_parser.cpp binarypath.cpp -lboost_iostreams-mt -lboost_system-mt  -lboost_program_options-mt -O3 -DDEBUG -DTIME
