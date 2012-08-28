
# which OS
UNAME := $(shell uname)

CC=g++
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
MACRO=


all: $(SOURCES) $(EXECUTABLE)

clean :
	rm $(EXECUTABLE) $(OBJECTS)

debug: CFLAGS += -DDEBUG
debug: all
	
$(EXECUTABLE): $(OBJECTS) 
	$(CC) -o $(EXECUTABLE)  $(OBJECTS)  $(LIBS) $(LFLAGS) $(LDFLAGS)

.cpp.o:
	$(CC) -c $(CFLAGS) $(INCLUDES) $(MACRO) $< -o $@  
#$(LIBS) 


#g++ -o clustgun main.cpp kmer_iterator.cpp fasta_parser.cpp binarypath.cpp -lboost_iostreams-mt -lboost_system-mt  -lboost_program_options-mt -O3 -DDEBUG -DTIME
