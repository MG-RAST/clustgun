

clustgun


INSTALL
-------
google hash map:
http://goog-sparsehash.sourceforge.net/

boost lib:
http://www.boost.org (you might also be able to install it from your systems software repository)
./bootstrap.sh --with-libraries=iostreams,program_options,system --prefix=/homes/wgerlach/
./b2 install

COMPILE
-------
Linux/UNIX:
BOOSTLIB=/homes/wgerlach/local ; g++ main.cpp kmer_iterator.cpp fasta_parser.cpp -o clustgun -lboost_iostreams -lboost_system  -lboost_program_options -I$BOOSTLIB/include -L$BOOSTLIB/lib  -m64 -Wl,-R $BOOSTLIB/lib -DDEBUG -O3

OSX:
g++ main.cpp kmer_iterator.cpp fasta_parser.cpp -o clustgun -lboost_iostreams-mt -lboost_system-mt  -lboost_program_options-mt -DDEBUG -O3


RUN
---
./clustgun input.fas
./clustgun --help

