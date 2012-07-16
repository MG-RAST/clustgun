

clustgun


INSTALLATION
------------
If you have admin rights on your system, you won't need to use --prefix option.

1) Sparsehash Package (formerly Google Sparsehash) (http://code.google.com/p/sparsehash/)
wget http://sparsehash.googlecode.com/files/sparsehash-2.0.2.tar.gz
tar xvfz sparsehash-2.0.2.tar.gz
cd sparsehash-2.0.2
./configure --prefix=$HOME
make

2) Boost C++ libraries (http://www.boost.org)
(you might also be able to install it from your systems software repository)
wget http://sourceforge.net/projects/boost/files/boost/1.50.0/boost_1_50_0.tar.gz/download
tar xvfz boost_1_50_0.tar.gz
cd boost_1_50_0
./bootstrap.sh --with-libraries=iostreams,program_options,system --prefix=$HOME
./b2 install

3) clustgun (https://github.com/MG-RAST/clustgun/)
git clone git://github.com/MG-RAST/clustgun.git
cd clustgun

Linux/UNIX:
BOOSTLIB=$HOME/local ; g++ main.cpp kmer_iterator.cpp fasta_parser.cpp -o clustgun -lboost_iostreams -lboost_system  -lboost_program_options -I$BOOSTLIB/include -L$BOOSTLIB/lib  -m64 -Wl,-R $BOOSTLIB/lib -O3 -DDEBUG

OSX:
g++ main.cpp kmer_iterator.cpp fasta_parser.cpp -o clustgun -lboost_iostreams-mt -lboost_system-mt  -lboost_program_options-mt -O3 -DDEBUG


RUN
---
./clustgun input.fas
./clustgun --help

AUTHOR
------
Wolfgang Gerlach
xwgerlach@mcs.anl.gov (remove x in front of my name)



