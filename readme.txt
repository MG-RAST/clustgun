

clustgun


INSTALLATION
------------
If you have admin rights on your system, you won't need to use --prefix option.

1) Sparsehash Package (formerly Google Sparsehash) (http://code.google.com/p/sparsehash/)
> wget http://sparsehash.googlecode.com/files/sparsehash-2.0.2.tar.gz
> tar xvfz sparsehash-2.0.2.tar.gz
> cd sparsehash-2.0.2
> ./configure --prefix=$HOME
> make
> make install

2) Boost C++ libraries (http://www.boost.org)
(you might also be able to install it from your systems software repository, e.g. ubuntu:  "sudo apt-get install libboost-dev")
> wget http://sourceforge.net/projects/boost/files/boost/1.50.0/boost_1_50_0.tar.gz
> tar xvfz boost_1_50_0.tar.gz
> cd boost_1_50_0
> ./bootstrap.sh --with-libraries=iostreams,program_options,system --prefix=$HOME
> ./b2 install

3) clustgun (https://github.com/MG-RAST/clustgun/)
download:
> git clone git://github.com/MG-RAST/clustgun.git
> cd clustgun

compile clustgun:
if you installed the libraries as above:
> make PREFIX=${HOME}
if your libraries are installed in systems default locations
> make
if your libraries are somewhere else:
> make INCLUDES='-I${HOME}/local/include' LFLAGS='-L${HOME}/local/lib ' LDFLAGS='-Wl,-R ${HOME}/local/lib'
and optionally you may want use macros like these:
> make CFLAGS='-DEBUG -DTIME'



old:
BOOSTLIB=$HOME/local ; g++ main.cpp kmer_iterator.cpp fasta_parser.cpp binarypath.cpp -o clustgun -lboost_iostreams -lboost_system  -lboost_program_options -I$BOOSTLIB/include -L$BOOSTLIB/lib  -m64 -Wl,-R $BOOSTLIB/lib -O3 -DDEBUG -DTIME

RUN
---
basic usage:
> ./clustgun input.fas

for help call clustgun without any parameters:
> ./clustgun


AUTHOR
------
Wolfgang Gerlach
xwgerlach@mcs.anl.gov (remove x in front of my name)



