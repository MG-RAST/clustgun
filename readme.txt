

clustgun


INSTALLATION
------------
If you have admin rights on your system, you won't need to use --prefix option.

1) Boost C++ libraries (http://www.boost.org)
in ubuntu:  "sudo apt-get install libboost-system-dev libboost-iostreams-dev libboost-program-options-dev"
or from source:
> wget http://sourceforge.net/projects/boost/files/boost/1.53.0/boost_1_53_0.tar.gz
> tar xvfz boost_1_53_0.tar.gz
> cd boost_1_53_0
> ./bootstrap.sh --with-libraries=iostreams,program_options,system --prefix=$HOME
> ./b2 --layout=system install

2) clustgun (https://github.com/MG-RAST/clustgun/)
download:
> git clone git://github.com/MG-RAST/clustgun.git
> cd clustgun

compile clustgun:
if you installed the libraries as above with --prefix:
> make PREFIX=${HOME}

or, if your libraries are installed in systems default locations
> make

or, if your libraries are somewhere else:
> make INCLUDES='-I${HOME}/local/include' LFLAGS='-L${HOME}/local/lib ' LDFLAGS='-Wl,-R ${HOME}/local/lib'
alternatively, if you need to compile more often, you may want to set the following environment variables within your .bashrc:
----------------------
CPLUS_INCLUDE_PATH=${HOME}/local/include/:$CPLUS_INCLUDE_PATH
export CPLUS_INCLUDE_PATH

LIBRARY_PATH=${HOME}/local/lib/:$LIBRARY_PATH
export LIBRARY_PATH

LD_LIBRARY_PATH=${HOME}/local/lib/:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH
----------------------

If the BOOST library is not found (although installed), check if the library file names are suffixed with "-mt" and change the makefile accordingly.

only for developing purposes:
> make CFLAGS='-DEBUG -DTIME -DDEBUG_DEADLOCK'

RUN
---
basic usage:
> ./clustgun input.faa

for help call clustgun without any parameters:
> ./clustgun


AUTHOR
------
Wolfgang Gerlach
xwgerlach@mcs.anl.gov (remove x in front of my name)



