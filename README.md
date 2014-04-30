
clustgun
========



installation
------------

If you have admin rights on your system, you won't need to use --prefix option.

**1 Boost C++ libraries** (http://www.boost.org)
ubuntu:
```
> sudo apt-get install libboost-system-dev libboost-iostreams-dev libboost-program-options-dev
```
or from source:
```
> wget http://sourceforge.net/projects/boost/files/boost/1.53.0/boost_1_53_0.tar.gz
> tar xvfz boost_1_53_0.tar.gz
> cd boost_1_53_0
> ./bootstrap.sh --with-libraries=iostreams,program_options,system --prefix=$HOME
> ./b2 --layout=system install
```
**2 clustgun** (https://github.com/MG-RAST/clustgun/)
download:
```
> git clone git://github.com/MG-RAST/clustgun.git
> cd clustgun
```
compile clustgun:
ubuntu: (or if your libraries are installed in systems default locations)
```
> make
```
or, if you installed the libraries as above with --prefix:
```
> make PREFIX=${HOME}
```
or, if your libraries are somewhere else:
```
> make INCLUDES='-I${HOME}/local/include' LFLAGS='-L${HOME}/local/lib ' LDFLAGS='-Wl,-R ${HOME}/local/lib'
```
alternatively, if you need to compile more often, you may want to set the following environment variables within your .bashrc:
```
CPLUS_INCLUDE_PATH=${HOME}/local/include/:$CPLUS_INCLUDE_PATH
export CPLUS_INCLUDE_PATH

LIBRARY_PATH=${HOME}/local/lib/:$LIBRARY_PATH
export LIBRARY_PATH

LD_LIBRARY_PATH=${HOME}/local/lib/:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH
```


If the BOOST library is not found (although installed on your system), one reason can be that the library file names are suffixed with "-mt" or something else. In that case you need to modify the BOOST libary file names in the makefile accordingly. The only way for me to detect the correct BOOST file names would be to use Autotools, but for that I was too lazy.

only for developers:
```
> make CFLAGS='-DDEBUG -DTIME -DDEBUG_DEADLOCK'
```

**3 FragGeneScan**

```
> git clone git://github.com/wltrimbl/FGS.git
> cd FGS
```
In case you need the exact version of FragGeneScan that was used in the evaluation of clustgun, use "git reset --hard e1aedead8512801f4324a049caa2f38eebbdfe8d"
```
> make
```

**4 GNU Parallel**

To run FGS in parallel using multiple CPU cores you can use the fgs.sh wrapper script that is included in the clustgun directory. This script uses a recent version of GNU Parallel.
Under ubuntu you might first need to remove an older version of GNU Parallel:
```
> sudo apt-get remove moreutils
```
To install GNU Parallel on a debian system (like ubuntu) use:
```
> wget http://launchpadlibrarian.net/104233941/parallel_20120422-1_all.deb && sudo dpkg -i ./parallel_20120422-1_all.deb
```
The version reported by   
```
> parallel --minversion 1 
```
should be "20120422".



run FragGeneScan
----------------
Copy fgs.sh into the FGS directory.
To run FGS with 4 CPUs call from within the FGS directory:
```
> ./fgs.sh 4 input.fasta
```
This should create a result.faa

To change the parameters of FGS please modify the fgs.sh script with a text editor accordingly. Note that for IO-efficiency this way of running FGS will give you only the protein fragment files (.faa), not the other output files that FGS produces by default.  


run clustgun
------------
basic usage:
```
> ./clustgun input.faa
```

for help call clustgun without any parameters:
```
> ./clustgun
```



author
------
Wolfgang Gerlach
xwgerlach@mcs.anl.gov (remove x in front of my name)
Please feel free to contact me if you have issues installing or running clustgun.


