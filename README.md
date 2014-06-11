clustgun
========



Installation
------------

If you have admin rights on your system, you won't need to use --prefix option.

**1 Boost C++ libraries** (http://www.boost.org)

in ubuntu:
```bash
> sudo apt-get install libboost-system-dev libboost-iostreams-dev libboost-program-options-dev
```
or from source:
```bash
> wget http://sourceforge.net/projects/boost/files/boost/1.53.0/boost_1_53_0.tar.gz
> tar xvfz boost_1_53_0.tar.gz
> cd boost_1_53_0
> ./bootstrap.sh --with-libraries=iostreams,program_options,system --prefix=$HOME
> ./b2 --layout=system install
```
**2 clustgun** (https://github.com/MG-RAST/clustgun/)
download:
```bash
> git clone git://github.com/MG-RAST/clustgun.git
> cd clustgun
```
compile clustgun:
ubuntu: (or if your libraries are installed in systems default locations)
```bash
> make
```
or, if you installed the libraries as above with --prefix:
```bash
> make PREFIX=${HOME}
```
or, if your libraries are somewhere else:
```bash
> make INCLUDES='-I${HOME}/local/include' LFLAGS='-L${HOME}/local/lib ' LDFLAGS='-Wl,-R ${HOME}/local/lib'
```
alternatively, if you need to compile more often, you may want to set the following environment variables within your .bashrc:
```bash
CPLUS_INCLUDE_PATH=${HOME}/local/include/:$CPLUS_INCLUDE_PATH
export CPLUS_INCLUDE_PATH

LIBRARY_PATH=${HOME}/local/lib/:$LIBRARY_PATH
export LIBRARY_PATH

LD_LIBRARY_PATH=${HOME}/local/lib/:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH
```


If the BOOST library is not found (although installed on your system), one reason can be that the library file names are suffixed with "-mt" or something else. In that case you need to modify the BOOST libary file names in the makefile accordingly. The only way for me to detect the correct BOOST file names would be to use Autotools, but for that I was too lazy.

only for developers:
```bash
> make CFLAGS='-DDEBUG -DTIME -DDEBUG_DEADLOCK'
```

**3 FragGeneScan**

```bash
> git clone git://github.com/wltrimbl/FGS.git
> cd FGS
```
In case you need the exact version of FragGeneScan that was used in the evaluation of clustgun, use
```bash
> git reset --hard e1aedead8512801f4324a049caa2f38eebbdfe8d"
```
and then
```bash
> make
```

**4 GNU Parallel**

To run FGS in parallel using multiple CPU cores you can use the fgs.sh wrapper script that is included in the clustgun directory. This script uses a recent version of GNU Parallel.
Under ubuntu you might first need to remove an older version of GNU Parallel:
```bash
> sudo apt-get remove moreutils
```
To install GNU Parallel on a debian system (like ubuntu) use:

```bash
> wget http://launchpadlibrarian.net/104233941/parallel_20120422-1_all.deb && sudo dpkg -i ./parallel_20120422-1_all.deb
```

The version reported by   

```bash
> parallel --minversion 1 
```

should be "20120422".



Run FragGeneScan
----------------
Copy fgs.sh into the FGS directory.
To run FGS with 4 CPUs call from within the FGS directory:

```bash
> ./fgs.sh 4 input.fasta
```

This should create a result.faa. To change the parameters of FGS please modify the fgs.sh script with a text editor accordingly. Note that for IO-efficiency this way of running FGS will give you only the protein fragment files (.faa), not the other output files that FGS produces by default.  


Run Clustgun
------------
basic usage:
```bash
> ./clustgun input.faa
```

for help call clustgun without any parameters:
```bash
> ./clustgun
```

Test Dataset
-------------
Protein fragments predicted with FragGeneScan from a soil sample:
* 10k fragments: [4477899.3.fgs.10k.faa.gz](http://shock.metagenomics.anl.gov:80/node/9f6ddf16-e297-42fc-b3a4-808c025c20bc?download&filename=4477899.3.fgs.10k.faa.gz)
* 100k fragments: [4477899.3.fgs.100k.faa.gz](http://shock.metagenomics.anl.gov:80/node/1f17b284-f046-41c3-b558-2ad378714b52?download&filename=4477899.3.fgs.100k.faa.gz)
* 1m fragments: [4477899.3.fgs.1m.faa.gz](http://shock.metagenomics.anl.gov:80/node/1c952b5e-b538-471e-b772-22c2ac878304?download&filename=4477899.3.fgs.1m.faa.gz)

Example:
```bash
> ./clustgun --threads=2 4477899.3.fgs.10k.faa.gz
```

Author
------
Wolfgang Gerlach
xwgerlach@mcs.anl.gov (remove x in front of my name)
Please feel free to contact me if you have issues installing or running clustgun.


