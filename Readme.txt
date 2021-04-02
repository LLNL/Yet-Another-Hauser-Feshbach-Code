Readme to the Hauser-Feshbach code YAHFC-MC (A poor choice of name thus 
far, but it is "Yet Another hauser-Feshbach Code - Monte Carlo")

Here are the build instructions:

YAHFC-MC currently only supports the optical model code FRESCO. 

Here are the primary steps needed to build YAHFC-MC using FRESCO.

The basic steps: 


1] You need to set the following environment variables. In bash, you can do 
this by including the following lines in your .bashrc file.

#---------   YAHFC enviornment variables  ---------------------

YAHFC_DIR=$HOME/Projects/YAHFC
export YAHFC_DIR

YAHFC_DATA=$YAHFC_DIR/Data
export YAHFC_DATA

YAHFC_BIN_DIR=$HOME/bin
export YAHFC_BIN_DIR

PATH=$PATH:.:$YAHFC_BIN_DIR
export PATH

#---------   YAHFC enviornment variables  ---------------------


Here, $YAHFC_DIR is the parent directory where you will install all the
files and $YAHFC_BIN_DIR is a 'bin' directory where the executables will go, 
which you are free to choose for yourself. The directory $YAHFC_DATA is where 
the data files reside. The last line appends $YAHFC_BIN_DIR to your path, 
which you don't need if YAHFC_BIN_DIR already is.

Note you may have to 'source' your .bashrc file to get the enviornment 
variables set after you edit the .bashrc file.

2] Go to the directory $YAHFC_DIR/Src and type 'make YAHFC'.
This builds and places YAHFC.x into the $YAHFC_BIN_DIR directory.

3] Go to the directory where FRESCO is installed. Change directories to source,
and type "mk compiler", where compiler is a compiler name, e.g., gfortran. You
will need to change directories to the complier, i.e., "cd gfortran", then copy
frescox and sfrescox to a directory in your path, for example YAHFC_BIN_DIR.

4] You should create a run directory in a convenient place. You might want to
create a separate area since these directories can become quite large. 
There are several examples to try. These are in the directory 
$YAHFC_DIR/Runs. There are neutron reaction calculations for U238, U235, Pu239, 
Am241, and Y87. There is also a population example in Y88/Pop.
$YAHFC_DIR/Runs/U238. Start with the U238 example, and examine the .com file.
Be sure to set these to executable with chmod +x.

See the file Commands.txt for a brief description of the commands you can
give it.
