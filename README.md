# FMM2D
Classic Fast Multipole Method for 2 Dimensions

## Purpose of Software
The FMM2D repository is an academic serial code of the Fast Multipole Method in two dimensions.  All target and source points are located in the unit square.  The main purpose of the repository is to obtain an understanding of the Fast Multipole Method by way of a documented, object-oriented C++ code.  

The repository code is based on the 2005 Master thesis code by Yang Wang at UMIACS (University of Maryland Advanced Computer Studies). The thesis was directed by Dr. Ramani Duraiswami and it can be found towards the bottom of Dr. Nail Gumerov's personal web page http://www.umiacs.umd.edu/~gumerov/.  Accompanying the thesis is a link to the open source and object-oriented Java code written by Yang Wang.   

As mentioned, the FMM2D repository code is object-oriented and written in C++.  Modifications to Yang Wang's serial code have been made that avoid a bottleneck mentioned in his thesis.  Yang Wang mentions code slowdown due to Java's "garbage collection".  The FMM2D repository avoids this problems by not using dynamic memory allocation in the C++ code.  The same object classes used by Yang Wang are also seen in the C++ code and the codes are effectively the same.  One major difference between the two codes is the amount of documentation, using Doxygen and code comments, for ease of understanding.  The FMM2D repository also documents the mathematical derivations of the formulas used in both codes.  These mathematical details provide a supplement to Yang Wang's thesis, and an effort has been made to utilize the same symbolic notations.   

## FMM2D version 1.0.0

## Requirements
* Linux (Ubuntu Linux 14.04)
* g++   (g++ 4.8.4)
* Doxygen (doxygen 1.8.6)

The FMM2D repository code has only been tested on Ubuntu Linux 14.04, g++ version 4.8.4, and doxygen 1.8.6.  However, running on a linux machine with the software listed hopefully does not see too much trouble.  The make file can be run in a Bash shell terminal.

## Setup
The directory structure for the repository code is shown below.  Running 'make -f Makefile' in a Bash terminal of the working directory (top directory where the makefile Makefile is located) will compile the code and create the executable.  The command will build the object files, dependencies and executable and place them in the directory build/, deps/, and bin/, respectively. After compiling, run the program by going to the directory bin/ where the executable 'hello' is located and type ./hello in the terminal.

* Makefile
* Doxyfile
* customdoxygen.css
* bilbio.bib
* contributors.txt
* README.md
* src/
  * Main.cc 
  * FmmTree.cc
  * Box.cc
  * Point.cc
  * Potential.cc
  * Util.cc
  * Example1.cc
* include/
  * Main.h 
  * FmmTree.h
  * Box.h
  * Point.h
  * Potential.h
  * Util.h
  * Example1.h
* docs/
* deps/
* build/
* bin/
* test/

## Details

### Obtaining the Documentation
The repository can be imported in Eclipse and converted to a C/C++ project.  The code was written using Eclipse 3.8.  However a Bash terminal works fine.  Entering and returning 'make doc -f Makefile' will generate the Doxygen documentation (found in the top portion of Main.cc).  The index.html file for viewing the html version of the documentation can be found in the directory docs/html/.  However, Doxygen generates many files in that directory and therefore a .PHONY commmand called 'htmlIndex' has been added to the makefile Makefile.  To make finding the file easier, run 'make htmlIndex -f Makefile' in a terminal in the working directory.  This results in a symbolic link to index.html in the docs/html/ folder being placed in the working directory - called htmlIndex.  The html documentation for the code can then be opened by opening the htmlIndex link with a web browser.

### Member Functions Explained
The details to how certain member functions of the classes in the repository work can be found in the comments above the particular member functions of those classes.  For example, getting the indices of certain boxes (cells) of the FMM tree utilizes bit manipulation.  An explanation of the bit manipulation used in those member functions as well explanations of the member functions are given (by way of examples) in the code commments (above and sometimes in the member function code). 

### Test Example
Only one test example is provided.  The example is instantiated in Main.cc by a call made to class Example1 located in Example1.cc.  Instantiations of Example1 have the same source and target points and the class was created for this purpose. 
Example1 objects can be modified to increase (the refinement level and therefore) the number of source and target points that are created.  Classes to test examples that do not have source points that are the same as the target points can also be created.  The code is also written to handle these examples in the unit square.     
