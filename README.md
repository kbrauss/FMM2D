# FMM2D
Classic Fast Multipole Method for 2 Dimensions

# # Purpose of Software
The FMM2D repository is an academic serial code of the Fast Multipole Method in two dimensions.  The main purpose of the repository is to obtain an understanding of the Fast Multipole Method by way of a documented, object-oriented C++ code.  
The repository code is based on the 2005 Master thesis code by Yang Wang at UMIACS (University of Maryland Advanced Computer Studies). The thesis was directed by Dr. Ramani Duraiswami and it can be found towards the bottom of Dr. Nail Gumerov's personal web page http://www.umiacs.umd.edu/~gumerov/.  Accompanying the thesis is a link to the open source and object-oriented Java code written by Yang Wang.   

As mentioned, the FMM2D repository code is also object-oriented and written in C++.  Modifications to Yang Wang's serial code have been made that avoid a bottleneck mentioned in his thesis.  Yang Wang mentions code slowdown due to Java's "garbage collection".  The FMM2D repository avoids this problems by not using dynamic memory allocation in the C++ code.  The same object classes used by Yang Wang are also seen in the C++ code and the codes are effectively the same.  One major difference between the two codes is the amount of documentation, using Doxygen and code comments, for ease of understanding.  The FMM2D repository also documents the mathematical derivations of the formulas used in both codes.  These mathematical details provide a supplement to Yang Wang's thesis, and an effort has been made to utilize the same symbolic notations.   
