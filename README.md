This project has code to verify Oppermann's and Legendre's conjectures.

The code uses MPI and GMP packages, and the main program is in "parmain.cpp".

The bigdawg folder is for the code and data from the run on our linux cluster (192 cores).
The phi folder is for the code and data from the run on our four Intel phi coprocessors (256 cores).

Each folder contains a makefile for compiling the code. Use "make parmain".
Make sure you have an EMPTY folder called "data" in the folder you run the code from (parmain).
The program will benchmark data into files called id.out, where id is replaced by the 3-digit id number of each process.

If you want to change the range of the computation, you must edit the main loop in "parmain.cpp".
