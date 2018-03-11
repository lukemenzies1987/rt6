# rt6
A DEMO tailored rate theory code , used for nuclear fusion material modelling research

This Fortran90 code was used for work in within the thesis 'Modelling Helium Embrittlement in Iron Based Metals Under DEMO Conditions' , Univerity of Manchester. 

Please note, this code is not recommended to be used by anyone who isn't compitent with Fortran90. 

The code requires the nag-fortran library, as well as blas and lapack. The makefile is locate in the main directory. In the makefile, type the nag directory on your CPU in 'NAGFDIR' (e.g NAGFDIR   ='/opt/NAG/fll...').

It is recommended that you run the make file using the optimization flags no higher than 02, though 03 is available. 
Debugging flags are written into the make file if needed (uncomment them). Change the compiler to the desired compiler in the 'COMPILE' and 'LOAD' part. If the intel fortran compiler is the default compiler for your computer then ignore this. 

Compiling-----
In the directory, to compile the code, type 'make'. 
Running--------
If the code compiles successfully, run the RT6 file created. 
Input file-----
The input file is located in the directory named 'param.in'. Change the parameter if to your desired values. 
Output--------
All output files are located in the Output directory. 
