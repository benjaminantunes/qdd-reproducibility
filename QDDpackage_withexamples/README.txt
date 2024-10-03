This collection contains the following of supplementary material:

sub-directory 'doc': 
  Contains the user manual 'User_manual.pdf' which explains input,
  output, and handling of the code in details.

sub-directory 'src/qdd':
  Contains the sources for the code with the makfiles.

sub-directory 'src/auxiliary':
  Contains a few auxiliary routines for post processing. The handling
  is explained in the user manual.

sub/directory 'bin':
  This is an empty sub-directory, provided to hold the executable
  after compiling and linking.

sub-directory 'examples':
  Contains input and releveant output of all examples discussed in
  paper and manual. The collection serves also as becnhmark for first
  own tests. Note that the output files are zipped with the 'xz'
  utility (which you may have to install on some systems). To unzip 
  use 'unxz'.


All experiments are located in /bin folder. 
There is 50 replications (measurment of time and consumption) for each experiment.
There is a jupyter notebook to guide you in the workflow.

The compilation is located in src/qdd, with makefile. This is where we compiled each experiment.

You need to install gfortran, ifort and ifx compiler.
You also need to install FFTW3 and MKL.

All are free.