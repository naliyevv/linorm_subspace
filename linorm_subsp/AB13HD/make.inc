#  Modify the FORTRAN and OPTS definitions to refer to the
#  compiler and desired compiler options for your machine.  NOOPT
#  refers to the compiler options desired when NO OPTIMIZATION is
#  selected.  Define LOADER and LOADOPTS to refer to the loader and 
#  desired load options for your machine.
#
#  IMPORTANT: Use the options -fPIC -fdefault-integer-8 for 64bit
#             architectures.
#
FORTRAN  = gfortran
OPTS     = -O2 -fPIC -fdefault-integer-8
NOOPT    = -O0 -fPIC -fdefault-integer-8
