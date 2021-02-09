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
LOADER   = gfortran
LOADOPTS = $(SLICOTLIB) $(LPKAUXLIB) $(LAPACKLIB) $(BLASLIB)
#
#  The location of the libraries to which you will link.  (The 
#  machine-specific, optimized BLAS library should be used whenever
#  possible.)
#
BLASLIB      = ../libs/blas64.a
LAPACKLIB    = ../libs/lapack64.a
SLICOTLIB    = ../libs/slicot64.a
LPKAUXLIB    = ../libs/lpkaux64.a
