#
#    CONFIGURED TO USE *REPO* LIBRARIES FOR LAPACK AND BLAS (SEE AT THE BOTTOM)
#
# https://stackoverflow.com/questions/18136918/how-to-get-current-relative-directory-of-your-makefile
#
#CURRENT_DIR := $(dir $(realpath $(firstword $(MAKEFILE_LIST))))
#MKFILE_PATH := $(abspath $(firstword $(MAKEFILE_LIST)))
#CURRENT_DIR := $(notdir $(patsubst %/,%,$(dir $(MKFILE_PATH))))
#CURRENT_DIR := $(CURDIR)
# LAPACK_DIR is defined by the caller and passed to "make" as argument
#
############################################################################
#
#  Program:         ScaLAPACK
#
#  Module:          SLmake.inc
#
#  Purpose:         Top-level Definitions
#
#  Creation date:   February 15, 2000
#
#  Modified:        October 13, 2011
#
#  Send bug reports, comments or suggestions to scalapack@cs.utk.edu
#
############################################################################
#
#  C preprocessor definitions:  set CDEFS to one of the following:
#
#     -DNoChange (fortran subprogram names are lower case without any suffix)
#     -DUpCase   (fortran subprogram names are upper case without any suffix)
#     -DAdd_     (fortran subprogram names are lower case with "_" appended)

CDEFS         = -DAdd_

#
#  The fortran and C compilers, loaders, and their flags
#

#FC            = mpif90
#CC            = mpicc 
FC            = mpiifort
CC            = mpiicc 
NOOPT         = -O0
FCFLAGS       = -O2
CCFLAGS       = -O2
FCLOADER      = $(FC)
CCLOADER      = $(CC)
FCLOADFLAGS   = $(FCFLAGS)
CCLOADFLAGS   = $(CCFLAGS)

#
#  The archiver and the flag(s) to use when building archive (library)
#  Also the ranlib routine.  If your system has no ranlib, set RANLIB = echo
#

ARCH          = ar
ARCHFLAGS     = cr
RANLIB        = ranlib

#
#  The name of the ScaLAPACK library to be created
#

SCALAPACKLIB  = libscalapack.a

#
#  BLAS, LAPACK (and possibly other) libraries needed for linking test programs
#

#BLASLIB       = -lblas
#LAPACKLIB     = -llapack
BLASLIB       = $(LAPACK_DIR)/librefblas.a
LAPACKLIB     = $(LAPACK_DIR)/liblapack.a
LIBS          = -nofor-main $(LAPACKLIB) $(BLASLIB) -lgfortran -lifcore -ldl
