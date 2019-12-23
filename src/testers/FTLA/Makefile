
#FFLAGS +=  -fsecond-underscore
CFLAGS += -Wall 
CFLAGS += -g -O0
CFLAGS += -DINJECT
#CFLAGS += -DUSE_CoF
#CFLAGS += -DNO_CHECKPOINT
#CFLAGS += -DVERIFY_CHK
#CFLAGS += -DNO_QLOCALCOPY
#CFLAGS += -DNO_EXTRAFLOPS
#CFLAGS += -DTIMING
#CFLAGS += -DPRINTMATRIX

# Compiler selection based on architecture, should not have to change from there
ifneq (,${MKLROOT})
  CFLAGS += -I${MKLROOT}/include -I$(MKL_INC) 
  BLACSLIB = -L$(MKL_LIB) -lmkl_scalapack_lp64 -lmkl_intel_lp64 -lmkl_blacs_intelmpi_lp64 -lmkl_core -lmkl_sequential -lmkl_gf_lp64 -lm
endif
# Mac OS-X
ifneq (,$(findstring darwin,$(shell echo $$OSTYPE)))
  #CC = openmpicc
  #FC = openmpif77
  CC = mpicc
  FC = mpif77
  BLACSLIB ?= ${HOME}/DEVEL/FTLA/ftla/build/lib/libscalapack.a /opt/local/lib/libopenblas.a ${HOME}/DEVEL/FTLA/lapack-3.4.2/liblapack.a 
  #CFLAGS += -std=c99
  #FFLAGS += -cpp
  #BLACSLIB ?= ${HOME}/DEVEL/FTLA/ftla/build/lib/libscalapack.a /opt/local/lib/libopenblas.a ${HOME}/DEVEL/FTLA/lapack-3.4.2/liblapack.a 
else 
# Cray XT
ifneq (,$(findstring linux,$(shell echo $$XTPE_COMPILE_TARGET)))
  CC = cc
  FC = ftn
  FFLAGS += -Mnomain
else
  CC = mpicc
  FC = mpif77
  BLACSLIB ?= -lscalapack -llapack -lgotoblas2
endif
endif
#LDFLAGS =  -L$(HOME)/lib/GotoBLAS2   \
		   -L$(HOME)/lib/scalapack-1.8.0 \
		   -L$(HOME)/lib/lapack-3.3.0 \
		   -L$(HOME)/lib/BLACS/LIB 

#LDFLAGS="-L/usr/local/opt/llvm/lib -L/usr/local/opt/lapack/lib -L/usr/local/opt/openblas/lib"

all: libftla.a testing_ftdtr.x testing_ftdqr.x 
#testing_ftdqr.x 
#testing_ftdtr.x   

.SECONDEXPANSION:
HELPERS = psmatgen.o pdmatgen.o pcmatgen.o pzmatgen.o pmatgeninc.o 
dqr_ORIG = pdgeqrrv.o #pdgetrrv.o #pdgeqrf.o pdlarfb.o pdgetrf.o 
dtr_ORIG = pdgetrrv.o pdlaprnt.o pdgeqrf.o pdlarfb.o pdgetrf.o pdgetf2.o

SFTLA=ftla_scsum.o ftla_psgeqrf.o ftla_pslarfb.o ftla_psgetrf.o
DFTLA=ftla_dcsum.o ftla_pdgeqrf.o ftla_pdlarfb.o ftla_pdgetrf.o
CFTLA=ftla_ccsum.o ftla_pcgeqrf.o ftla_pclarfb.o ftla_pcgetrf.o
ZFTLA=ftla_zcsum.o ftla_pzgeqrf.o ftla_pzlarfb.o ftla_pzgetrf.o

libftla.a: libftla.a($(SFTLA) $(DFTLA) $(CFTLA) $(ZFTLA) ftla_cof.o ftla_ftwork.o util_matrix.o util_inject.o)

testing_ft%.x: $(HELPERS) $$($$*_ORIG) ft%_main.o libftla.a 
	$(FC) $(CFLAGS) -o $@ $^ $(LDFLAGS) $(BLACSLIB)  


clean:
	rm -f *.a *.o *.x

#.PRECIOUS: %.o

