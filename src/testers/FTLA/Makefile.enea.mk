#FFLAGS +=  -fsecond-underscore
#CFLAGS += -Wall 
CFLAGS += -g -O3
CFLAGS += -DINJECT
#CFLAGS += -DUSE_CoF
#CFLAGS += -DNO_CHECKPOINT
#CFLAGS += -DVERIFY_CHK
#CFLAGS += -DNO_QLOCALCOPY
#CFLAGS += -DNO_EXTRAFLOPS
#CFLAGS += -DTIMING
#CFLAGS += -DPRINTMATRIX

CFLAGS += -I${MKLROOT}/include
BLACSLIB = -L/opt/intel/compilers_and_libraries/linux/mkl/lib/intel64 -lmkl_scalapack_lp64 -lmkl_intel_lp64 -lmkl_blacs_intelmpi_lp64 -lmkl_lapack95_lp64 -lmkl_core -mkl -lm -ldl

CC = mpiicc
FC = mpiifort

all: libftla.a helpersftla.a testing_ftdtr.x testing_ftdqr.x 

.SECONDEXPANSION:
HELPERS = psmatgen.o pdmatgen.o pcmatgen.o pzmatgen.o pmatgeninc.o 
dqr_ORIG = pdgeqrrv.o #pdgetrrv.o #pdgeqrf.o pdlarfb.o pdgetrf.o 
dtr_ORIG = pdgetrrv.o pdlaprnt.o pdgeqrf.o pdlarfb.o pdgetrf.o pdgetf2.o

SFTLA=ftla_scsum.o ftla_psgeqrf.o ftla_pslarfb.o ftla_psgetrf.o
DFTLA=ftla_dcsum.o ftla_pdgeqrf.o ftla_pdlarfb.o ftla_pdgetrf.o
CFTLA=ftla_ccsum.o ftla_pcgeqrf.o ftla_pclarfb.o ftla_pcgetrf.o
ZFTLA=ftla_zcsum.o ftla_pzgeqrf.o ftla_pzlarfb.o ftla_pzgetrf.o

libftla.a: libftla.a($(SFTLA) $(DFTLA) $(CFTLA) $(ZFTLA) ftla_cof.o ftla_ftwork.o util_matrix.o util_inject.o )

helpersftla.a: helpersftla.a( $(HELPERS) $(dqr_ORIG) $(dtr_ORIG) )

testing_ft%.x: $$($$*_ORIG) ft%_main.o libftla.a helpersftla.a
	$(FC) $(CFLAGS) -nofor-main -o $@ $^  $(BLACSLIB)


clean:
	rm -f *.a *.o *.x

#.PRECIOUS: %.o
