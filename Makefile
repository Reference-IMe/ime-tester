#ROOT_DIR := $(dir $(realpath $(lastword $(MAKEFILE_LIST))))
#$(info $(ROOT_DIR))

MACHINETYPE = cineca|enea|ubuntu
MPIFLAVOUR = intel|ch|open|spectrum
LIBTYPE  = mkl|src|sys
ENERGYREADING = cresco6

PROJECT_DIR = $(CURDIR)
BIN_DIR = $(PROJECT_DIR)/bin
SRC_DIR = $(PROJECT_DIR)/src
TST_DIR = $(SRC_DIR)/testers

IME_VERSION       = latest
IME_REPO          = https://gitlab.com/marcello.artioli/ime-lib.git
IME_TAG           = $(IME_VERSION)
IME_LIB_DIR       = $(TST_DIR)/IMe/$(IME_TAG)

LAPACK_VERSION    = 3.9.0
LAPACK_TAG        = lapack-$(LAPACK_VERSION)
LAPACK_REPO       = https://github.com/Reference-LAPACK/lapack-pre-github-historical-releases.git
LAPACK_LIB_DIR    = $(TST_DIR)/LAPACK/$(LAPACK_TAG)

SCALAPACK_VERSION = 2.1.0
SCALAPACK_TAG     = v$(SCALAPACK_VERSION)
SCALAPACK_REPO    = https://github.com/Reference-ScaLAPACK/scalapack.git
SCALAPACK_LIB_DIR = $(TST_DIR)/ScaLAPACK/scalapack-$(SCALAPACK_VERSION)

FTLA_VERSION      = rSC13
FTLA_TAG          = ftla-$(FTLA_VERSION)
FTLA_ARCHIVE      = $(FTLA_TAG).tgz
FTLA_REPO         = https://icl.utk.edu/projectsfiles/ft-la/software/$(FTLA_ARCHIVE)
FTLA_LIB_DIR      = $(TST_DIR)/FTLA/$(FTLA_TAG)

SDS_VERSION       = 2.0.0
SDS_TAG           = $(SDS_VERSION)
SDS_REPO          = https://github.com/antirez/sds.git
SDS_LIB_DIR       = $(SRC_DIR)/helpers/simple_dynamic_strings

OPTIMIZATION = -O3
DEBUG = -g
NO_WARN_UNUSED = -Wno-unused-but-set-variable -Wno-unused-variable

CFLAGS_cineca_intel = -lifcore -w3 -wd1418 -wd2259
CFLAGS_cineca_open  = -lgfortran
CFLAGS_cineca_spectrum = -I$(MPI_ROOT)/include

CFLAGS_enea_intel   = -lifcore -w3 -wd1418 -wd2259
CFLAGS_enea_open    = -lgfortran #-Wno-unknown-pragmas

CFLAGS_ubuntu_open    = -lgfortran #-Wno-unknown-pragmas


CFLAGS = $(OPTIMIZATION) $(DEBUG) -DINJECT -Wall $(CFLAGS_$(machine)_$(mpi))
ifeq ($(energy),cresco6)
CFLAGS += -DCRESCO_POWERCAP
endif

FFLAGS_cineca_intel = -nofor-main
FFLAGS_cineca_open  = 
FFLAGS_cineca_spectrum  = -I$(MPI_ROOT)/include -qextname

FFLAGS_enea_intel = 
FFLAGS_enea_open = 

FFLAGS_ubuntu_open = -fallow-argument-mismatch

FFLAGS = $(OPTIMIZATION) $(DEBUG) $(FFLAGS_$(machine)_$(mpi))

## check machine/platform, library type and mpi flavour
ifeq ($(machine),)
  $(error Unspecified machine, please specify 'machine = $(MACHINETYPE)')
else
  ifneq ($(machine),$(filter $(machine),$(subst |, ,$(MACHINETYPE))))
    $(error Unknown machine '$(machine)', please specify 'machine = $(MACHINETYPE)')
  endif
endif

ifeq ($(library),)
  $(error Unspecified library type, please specify 'library = $(LIBTYPE)')
else
  ifneq ($(library),$(filter $(library),$(subst |, ,$(LIBTYPE))))
    $(error Unknown machine '$(library)', please specify 'library = $(LIBTYPE)')
  endif
endif

ifeq ($(mpi),)
  $(error Unspecified mpi flavour, please specify 'mpi = $(MPIFLAVOUR)')
else
  ifneq ($(mpi),$(filter $(mpi),$(subst |, ,$(MPIFLAVOUR))))
    $(error Unknown mpi flavour '$(mpi)', please specify 'mpi = $(MPIFLAVOUR)')
  endif
endif


## set compilers, depending on MPI flavour
ifeq ($(mpi),intel)
  MPICC = mpiicc
  MPIFC = mpiifort
endif
ifeq ($(mpi),open)
  MPICC = mpicc
  MPIFC = mpif77
endif
ifeq ($(mpi),ch)
  MPICC = mpicc
  MPIFC = mpif77
endif
ifeq ($(mpi),spectrum)
  MPICC = mpicc
  MPIFC = mpif77
endif

## machine/platform flags
##
## how to make a smart switch: https://stackoverflow.com/questions/200205/good-way-to-do-a-switch-in-a-makefile

# marconi|galileo -> cineca
SEQ_MFLAGS_cineca_intel_mkl = -I$(MKL_INC) -L$(MKL_LIB) -lmkl_intel_lp64 -lmkl_core -lmkl_sequential -lm
PAR_MFLAGS_cineca_intel_mkl = -I$(MKL_INC) -L$(MKL_LIB) -lmkl_scalapack_lp64 -lmkl_intel_lp64 -lmkl_blacs_intelmpi_lp64 -lmkl_gf_lp64 -lmkl_core -lmkl_sequential -lm

SEQ_MFLAGS_cineca_intel_src = $(LAPACK_LIB_DIR)/liblapack.a $(LAPACK_LIB_DIR)/librefblas.a -lm
PAR_MFLAGS_cineca_intel_src = $(SCALAPACK_LIB_DIR)/libscalapack.a $(LAPACK_LIB_DIR)/liblapack.a $(LAPACK_LIB_DIR)/librefblas.a

SEQ_MFLAGS_cineca_open_src = $(LAPACK_LIB_DIR)/liblapack.a $(LAPACK_LIB_DIR)/librefblas.a -lm -ldl
PAR_MFLAGS_cineca_open_src = $(SCALAPACK_LIB_DIR)/libscalapack.a $(LAPACK_LIB_DIR)/liblapack.a $(LAPACK_LIB_DIR)/librefblas.a -lmpi_mpifh -lmpi -lm -ldl

SEQ_MFLAGS_cineca_spectrum_src = $(LAPACK_LIB_DIR)/liblapack.a $(LAPACK_LIB_DIR)/librefblas.a -lm -ldl
PAR_MFLAGS_cineca_spectrum_src = $(SCALAPACK_LIB_DIR)/libscalapack.a $(LAPACK_LIB_DIR)/liblapack.a $(LAPACK_LIB_DIR)/librefblas.a -L/m100/prod/opt/compilers/xl/16.1.1/binary/xlf/16.1.1/lib -lxlf90_r -lxlfmath -lmpi_ibm_mpifh -lm -ldl

PAR_MFLAGS_cineca_ch_mkl    = -I$(MKL_INC) -L$(MKL_LIB) -lmkl_scalapack_lp64 -lmkl_intel_lp64 -lmkl_blacs_intelmpi_lp64 -lmkl_core -lmkl_sequential -lm
PAR_MFLAGS_cineca_open_sys  = -L$(SCALAPACK_LIB) -lscalapack -L$(LAPACK_LIB) -llapack -L$(BLAS_LIB) -lblas -lifcore -lm

FTLAMAKEFILE_cineca_intel = Makefile.cineca.intel.mk
FTLAMAKEFILE_cineca_open  = Makefile.cineca.open.mk
FTLAMAKEFILE_cineca_spectrum  = Makefile.cineca.spectrum.mk

# cresco6/cresco4 -> enea
SEQ_MFLAGS_enea_intel_mkl = -I$(MKLROOT)/include -L$(MKLROOT)/lib -mkl -ldl -lm
PAR_MFLAGS_enea_intel_mkl = -I$(MKLROOT)/include -L$(MKLROOT)/lib -mkl -lmkl_scalapack_lp64 -lmkl_blacs_intelmpi_lp64 -ldl -lm

SEQ_MFLAGS_enea_intel_src = $(LAPACK_LIB_DIR)/liblapack.a $(LAPACK_LIB_DIR)/librefblas.a -lm
PAR_MFLAGS_enea_intel_src = $(SCALAPACK_LIB_DIR)/libscalapack.a $(LAPACK_LIB_DIR)/liblapack.a $(LAPACK_LIB_DIR)/librefblas.a
FTLAMAKEFILE_enea_intel = Makefile.enea.intel.mk

SEQ_MFLAGS_enea_open_src = $(LAPACK_LIB_DIR)/liblapack.a $(LAPACK_LIB_DIR)/librefblas.a -ldl -lm
PAR_MFLAGS_enea_open_src = $(SCALAPACK_LIB_DIR)/libscalapack.a $(LAPACK_LIB_DIR)/liblapack.a $(LAPACK_LIB_DIR)/librefblas.a -lmpi_mpifh -lmpi -lm -ldl
FTLAMAKEFILE_enea_open = Makefile.enea.open.mk

# ubuntu
SEQ_MFLAGS_ubuntu_open_src = $(LAPACK_LIB_DIR)/liblapack.a $(LAPACK_LIB_DIR)/librefblas.a -ldl -lm
PAR_MFLAGS_ubuntu_open_src = $(SCALAPACK_LIB_DIR)/libscalapack.a $(LAPACK_LIB_DIR)/liblapack.a $(LAPACK_LIB_DIR)/librefblas.a -lmpi_mpifh -lmpi -lgfortran -lm -ldl
FTLAMAKEFILE_ubuntu_open = Makefile.ubuntu.open.mk

PAR_MACHINEFLAGS = $(PAR_MFLAGS_$(machine)_$(mpi)_$(library))
SEQ_MACHINEFLAGS = $(SEQ_MFLAGS_$(machine)_$(mpi)_$(library))
FTLAMAKEFILE     = $(FTLAMAKEFILE_$(machine)_$(mpi))

## set targets
SEQ_EXE = # compare_solve
PAR_EXE = tester
EXE = $(addprefix $(BIN_DIR)/, $(SEQ_EXE) $(PAR_EXE) )

PAR_STD_DEP = $(TST_DIR)/test_*.h $(SRC_DIR)/helpers/*.h
SEQ_STD_DEP = $(TST_DIR)/tester_head.c $(TST_DIR)/tester_shoulder.c $(TST_DIR)/tester_tail.c $(SRC_DIR)/helpers/*.h

all: $(LAPACK_LIB_DIR)/librefblas.a \
		$(LAPACK_LIB_DIR)/liblapack.a \
		$(SCALAPACK_LIB_DIR)/libscalapack.a \
		$(FTLA_LIB_DIR)/libftla.a \
		$(SDS_LIB_DIR)/sds.o \
		$(EXE) \
		| $(BIN_DIR)

.PHONY: all clean clean_all 

$(BIN_DIR):
	mkdir -p $(BIN_DIR)

### Modules
#
# https://stackoverflow.com/questions/16315089/how-to-get-exit-status-of-a-shell-command-used-in-gnu-makefile

## IMe

.PHONY: clean_ime clone_ime remove_ime

# prereq.

$(IME_LIB_DIR):
	mkdir -p $(IME_LIB_DIR)

remove_ime:
	rm -rf $(IME_LIB_DIR)

clone_ime: $(IME_LIB_DIR)
	@if [ ! -d $(IME_LIB_DIR)/src ]; then \
		git clone --depth 1 --branch master $(IME_REPO) $(IME_LIB_DIR) ; \
	fi
#	@if [ ! -f $(IME_LIB_DIR)/README.md ]; then \
		cd $(IME_LIB_DIR) && git checkout master; \
		CLONED=$$?; \
		if [ $$CLONED -ne 0 ]; then \
			git clone --depth 1 --branch master $(IME_REPO) $(IME_LIB_DIR) ; \
		fi; \
	fi

$(IME_LIB_DIR)/src/*: clone_ime


## LAPACK

.PHONY: clean_lapack clone_lapack remove_lapack

# prereq.

$(LAPACK_LIB_DIR):
	mkdir -p $(LAPACK_LIB_DIR)

remove_lapack:
	rm -rf $(LAPACK_LIB_DIR)

clone_lapack: $(LAPACK_LIB_DIR)
	@if [ ! -f $(LAPACK_LIB_DIR)/make.inc ]; then \
		cd $(LAPACK_LIB_DIR) && git checkout $(LAPACK_TAG); \
		CLONED=$$?; \
		if [ $$CLONED -ne 0 ]; then \
			git clone --depth 1 --branch $(LAPACK_TAG) $(LAPACK_REPO) $(LAPACK_LIB_DIR) ; \
		fi; \
		cp -v $(LAPACK_LIB_DIR)/../make.inc.copy $(LAPACK_LIB_DIR)/make.inc; \
	fi
#	if [ ! -f $(LAPACK_LIB_DIR)/make.inc ]; then cp -v $(LAPACK_LIB_DIR)/../make.inc.copy $(LAPACK_LIB_DIR)/make.inc; fi

$(LAPACK_LIB_DIR)/make.inc: clone_lapack

# build
$(LAPACK_LIB_DIR)/librefblas.a: $(LAPACK_LIB_DIR)/make.inc
	cd $(LAPACK_LIB_DIR) && $(MAKE) CC=$(MPICC) FC=$(MPIFC) CFLAGS="$(CFLAGS)" FFLAGS="$(FFLAGS)" blaslib

$(LAPACK_LIB_DIR)/liblapack.a: $(LAPACK_LIB_DIR)/make.inc
#$(LAPACK_LIB_DIR)/librefblas.a
	cd $(LAPACK_LIB_DIR) && $(MAKE) CC=$(MPICC) FC=$(MPIFC) CFLAGS="$(CFLAGS)" FFLAGS="$(FFLAGS)" lapacklib

lapack: $(LAPACK_LIB_DIR)/librefblas.a $(LAPACK_LIB_DIR)/liblapack.a


## ScaLAPACK

.PHONY: clean_scalapack clone_scalapack remove_scalapack

# prereq.

$(SCALAPACK_LIB_DIR):
	mkdir -p $(SCALAPACK_LIB_DIR)

remove_scalapack:
	rm -rf $(SCALAPACK_LIB_DIR)

clone_scalapack: $(SCALAPACK_LIB_DIR)
	@if [ ! -f $(SCALAPACK_LIB_DIR)/SLmake.inc ]; then \
		cd $(SCALAPACK_LIB_DIR) && git checkout $(SCALAPACK_TAG); \
		CLONED=$$?; \
		if [ $$CLONED -ne 0 ]; then \
			git clone --depth 1 --branch $(SCALAPACK_TAG) $(SCALAPACK_REPO) $(SCALAPACK_LIB_DIR) ; \
		fi; \
		cp -v $(SCALAPACK_LIB_DIR)/../SLmake.inc.copy $(SCALAPACK_LIB_DIR)/SLmake.inc; \
	fi
#	if [ ! -f $(SCALAPACK_LIB_DIR)/SLmake.inc ]; then cp -v $(SCALAPACK_LIB_DIR)/../SLmake.inc.copy $(SCALAPACK_LIB_DIR)/SLmake.inc; fi

$(SCALAPACK_LIB_DIR)/SLmake.inc: clone_scalapack

# build
#
# do not use "-j" flag: compilation inconsistency!
# ScaLAPACK's makefile has been modified to accept a variable for pointing to the local LAPACK lib in this repository
$(SCALAPACK_LIB_DIR)/libscalapack.a: $(SCALAPACK_LIB_DIR)/SLmake.inc $(LAPACK_LIB_DIR)/librefblas.a $(LAPACK_LIB_DIR)/liblapack.a
	cd $(SCALAPACK_LIB_DIR) && $(MAKE) CC=$(MPICC) FC=$(MPIFC) CCFLAGS="$(CFLAGS)" FCFLAGS="$(FFLAGS)" LAPACK_DIR=$(LAPACK_LIB_DIR) lib

scalapack: $(SCALAPACK_LIB_DIR)/libscalapack.a

$(TST_DIR)/ScaLAPACK/%.o: $(TST_DIR)/ScaLAPACK/%.f
#	$(MPIFC) $(FFLAGS) $< -o $@ $(PAR_MACHINEFLAGS)
	$(MPIFC) $(FFLAGS) -c $< -o $@ #$(PAR_MACHINEFLAGS)

## FTLA

.PHONY: clean_ftla clone_ftla remove_ftla

# prereq.

$(FTLA_LIB_DIR):
	mkdir -p $(FTLA_LIB_DIR)

remove_ftla:
	rm -rf $(FTLA_LIB_DIR)

clone_ftla: $(FTLA_LIB_DIR)
	@if [ ! -f $(FTLA_LIB_DIR)/slp.h ]; then \
		cd $(FTLA_LIB_DIR); \
		if [ ! -f $(FTLA_ARCHIVE) ]; then \
			wget -v $(FTLA_REPO); \
		fi; \
		tar zxvf $(FTLA_ARCHIVE) --strip-components=1; \
		for file in pdgeqrf pdlarfb pdgetrf pdgetf2 ; do \
			cp -v $(SCALAPACK_LIB_DIR)/SRC/$$file.f $(FTLA_LIB_DIR)/"$$file".f; \
		done; \
		cp -v $(SCALAPACK_LIB_DIR)/TESTING/LIN/pdgetrrv.f $(FTLA_LIB_DIR)/pdgetrrv.f; \
		cp -v $(SCALAPACK_LIB_DIR)/TOOLS/pdlaprnt.f $(FTLA_LIB_DIR)/pdlaprnt.f; \
		for file in ftdqr_main.c ftdtr_main.c ftla_dcsum.c ftla_ftwork.c ; do \
			mv -v $(FTLA_LIB_DIR)/$$file $(FTLA_LIB_DIR)/"$$file".orig; \
			cp -v $(FTLA_LIB_DIR)/../"$$file".copy $(FTLA_LIB_DIR)/$$file; \
		done; \
		for mf in $$(ls $(FTLA_LIB_DIR)/../*.copy); do \
			cp -v $$mf $(FTLA_LIB_DIR)/$$(basename $$mf .copy); \
		done; \
	fi

#clone_ftla: $(FTLA_LIB_DIR)
#	@if [ ! -f $(FTLA_LIB_DIR)/slp.h ]; then \
		cd $(FTLA_LIB_DIR); \
		if [ ! -f $(FTLA_ARCHIVE) ]; then \
			 wget -v $(FTLA_REPO); \
		fi; \
		tar zxvf $(FTLA_ARCHIVE) --strip-components=1 --skip-old-files; \
		for file in pdgeqrf pdlarfb pdgetrf pdgetf2 ; do \
			if [ ! -f $(FTLA_LIB_DIR)/"$$file".f ]; then \
				cp -v $(SCALAPACK_LIB_DIR)/SRC/$$file.f $(FTLA_LIB_DIR)/"$$file".f; \
			fi; \
		done; \
		if [ ! -f $(FTLA_LIB_DIR)/pdgetrrv.f ]; then cp -v $(SCALAPACK_LIB_DIR)/TESTING/LIN/pdgetrrv.f $(FTLA_LIB_DIR)/pdgetrrv.f; fi; \
		if [ ! -f $(FTLA_LIB_DIR)/pdlaprnt.f ]; then cp -v $(SCALAPACK_LIB_DIR)/TOOLS/pdlaprnt.f $(FTLA_LIB_DIR)/pdlaprnt.f; fi; \
		for file in ftdqr_main.c ftdtr_main.c ftla_dcsum.c ftla_ftwork.c slp.h ; do \
			if [ ! -f $(FTLA_LIB_DIR)/"$$file".orig ]; then \
				mv -v $(FTLA_LIB_DIR)/$$file $(FTLA_LIB_DIR)/"$$file".orig; \
				cp -v $(FTLA_LIB_DIR)/../"$$file".copy $(FTLA_LIB_DIR)/$$file; \
			fi; \
		done; \
		for file in $$(ls $(FTLA_LIB_DIR)/../*.copy); do \
			mf=$$(basename $$file .copy); \
			if [ ! -f $(FTLA_LIB_DIR)/$$mf ]; then \
				cp -v $(FTLA_LIB_DIR)/../$$mf.copy $(FTLA_LIB_DIR)/$$mf; \
			fi; \
		done; \
	fi

$(FTLA_LIB_DIR)/slp.h: clone_ftla

# build
	
$(FTLA_LIB_DIR)/libftla.a: $(FTLA_LIB_DIR)/slp.h $(LAPACK_LIB_DIR)/librefblas.a $(LAPACK_LIB_DIR)/liblapack.a $(SCALAPACK_LIB_DIR)/libscalapack.a
	cd $(FTLA_LIB_DIR) && $(MAKE) FC=$(MPIFC) CC=$(MPICC) CFLAGS="$(CFLAGS)" FFLAGS="$(FFLAGS)" LAPACK_DIR=$(LAPACK_LIB_DIR) SCALAPACK_DIR=$(SCALAPACK_LIB_DIR) -f $(FTLAMAKEFILE)

ftla: $(FTLA_LIB_DIR)/libftla.a

## SDS

.PHONY: clean_sds clone_sds remove_sds

# prereq.

$(SDS_LIB_DIR):
	mkdir -p $(SDS_LIB_DIR)

remove_sds:
	rm -rf $(SDS_LIB_DIR)

clone_sds: $(SCALAPACK_LIB_DIR)
	@if [ ! -f $(SDS_LIB_DIR)/sds.h ]; then \
		cd $(SDS_LIB_DIR) && git checkout $(SDS_TAG); \
		CLONED=$$?; \
		if [ $$CLONED -ne 0 ]; then \
			git clone --depth 1 --branch $(SDS_TAG) $(SDS_REPO) $(SDS_LIB_DIR) ; \
		fi; \
	fi

$(SDS_LIB_DIR)/sds.c $(SDS_LIB_DIR)/sds.h: clone_sds

# build

$(SDS_LIB_DIR)/sds.o: $(SDS_LIB_DIR)/sds.c $(SDS_LIB_DIR)/sds.h
	cd $(SDS_LIB_DIR) && $(CC) $(CFLAGS) -std=c99 -pedantic -c sds.c

sds: $(SDS_LIB_DIR)/sds.o


# static linking experiment
#
#$(BIN_DIR)/run_SPK-SV_mkl_stat: $(TST_DIR)/run_SPK-SV.c \
				$(TST_DIR)/test_ScaLAPACK_pDGESV.h \
				$(TST_DIR)/ScaLAPACK/ScaLAPACK_pDGESV.h \
				| $(BIN_DIR)
#	$(MPICC) $(CFLAGS)  -DMKL_ILP64 -I${MKLROOT}/include -lifcore -o $(BIN_DIR)/run_SPK-SV_mkl_stat $(TST_DIR)/run_SPK-SV.c $(PAR_MFLAGS_$(machine)_$(mpi)_mkl) ${MKLROOT}/lib/intel64/libmkl_scalapack_ilp64.a -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_ilp64.a ${MKLROOT}/lib/intel64/libmkl_sequential.a ${MKLROOT}/lib/intel64/libmkl_core.a ${MKLROOT}/lib/intel64/libmkl_blacs_intelmpi_ilp64.a -Wl,--end-group -lpthread -lm -ldl

#$(BIN_DIR)/run_SPK-SV_mkl_stat: $(TST_DIR)/run_SPK-SV.c \
				$(TST_DIR)/test_ScaLAPACK_pDGESV.h \
				$(TST_DIR)/ScaLAPACK/ScaLAPACK_pDGESV.h \
				| $(BIN_DIR)
#	$(MPICC) $(CFLAGS) -lifcore -o $(BIN_DIR)/run_SPK-SV_mkl_stat $(TST_DIR)/run_SPK-SV.c -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_scalapack_lp64.a ${MKLROOT}/lib/intel64/libmkl_blacs_intelmpi_lp64.a ${MKLROOT}/lib/intel64/libmkl_intel_lp64.a ${MKLROOT}/lib/intel64/libmkl_sequential.a ${MKLROOT}/lib/intel64/libmkl_core.a -static_mpi -Wl,--end-group -lpthread -lm
#
# => static linking: ueseless on cresco6 and galileo


$(BIN_DIR)/tester: $(TST_DIR)/tester.c \
				$(TST_DIR)/tester*.h \
				$(TST_DIR)/IMe/latest/src/* \
				$(PAR_STD_DEP) \
				$(SRC_DIR)/helpers/simple_dynamic_strings/sds.o \
				$(TST_DIR)/ScaLAPACK/*.h \
				$(TST_DIR)/ScaLAPACK/psgetrf_cp.o \
				$(TST_DIR)/ScaLAPACK/psgetrf_cpx.o \
				$(TST_DIR)/ScaLAPACK/pdgetrf_cp.o \
				$(TST_DIR)/ScaLAPACK/pdgetrf_cpx.o \
				$(TST_DIR)/ScaLAPACK/pdgeqrf_cp.o \
				$(TST_DIR)/ScaLAPACK/pdgetrf_cs.o \
				$(TST_DIR)/ScaLAPACK/pdgesv_nopivot.o \
				$(TST_DIR)/ScaLAPACK/pdgetf2_nopivot.o \
				$(TST_DIR)/ScaLAPACK/pdgetrf_nopivot.o \
				$(TST_DIR)/FTLA/*.h \
				| $(FTLA_LIB_DIR)/libftla.a \
				$(BIN_DIR)
#	$(MPICC) $(CFLAGS) $(TST_DIR)/tester.c -o $(BIN_DIR)/tester $(SRC_DIR)/helpers/simple_dynamic_strings/sds.o $(TST_DIR)/ScaLAPACK/pdgetrf_cp.o $(TST_DIR)/ScaLAPACK/pdgeqrf_cp.o $(FTLA_LIB_DIR)/libftla.a $(FTLA_LIB_DIR)/helpersftla.a $(SCALAPACK_LIB_DIR)/libscalapack.a $(LAPACK_LIB_DIR)/librefblas.a $(LAPACK_LIB_DIR)/liblapack.a -L$(TST_DIR)/ScaLAPACK $(PAR_MACHINEFLAGS)
#	$(MPICC) $(CFLAGS) $< -o $(BIN_DIR)/tester $(SRC_DIR)/helpers/simple_dynamic_strings/sds.o $(TST_DIR)/ScaLAPACK/pdgetrf_cp.o $(TST_DIR)/ScaLAPACK/pdgeqrf_cp.o $(FTLA_LIB_DIR)/libftla.a $(FTLA_LIB_DIR)/helpersftla.a $(SCALAPACK_LIB_DIR)/libscalapack.a $(LAPACK_LIB_DIR)/librefblas.a $(LAPACK_LIB_DIR)/liblapack.a -L$(TST_DIR)/ScaLAPACK $(PAR_MACHINEFLAGS)
#	$(MPICC) $(CFLAGS) $(TST_DIR)/tester.c -o $(BIN_DIR)/tester $(SRC_DIR)/helpers/simple_dynamic_strings/sds.o $(TST_DIR)/ScaLAPACK/pdgetrf_cs.o $(TST_DIR)/ScaLAPACK/pdgetrf_cp.o $(TST_DIR)/ScaLAPACK/pdgeqrf_cp.o $(FTLA_LIB_DIR)/libftla.a $(FTLA_LIB_DIR)/helpersftla.a $(SCALAPACK_LIB_DIR)/libscalapack.a $(LAPACK_LIB_DIR)/librefblas.a $(LAPACK_LIB_DIR)/liblapack.a -L$(TST_DIR)/ScaLAPACK $(PAR_MACHINEFLAGS)
	$(MPICC) $(CFLAGS) -o $(BIN_DIR)/tester \
		$(TST_DIR)/tester.c \
		$(SRC_DIR)/helpers/simple_dynamic_strings/sds.o \
		$(TST_DIR)/ScaLAPACK/psgetrf_cp.o \
		$(TST_DIR)/ScaLAPACK/psgetrf_cpx.o \
		$(TST_DIR)/ScaLAPACK/pdgetrf_cs.o \
		$(TST_DIR)/ScaLAPACK/pdgetrf_cp.o \
		$(TST_DIR)/ScaLAPACK/pdgetrf_cpx.o \
		$(TST_DIR)/ScaLAPACK/pdgeqrf_cp.o \
		$(TST_DIR)/ScaLAPACK/pdgesv_nopivot.o \
		$(TST_DIR)/ScaLAPACK/pdgetf2_nopivot.o \
		$(TST_DIR)/ScaLAPACK/pdgetrf_nopivot.o \
		$(FTLA_LIB_DIR)/libftla.a \
		$(FTLA_LIB_DIR)/helpersftla.a \
		$(SCALAPACK_LIB_DIR)/libscalapack.a \
		$(LAPACK_LIB_DIR)/librefblas.a \
		$(LAPACK_LIB_DIR)/liblapack.a \
		-L$(TST_DIR)/ScaLAPACK \
		$(PAR_MACHINEFLAGS)


# cleanup
#
clean:
	rm -f $(EXE)
	cd $(TST_DIR)/ScaLAPACK && rm -f *.o

clean_ftla:
	cd $(FTLA_LIB_DIR) && $(MAKE) clean

clean_lapack:
	cd $(LAPACK_LIB_DIR) && $(MAKE) clean

clean_scalapack:
	cd $(SCALAPACK_LIB_DIR) && $(MAKE) clean

clean_all: clean clean_ftla clean_lapack clean_scalapack
