ENERGYREADING = cresco6

PROJECT_DIR = $(CURDIR)
BIN_DIR = $(PROJECT_DIR)/bin
SRC_DIR = $(PROJECT_DIR)/src
TST_DIR = $(SRC_DIR)/routines

IME_VERSION       = latest
IME_TAG           = $(IME_VERSION)
IME_REPO          = https://github.com/Reference-IMe/ime-lib.git
IME_REPO_DIR      = lib
IME_LIB_DIR       = $(TST_DIR)/IMe/$(IME_REPO_DIR)

LAPACK_VERSION    = 3.9.0
LAPACK_TAG        = lapack-$(LAPACK_VERSION)
LAPACK_REPO       = https://github.com/Reference-LAPACK/lapack-pre-github-historical-releases.git
LAPACK_REPO_DIR   = lib
LAPACK_LIB_DIR    = $(TST_DIR)/LAPACK/$(LAPACK_REPO_DIR)

SCALAPACK_VERSION = 2.1.0
SCALAPACK_TAG     = v$(SCALAPACK_VERSION)
SCALAPACK_REPO    = https://github.com/Reference-ScaLAPACK/scalapack.git
SCALAPACK_REPO_DIR = lib
SCALAPACK_LIB_DIR = $(TST_DIR)/ScaLAPACK/$(SCALAPACK_REPO_DIR)

FTLA_VERSION      = rSC13
FTLA_TAG          = ftla-$(FTLA_VERSION)
FTLA_ARCHIVE      = $(FTLA_TAG).tgz
FTLA_REPO         = https://icl.utk.edu/projectsfiles/ft-la/software/$(FTLA_ARCHIVE)
FTLA_REPO_DIR     = lib
FTLA_LIB_DIR      = $(TST_DIR)/FTLA/$(FTLA_REPO_DIR)
FTLA_MAKEFILE     = Makefile.mod

SDS_VERSION       = 2.0.0
SDS_TAG           = $(SDS_VERSION)
SDS_REPO          = https://github.com/antirez/sds.git
SDS_LIB_DIR       = $(SRC_DIR)/helpers/simple_dynamic_strings

OPTIMIZATION = -O3
DEBUG = -g
NO_WARN_UNUSED = -Wno-unused-but-set-variable -Wno-unused-variable -Wno-unknown-pragmas

CFLAGS = $(OPTIMIZATION) $(DEBUG) -DINJECT -Wall -lgfortran #$(NO_WARN_UNUSED)

ifeq ($(energy),cresco6)
CFLAGS += -DCRESCO_POWERCAP
endif

FFLAGS = $(OPTIMIZATION) $(DEBUG) -fallow-argument-mismatch

MPICC = mpicc
MPIFC = mpif77

## machine/platform flags
# ubuntu
SEQ_FLAGS   = $(LAPACK_LIB_DIR)/liblapack.a $(LAPACK_LIB_DIR)/librefblas.a -ldl -lm
PAR_FLAGS   = $(SCALAPACK_LIB_DIR)/libscalapack.a $(LAPACK_LIB_DIR)/liblapack.a $(LAPACK_LIB_DIR)/librefblas.a -lmpi_mpifh -lmpi -lgfortran -lm -ldl

## set targets
SEQ_EXE = # compare_solve
PAR_EXE = tester
EXE = $(addprefix $(BIN_DIR)/, $(SEQ_EXE) $(PAR_EXE) )

#PAR_STD_DEP = $(TST_DIR)/*/test_*.h $(SRC_DIR)/helpers/*.h
#SEQ_STD_DEP = $(TST_DIR)/tester_head.c $(TST_DIR)/tester_shoulder.c $(TST_DIR)/tester_tail.c $(SRC_DIR)/helpers/*.h

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

$(IME_LIB_DIR)/src/* $(IME_LIB_DIR)/src/helpers/*: clone_ime


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
		cp -v $(LAPACK_LIB_DIR)/../mods/make.inc.copy $(LAPACK_LIB_DIR)/make.inc; \
	fi


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
		cp -v $(SCALAPACK_LIB_DIR)/../mods/SLmake.inc.copy $(SCALAPACK_LIB_DIR)/SLmake.inc; \
	fi
#	if [ ! -f $(SCALAPACK_LIB_DIR)/SLmake.inc ]; then cp -v $(SCALAPACK_LIB_DIR)/../SLmake.inc.copy $(SCALAPACK_LIB_DIR)/SLmake.inc; fi

$(SCALAPACK_LIB_DIR)/SLmake.inc: clone_scalapack

# build
#
# do not use "-j" flag: compilation inconsistency!
# ScaLAPACK's makefile has been modified to accept a variable for pointing to the local LAPACK lib in this repository
$(SCALAPACK_LIB_DIR)/libscalapack.a: $(SCALAPACK_LIB_DIR)/SLmake.inc $(LAPACK_LIB_DIR)/librefblas.a $(LAPACK_LIB_DIR)/liblapack.a
	cd $(SCALAPACK_LIB_DIR) && $(MAKE) CC=$(MPICC) FC=$(MPIFC) CCFLAGS="$(CFLAGS)" FCFLAGS="$(FFLAGS)" LAPACK_DIR=$(LAPACK_LIB_DIR) lib

scalapack: $(SCALAPACK_LIB_DIR)/libscalapack.a $(TST_DIR)/ScaLAPACK/*.o

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
			wget -v --no-check-certificate $(FTLA_REPO); \
		fi; \
		tar zxvf $(FTLA_ARCHIVE) --strip-components=1; \
		for file in pdgeqrf pdlarfb pdgetrf pdgetf2 ; do \
			cp -v $(SCALAPACK_LIB_DIR)/SRC/$$file.f $(FTLA_LIB_DIR)/"$$file".f; \
		done; \
		cp -v $(SCALAPACK_LIB_DIR)/TESTING/LIN/pdgetrrv.f $(FTLA_LIB_DIR)/pdgetrrv.f; \
		cp -v $(SCALAPACK_LIB_DIR)/TOOLS/pdlaprnt.f $(FTLA_LIB_DIR)/pdlaprnt.f; \
		for file in ftdqr_main.c ftdtr_main.c ftla_dcsum.c ftla_ftwork.c ; do \
			mv -v $(FTLA_LIB_DIR)/$$file $(FTLA_LIB_DIR)/"$$file".orig; \
			cp -v $(FTLA_LIB_DIR)/../mods/"$$file".copy $(FTLA_LIB_DIR)/$$file; \
		done; \
		for mf in $$(ls $(FTLA_LIB_DIR)/../mods/*.copy); do \
			cp -v $$mf $(FTLA_LIB_DIR)/$$(basename $$mf .copy); \
		done; \
	fi

$(FTLA_LIB_DIR)/slp.h: clone_ftla

# build
	
$(FTLA_LIB_DIR)/libftla.a: $(FTLA_LIB_DIR)/slp.h $(LAPACK_LIB_DIR)/librefblas.a $(LAPACK_LIB_DIR)/liblapack.a $(SCALAPACK_LIB_DIR)/libscalapack.a
	cd $(FTLA_LIB_DIR) && $(MAKE) FC=$(MPIFC) CC=$(MPICC) CFLAGS="$(CFLAGS)" FFLAGS="$(FFLAGS)" LAPACK_DIR=$(LAPACK_LIB_DIR) SCALAPACK_DIR=$(SCALAPACK_LIB_DIR) -f $(FTLA_MAKEFILE) all

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


$(BIN_DIR)/tester: $(SRC_DIR)/tester.c $(SRC_DIR)/tester_*.h \
				$(SRC_DIR)/helpers/*.h \
				$(SRC_DIR)/helpers/simple_dynamic_strings/sds.o \
				$(IME_LIB_DIR)/src/* \
				$(IME_LIB_DIR)/src/helpers/* \
				$(LAPACK_LIB_DIR)/librefblas.a \
				$(LAPACK_LIB_DIR)/liblapack.a \
				$(FTLA_LIB_DIR)/libftla.a \
				$(FTLA_LIB_DIR)/helpersftla.a \
				$(SCALAPACK_LIB_DIR)/libscalapack.a \
				$(TST_DIR)/dummy/*.h \
				$(TST_DIR)/FTLA/*.h \
				$(TST_DIR)/IMe/*.h \
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
				$(BIN_DIR)

	$(MPICC) $(CFLAGS) -o $(BIN_DIR)/tester \
		$(SRC_DIR)/tester.c \
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
		$(PAR_FLAGS)


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
