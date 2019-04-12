CC = gcc
MPICC = mpicc
DEBUG = -g -O3
CFLAGS = -Wall $(NO_WARN_UNUSED) $(DEBUG)

all : testers

testers : SLGEWOS_compare SLGEWOPV_compare
  
SLGEWOS_compare : SLGEWOS_compare.c
	$(CC) $(CFLAGS) SLGEWOS_compare.c -o $(BIN_DIR)/SLGEWOS_compare -llapack -lblas

SLGEWOPV_compare : SLGEWOPV_compare.c
#	$(MPICC) $(CFLAGS) SLGEWOPV_compare.c -o $(BIN_DIR)/SLGEWOPV_compare /usr/lib/x86_64-linux-gnu/libmpi.so /usr/lib/x86_64-linux-gnu/libscalapack-openmpi.so.2.0 /usr/lib/x86_64-linux-gnu/libblacsCinit-openmpi.a /usr/lib/x86_64-linux-gnu/libblacs-openmpi.a
	$(MPICC) $(CFLAGS) SLGEWOPV_compare.c -o $(BIN_DIR)/SLGEWOPV_compare /usr/lib/libmpi.so /usr/lib/libscalapack-openmpi.so /usr/lib/libblacsCinit-openmpi.so /usr/lib/libblacs-openmpi.so

clean:
	\rm *.o *.a *~ *.so
