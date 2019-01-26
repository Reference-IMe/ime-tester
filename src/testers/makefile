CC = gcc
DEBUG = -g -O3
CFLAGS = -Wall $(NO_WARN_UNUSED) $(DEBUG)

all: GJE-compare-seq

GJE-compare-seq : GJE-compare-seq.c
	$(CC) $(CFLAGS) GJE-compare-seq.c -o $(BIN_DIR)/GJE-compare-seq

clean:
	\rm *.o *.a *~ *.so
