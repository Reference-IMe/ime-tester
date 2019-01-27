CC = gcc
DEBUG = -g -O3
CFLAGS = -Wall $(NO_WARN_UNUSED) $(DEBUG)

all: SLGEWOS-compare

SLGEWOS-compare : SLGEWOS-compare.c
	$(CC) $(CFLAGS) SLGEWOS-compare.c -o $(BIN_DIR)/SLGEWOS-compare -llapack -lblas

clean:
	\rm *.o *.a *~ *.so
