
CC=gcc 
CFLAGS= -march=native -Wall -Wextra -std=c99 -g
LIB= -lm

DEPS=Makefile Thirring.h mersenne.h

default: Thirring Thirring_exp


.PHONY: tests
tests: tests/test_worldline
	tests/test_worldline

%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

fermionbag: fermionbag.o mersenne_inline.o vec_ops.o fermion_matrix.o full_determinant.o fluctuation_determinant.o $(DEPS)
	$(CC) $(CFLAGS) -o fermionbag fermionbag.o mersenne_inline.o vec_ops.o fermion_matrix.o full_determinant.o fluctuation_determinant.o -llapack $(LIB)

worldline: worldline.o mersenne_inline.o measurements.o $(DEPS)
	$(CC) $(CFLAGS) -o worldline worldline.o mersenne_inline.o $(LIB)

gauged: gauged.o mersenne_inline.o $(DEPS)
	$(CC) $(CFLAGS) -o gauged gauged.o mersenne_inline.o -llapack $(LIB)

dimer: dimer.o mersenne_inline.o $(DEPS)
	$(CC) $(CFLAGS) -o dimer dimer.o mersenne_inline.o -llapack $(LIB)


clean:
	rm -f *.o  
