
CC=gcc 
CFLAGS=-march=native -Wall -Wextra -std=c99 -O3 -llapack -lm 

DEPS=Makefile Thirring.h mersenne.h

default: Thirring Thirring_exp

%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

fermionbag: fermionbag.o mersenne_inline.o vec_ops.o fermion_matrix.o full_determinant.o fluctuation_determinant.o $(DEPS)
	$(CC) -o fermionbag fermionbag.o mersenne_inline.o vec_ops.o fermion_matrix.o full_determinant.o fluctuation_determinant.o $(CFLAGS)

worldline: worldline.o mersenne_inline.o measurements.o vec_ops.o fermion_matrix.o full_determinant.o fluctuation_determinant.o $(DEPS)
	$(CC) -o worldline worldline.o mersenne_inline.o measurements.o vec_ops.o fermion_matrix.o fluctuation_determinant.o $(CFLAGS)

gauged: gauged.o mersenne_inline.o $(DEPS)
	$(CC) -o gauged gauged.o mersenne_inline.o $(CFLAGS)

dimer: dimer.o mersenne_inline.o $(DEPS)
	$(CC) -o dimer dimer.o mersenne_inline.o $(CFLAGS)

clean:
	rm -f *.o  
