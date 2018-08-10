
CC=gcc 
CFLAGS=-march=native -Wall -Wextra -std=c99 -O3 -llapack -lm #-DDEBUG

DEPS=Makefile Thirring.h

default: Thirring Thirring_exp

%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

Thirring: hmc.o mersenne_inline.o vec_ops.o fermion_matrix.o full_determinant.o fluctuation_determinant.o $(DEPS)
	$(CC) -o hmc hmc.o mersenne_inline.o vec_ops.o fermion_matrix.o full_determinant.o fluctuation_determinant.o $(CFLAGS)

worldline: worldline.o mersenne_inline.o measurements.o vec_ops.o fermion_matrix.o full_determinant.o fluctuation_determinant.o $(DEPS)
	$(CC) -o worldline worldline.o mersenne_inline.o measurements.o vec_ops.o fermion_matrix.o fluctuation_determinant.o $(CFLAGS)

gauged: 

dimer:

clean:
	rm -f *.o  
