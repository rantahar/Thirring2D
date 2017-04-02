
CC=gcc 
CFLAGS=-march=native -Wall -Wextra -std=c99 -O3 -llapack -lm #-DDEBUG

DEPS=Makefile Thirring.h

default: Thirring Thirring_exp

%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

Thirring: Thirring.o mersenne_inline.o vec_ops.o fermion_matrix.o full_determinant.o fluctuation_determinant.o $(DEPS)
	$(CC) -o Thirring Thirring.o mersenne_inline.o vec_ops.o fermion_matrix.o full_determinant.o fluctuation_determinant.o $(CFLAGS)

Thirring_exp: Thirring_exp.o mersenne_inline.o vec_ops.o fermion_matrix.o full_determinant.o fluctuation_determinant.o $(DEPS)
	$(CC) -o Thirring_exp Thirring_exp.o mersenne_inline.o vec_ops.o fermion_matrix.o full_determinant.o fluctuation_determinant.o $(CFLAGS)

clean:
	rm -f *.o  
