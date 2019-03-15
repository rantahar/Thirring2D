
CC=gcc 
CFLAGS= -march=native -Wall -Wextra -std=c99 -O3
LIB= -lm

DEPS=Makefile Thirring.h mersenne.h

default: Thirring Thirring_exp

%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

fermionbag: fermionbag.o mersenne_inline.o vec_ops.o fermion_matrix.o full_determinant.o fluctuation_determinant.o $(DEPS)
	$(CC) $(CFLAGS) -o fermionbag fermionbag.o mersenne_inline.o vec_ops.o fermion_matrix.o full_determinant.o fluctuation_determinant.o -llapack $(LIB)

worldline: worldline.o mersenne_inline.o measurements.o vec_ops.o fermion_matrix.o full_determinant.o fluctuation_determinant.o $(DEPS)
	$(CC) $(CFLAGS) -o worldline worldline.o mersenne_inline.o measurements.o vec_ops.o fermion_matrix.o fluctuation_determinant.o -llapack $(LIB)

gauged: gauged.o mersenne_inline.o $(DEPS)
	$(CC) $(CFLAGS) -o gauged gauged.o mersenne_inline.o -llapack $(LIB)

dimer: dimer.o mersenne_inline.o $(DEPS)
	$(CC) $(CFLAGS) -o dimer dimer.o mersenne_inline.o -llapack $(LIB)

worldline_llr.o: worldline.c $(DEPS)
	$(CC) $(CFLAGS) -DLLR -c -o $@ $< $(CFLAGS)

worldline_llr: worldline_llr.o mersenne_inline.o  $(DEPS)
	$(CC) $(CFLAGS) -DLLR -o worldline_llr worldline_llr.o mersenne_inline.o -llapack $(LIB)


clean:
	rm -f *.o  
