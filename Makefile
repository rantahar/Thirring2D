
CC=gcc 
CFLAGS= -march=native -Wall -Wextra -std=c99 -O3
LIB= -lm

DEPS=Makefile Thirring.h mersenne.h

default: Thirring Thirring_exp


partitioned: partitioned.o mersenne_inline.o fermion_matrix.o $(DEPS)
	$(CC) $(CFLAGS) -o partitioned partitioned.o mersenne_inline.o fermion_matrix.o -llapack $(LIB)


tests/test_worldline: mersenne_inline.o tests/test_worldline.c worldline.c
	$(CC) $(CFLAGS) -DTESTING tests/test_worldline.c worldline.c mersenne_inline.o -o tests/test_worldline -lcmocka -lm

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

worldline_llr.o: worldline.c $(DEPS)
	$(CC) $(CFLAGS) -DLLR -c -o $@ $< $(CFLAGS)

worldline_llr: worldline_llr.o mersenne_inline.o  $(DEPS)
	$(CC) $(CFLAGS) -DLLR -o worldline_llr worldline_llr.o mersenne_inline.o $(LIB)

worldline_WL.o: worldline.c $(DEPS)
	$(CC) $(CFLAGS) -DWANGLANDAU -c -o $@ $< $(CFLAGS)

worldline_WL: worldline_WL.o mersenne_inline.o  $(DEPS)
	$(CC) $(CFLAGS) -DWANGLANDAU -o worldline_WL worldline_WL.o mersenne_inline.o $(LIB)

worldline_sector.o: worldline.c $(DEPS)
	$(CC) $(CFLAGS) -DMEASURE_SECTOR -c -o $@ $< $(CFLAGS)

worldline_sector: worldline_sector.o mersenne_inline.o  $(DEPS)
	$(CC) $(CFLAGS) -DMEASURE_SECTOR -o worldline_sector worldline_sector.o mersenne_inline.o $(LIB)



clean:
	rm -f *.o  
