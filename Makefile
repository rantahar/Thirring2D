
CC=gcc 
CFLAGS=-march=native -Wall -Wextra -std=c99 -O3 -llapack -lm

DEPS=Makefile Thirring.h

%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

Thirring: Thirring.o mersenne_inline.o vec_ops.o $(DEPS)
	$(CC) -o Thirring Thirring.o mersenne_inline.o vec_ops.o $(CFLAGS)

clean:
	rm -f *.o  
