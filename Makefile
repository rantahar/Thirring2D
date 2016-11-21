
CC=gcc 
CFLAGS=-march=native -Wall -Wextra -std=c99 -O3 -fgcse-after-reload -ffast-math -fassociative-math -freciprocal-math -fno-signed-zeros -funroll-loops -funswitch-loops -llapack -lm

DEPS=Makefile Thirring.h

%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

Thirring: Thirring.o mersenne_inline.o $(DEPS)
	$(CC) -o Thirring Thirring.o mersenne_inline.o $(CFLAGS)

clean:
	rm -f *.o  
