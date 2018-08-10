#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mersenne.h"
#include <lapacke.h>
#include <time.h>
#include <sys/time.h>

/* Lattice size, adjust */
#define NT 128
#define NX 128
#define ND 2
#define NDIRS (2*ND)

#define VOLUME (NT*NX)

#define MAX_CHANGES 40


/* In vec_ops.c */
void vec_zero( double a[NT][NX] );
void cg_propagator( double propagator[NT][NX], double source[NT][NX] );
void fM( double chi[NT][NX], double psi[NT][NX] );

