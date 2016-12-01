#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mersenne.h"
#include <lapacke/lapacke.h>
#include <time.h>
#include <sys/time.h>

/* Lattice size and dimensions */
#define NT 24
#define NX 24
#define ND 2
#define NDIRS (2*ND)

#define TUP 0
#define XUP 1
#define TDN 2
#define XDN 3

#define VOLUME (NT*NX)


/* Enumerate possible values for a field */
#define MONOMER 1
#define LINK_TUP (2+TUP)   // Links enumerated as 2+direction 
#define LINK_XUP (2+XUP)
#define LINK_TDN (2+TDN)
#define LINK_XDN (2+XDN)
#define SOURCE_MONOMER (2+NDIRS)


/* Maximum number of fluctutaions from the background configuration */
#define MAX_CHANGES 20

/* Propability of exiting in the monomer moving worm update */
#define flip_exit_propability 0.2


/* Functions in vec_ops.c */
void vec_zero( double a[NT][NX] );
void cg_propagator( double propagator[NT][NX], double source[NT][NX] );
void fM( double chi[NT][NX], double psi[NT][NX] );

