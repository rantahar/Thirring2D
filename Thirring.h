#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mersenne.h"
#include <lapacke/lapacke.h>
#include <time.h>
#include <sys/time.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846264338327 
#endif

/* Lattice size and dimensions */
#define NT 48
#define NX 48
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

#define CG_ACCURACY 1e-16
#define CG_MAX_ITER 10000

/* Propability of exiting in the monomer moving worm update */
#define flip_exit_propability 0.2


/* Functions in vec_ops.c */
void vec_zero( double **a );
void vec_one(double **a);
void vec_add(double **a, double **b);
void vec_d_mul(double **a, double d);
void vec_zero_occupied(double **a);
void free_vector(double ** a);
double ** alloc_vector();
void cg_propagator( double **propagator, double **source );

void fM( double **chi, double **psi );
void fM_transpose(double **chi, double **psi );


int cg_MdM_occupied( double *chi, double *psi );
void fM_occupied( double *chi, double *psi );
double action(double *psi);
void vec_gaussian(double *a);
double * alloc_field();


