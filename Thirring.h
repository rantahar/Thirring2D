#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mersenne.h"
#include <time.h>
#include <sys/time.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846264338327 
#endif

/* Lattice size and dimensions */
#define NT 128
#define NX 128

#define ND 2
#define NDIRS (2*ND)

#define TUP 0
#define XUP 1
#define TDN 2
#define XDN 3

#define VOLUME (NT*NX)

#define ANTISYMMETRIC //Antisymmetric boundaries
//#define SYMMETRIC   //implemented in Thirring_hop
//#define OPENX       //open in space, (anti)symmetric in time


/* Enumerate possible values for a field */
#define MONOMER 1
#define LINK_TUP (2+TUP)   // Links enumerated as 2+direction 
#define LINK_XUP (2+XUP)
#define LINK_TDN (2+TDN)
#define LINK_XDN (2+XDN)
#define SOURCE_MONOMER (2+NDIRS)
#define EMPTY -100   //A meta value for sites that don't exist

#define CG_ACCURACY 1e-30
#define CG_MAX_ITER 10000

/* Propability of exiting in the monomer moving worm update */
#define flip_exit_propability 0.2

#define FLUCTUATION_MATRIX
#define WITH_MASS_MONOMERS
#define PROPAGATOR_MATRIX

#define MAX_SECTOR 301

/* Thermalise without accept/reject */
void thermalise( int nsteps );


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

/* In fermion_matrix.c */
double fM_index( int t1, int x1, int t2, int x2 );
void calc_Dinv( );


int cg_MdM_occupied( double *chi, double *psi );
void fM_occupied( double *chi, double *psi );
double action(double *psi);
void vec_gaussian(double *a);
double * alloc_field();


/* In measurements.c */
void measure_propagator( int sign );
void measure_susceptibility();



