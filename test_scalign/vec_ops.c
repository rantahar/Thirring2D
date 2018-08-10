#include "Thirring.h"

#define CG_MAX_ITER 10000

double chi[NT][NX];   //2 for complex number

extern int    eta[NT][NX][ND];   //Staggered eta matrix
extern double propagator[NT][NX][NT][NX];  //Storage for Vol to Vol propagator
extern double m;
extern double U;
extern double mu;

extern int field[NT][NX];

/* Neighbour index arrays, to be filled at the beginning
 */
extern int tup[NT],xup[NX],tdn[NT],xdn[NX];

/* Library functions
 */
void vec_neg(double a[NT][NX])
{
  for (int t=0; t<NT; t++) for (int x=0; x<NX; x++) 
      a[t][x] = -a[t][x];
}
void vec_zero(double a[NT][NX])
{
  for (int t=0; t<NT; t++) for (int x=0; x<NX; x++)
      a[t][x] = 0;
}
void vec_assign(double a[NT][NX], double b[NT][NX])
{
  for (int t=0; t<NT; t++) for (int x=0; x<NX; x++) 
      a[t][x] = b[t][x];
}
void vec_dmul_add(double a[NT][NX], double b[NT][NX], double d[NT][NX], double e)
{
  for (int t=0; t<NT; t++) for (int x=0; x<NX; x++)  
      a[t][x] = b[t][x] + e * d[t][x] ;
}
double vec_rdot(double a[NT][NX], double b[NT][NX])
{
  double s = 0 ;
  for (int t=0; t<NT; t++) for (int x=0; x<NX; x++)  
      s += a[t][x]*b[t][x];
  return( s );
}
void vec_print_lat(double a[NT][NX])
{
  for (int t=0; t<NT; t++) { 
    for (int x=0; x<NX; x++)
      printf( " ( %f , %f )", a[t][x], a[t][x] );
    printf(" \n");
  }
  printf(" \n");
}


/* The fermion matrix
 */
void fM(double chi[NT][NX], double psi[NT][NX] )
{
  for (int t=0; t<NT; t++) for (int x=0; x<NX; x++) {
    chi[t][x] = 0;
    if( field[t][x] == 0 ) {

      float expmu = exp(mu);
      float expmmu = exp(-mu);
      if( field[tup[t]][x] ==0 && tup[t] > t )
        chi[t][x] += 0.5 * expmu * eta[t][x][0] * psi[tup[t]][x] ;
      if( field[tup[t]][x] ==0 && tup[t] < t )
        chi[t][x] -= 0.5 * expmu * eta[t][x][0] * psi[tup[t]][x] ;
      if( field[tdn[t]][x] ==0 && tdn[t] > t )
        chi[t][x] += 0.5 * expmmu * eta[t][x][0] * psi[tdn[t]][x] ;
      if( field[tdn[t]][x] ==0 && tdn[t] < t )
        chi[t][x] -= 0.5 * expmmu * eta[t][x][0] * psi[tdn[t]][x] ;

      if( field[t][xup[x]] ==0 && xup[x] > x )
        chi[t][x] += 0.5 * eta[t][x][1] * psi[t][xup[x]] ;
      if( field[t][xup[x]] ==0 && xup[x] < x )
        chi[t][x] -= 0.5 * eta[t][x][1] * psi[t][xup[x]] ;
      if( field[t][xdn[x]] ==0 && xdn[x] > x )
        chi[t][x] += 0.5 * eta[t][x][1] * psi[t][xdn[x]] ;
      if( field[t][xdn[x]] ==0 && xdn[x] < x )
        chi[t][x] -= 0.5 * eta[t][x][1] * psi[t][xdn[x]] ;
    }
  }
}


/* Conjugate gradient version
 */
void cg_propagator( double propagator[NT][NX], double source[NT][NX] )
{
  double rr;
  double r[NT][NX];
  double p[NT][NX];
  double Mp[NT][NX];
  double MMp[NT][NX];
  double inv[NT][NX];
  vec_zero( inv );

  vec_assign(r,source);
  vec_assign( p, r );
  double rr_old = vec_rdot( r, r );
  //printf("  %g %d\n", rr_old, 0 );

  if( rr_old == 0 ){
    // No free sites, propagator is 0
    vec_zero(propagator);
    //printf("CG done, %d iterations, r=%g\n", 0, rr_old );
    return ;
  }

  int k;
  for( k = 1; k < CG_MAX_ITER; k++ )
  {
    fM( Mp, p );
    fM( MMp, Mp );
    double pMp = vec_rdot( p, MMp ) ;
    double a = rr_old / pMp ;
    vec_dmul_add( inv, inv, p, a );
    vec_dmul_add( r, r, MMp, -a );

    rr = vec_rdot( r, r );
    //printf("  %g %g %d\n", rr, pMp, k );
    //if(k%1==0) usleep(100);
    if( rr < 1e-12 ) break;

    double b = rr / rr_old ;
    vec_dmul_add( p, r, p, b );

    rr_old = rr;
  }
  fM( propagator, inv );
  //printf("CG done, %d iterations, r=%g\n", k, rr );
}



























