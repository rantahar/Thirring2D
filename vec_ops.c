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
  static double expmu, expmmu;
  static int init=1;
  if(init){
    expmu = exp(mu);
    expmmu = exp(-mu);
  }
  for (int t=0; t<NT; t++) for (int x=0; x<NX; x++) {
    chi[t][x] = 0;
    if( field[t][x] == 0 ) {
      int t2=tup[t];
      if( field[t2][x] == 0 ) {
        if( t2 > t ) chi[t][x] += 0.5 * expmu * psi[t2][x] ;
        else chi[t][x] -= 0.5 * expmu * psi[t2][x] ;
      }
      t2 = tdn[t];
      if( field[t2][x] == 0 ) {
        if( t2 > t ) chi[t][x] += 0.5 * expmmu * psi[t2][x] ;
        else chi[t][x] -= 0.5 * expmmu * psi[t2][x] ;
      }
      int x2 = xup[x];
      if( field[t][x2] == 0 ){
        if( x2 > x ) chi[t][x] += 0.5 * eta[t][x][1] * psi[t][x2] ;
        else chi[t][x] -= 0.5 * eta[t][x][1] * psi[t][x2] ;
      }
      x2 = xdn[x];
      if( field[t][x2] == 0 ){
        if( x2 > x ) chi[t][x] += 0.5 * eta[t][x][1] * psi[t][x2] ;
        else chi[t][x] -= 0.5 * eta[t][x][1] * psi[t][x2] ;
      }
    }
  }
}


/* Conjugate gradient for inverting the fermion matrix
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

  if( rr_old == 0 ){
    // If the source is 0, the solution is 0
    vec_zero(propagator);
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
    if( rr < 1e-12 ) break;

    double b = rr / rr_old ;
    vec_dmul_add( p, r, p, b );

    rr_old = rr;
  }
  fM( propagator, inv );
}



























