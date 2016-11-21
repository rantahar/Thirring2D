#include "Thirring.h"

#define CG_MAX_ITER 10000

double chi[NT][NX];   //2 for complex number

extern int    eta[NT][NX][ND];   //Staggered eta matrix
extern double propagator[NT][NX][NT][NX];  //Storage for Vol to Vol propagator
extern double m;
extern double U;
extern double mu;

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
    
    float expmu = exp(mu);
    float expmmu = exp(-mu);
    chi[t][x] = 0.5 * expmu * eta[t][x][0] * psi[tup[t]][x] ;
    chi[t][x] = 0.5 * expmu * eta[t][x][0] * psi[tup[t]][x] ;
    chi[t][x] -= 0.5 * expmmu * eta[t][x][0] * psi[tdn[t]][x] ;
    chi[t][x] -= 0.5 * expmmu * eta[t][x][0] * psi[tdn[t]][x] ;

    chi[t][x] += 0.5 * eta[t][x][1] * psi[t][xup[x]] ;
    chi[t][x] += 0.5 * eta[t][x][1] * psi[t][xup[x]] ;
    chi[t][x] -= 0.5 * eta[t][x][1] * psi[t][xdn[x]] ;
    chi[t][x] -= 0.5 * eta[t][x][1] * psi[t][xdn[x]] ;

    chi[t][x] += m * psi[t][x] ;
    chi[t][x] += m * psi[t][x] ;

  }
}


/* Conjugate gradient version
 */
void calc_propagator_cg(double propagator[NT][NX][NT][NX] )
{
  double psi[NT][NX];
  double r[NT][NX];
  double p[NT][NX];
  double Mp[NT][NX];
  for (int t=0; t<NT; t++) for (int x=0; x<NX; x++){
    vec_zero( psi );
    psi[t][x] = 1;
    psi[tup[t]][x] = 1;

    // Invert the matrix at t,x using conjugate gradient
    //vec_assign( r, psi );
    fM( r, psi );
    vec_neg(r);
    vec_assign( p, r );
    double rr_old = vec_rdot( r, r );
    printf("  %g %d\n", rr_old, 0 );

    int k;
    for( k = 1; k < CG_MAX_ITER; k++ )
    {
      fM( Mp, p );
      double pMp = vec_rdot( p, Mp ) ;
      double a = rr_old / pMp ;
      vec_dmul_add( propagator[t][x], propagator[t][x], p, a );
      vec_dmul_add( r, r, Mp, -a );

      double rr= vec_rdot( r, r );
      printf("  %g %g %d\n", rr, pMp, k );
      if(k%1==0) sleep(1);
      if( rr < 1e-10 ) break;

      double b = rr / rr_old ;
      vec_dmul_add( p, r, p, b );

      rr_old = rr;
    }
    printf(" ( %d, %d ) %g %d\n", t,x, rr_old, k );
  }
}



























