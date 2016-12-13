#include "Thirring.h"

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
void vec_one(double a[NT][NX])
{
  for (int t=0; t<NT; t++) for (int x=0; x<NX; x++)
      a[t][x] = 1;
}
void vec_set(double a[NT][NX], double d)
{
  for (int t=0; t<NT; t++) for (int x=0; x<NX; x++)
      a[t][x] = d;
}
void vec_d_mul(double a[NT][NX], double d)
{
  for (int t=0; t<NT; t++) for (int x=0; x<NX; x++)
      a[t][x] = a[t][x]*d;
}
void vec_assign(double a[NT][NX], double b[NT][NX])
{
  for (int t=0; t<NT; t++) for (int x=0; x<NX; x++) 
      a[t][x] = b[t][x];
}
void vec_add(double a[NT][NX], double b[NT][NX])
{
  for (int t=0; t<NT; t++) for (int x=0; x<NX; x++)  
      a[t][x] += b[t][x];
}
void vec_dmul_add(double a[NT][NX], double b[NT][NX], double d[NT][NX], double e)
{
  for (int t=0; t<NT; t++) for (int x=0; x<NX; x++)  
      a[t][x] = b[t][x] + e * d[t][x] ;
}
inline double vec_dot(double a[NT][NX], double b[NT][NX])
{
  double s = 0 ;
  for (int t=0; t<NT; t++) for (int x=0; x<NX; x++)
      s += a[t][x]*b[t][x];
  return( s );
}

void vec_zero_occupied(double a[NT][NX]) {
  for (int t=0; t<NT; t++) for (int x=0; x<NX; x++)
   if( field[t][x] > 0 ) a[t][x] = 0; 
}

void vec_gaussian(double a[NT][NX])
{
  for (int t=0; t<NT; t++) for (int x=0; x<NX; x++){ 
    double x1,x2,r=2;
    x1 = mersenne(), x2=mersenne();
    a[t][x] = sqrt( -2*log(x1) )*cos(2*M_PI*x2);
    if(x<NX-1){
      a[t][++x] = sqrt( -2*log(x1) )*sin(2*M_PI*x2);
    }
  }
}

void vec_print_lat(double a[NT][NX])
{
  for (int t=0; t<NT; t++) { 
    for (int x=0; x<NX; x++)
      printf( " %8.2f ", a[t][x] );
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
    } else {
      chi[t][x] = psi[t][x]; //An occupied site, matrix becomes identity
    }
  }
}

void fM_transpose(double chi[NT][NX], double psi[NT][NX] )
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
        if( t2 > t ) chi[t][x] -= 0.5 * expmmu * psi[t2][x] ;
        else chi[t][x] += 0.5 * expmmu * psi[t2][x] ;
      }
      t2 = tdn[t];
      if( field[t2][x] == 0 ) {
        if( t2 > t ) chi[t][x] -= 0.5 * expmu * psi[t2][x] ;
        else chi[t][x] += 0.5 * expmu * psi[t2][x] ;
      }
      int x2 = xup[x];
      if( field[t][x2] == 0 ){
        if( x2 > x ) chi[t][x] -= 0.5 * eta[t][x][1] * psi[t][x2] ;
        else chi[t][x] += 0.5 * eta[t][x][1] * psi[t][x2] ;
      }
      x2 = xdn[x];
      if( field[t][x2] == 0 ){
        if( x2 > x ) chi[t][x] -= 0.5 * eta[t][x][1] * psi[t][x2] ;
        else chi[t][x] += 0.5 * eta[t][x][1] * psi[t][x2] ;
      }
    } else {
      chi[t][x] = psi[t][x]; //An occupied site, matrix becomes identity
    }
  }
}

/* Conjugate gradient for inverting the fermion matrix
 */
void cg_MdM( double inv[NT][NX], double source[NT][NX] )
{
  double rr, pMp, a;
  double r[NT][NX];
  double p[NT][NX];
  double Mp[NT][NX];
  double MMp[NT][NX];
  //vec_zero( inv );

  fM_transpose( Mp, inv );
  fM( MMp, Mp );
  vec_dmul_add( r, source, MMp, -1 );
  //vec_assign( r,source );
  vec_assign( p, r );
  double rr_old = vec_dot( r, r );
  double rr_init = rr_old;
  //printf("CG: %d  %g \n",0,rr);

  if( rr_old < CG_ACCURACY ){
    return ;
  }

  int k;
  for( k = 1; k < CG_MAX_ITER; k++ )
  {
    fM_transpose( Mp, p );
    fM( MMp, Mp );
    pMp = vec_dot( p, MMp ) ;
    a = rr_old / pMp ;
    vec_dmul_add( inv, inv, p, a );
    vec_dmul_add( r, r, MMp, -a );

    rr = vec_dot( r, r );
    //printf("CG: %d  %g %g %g %g\n",k,rr,rr_old,pMp,a);
    if( rr < CG_ACCURACY )
      break;
    if( rr/rr_init > 1e10 ) {
      /* There is an (somewhat) exact zero mode, bail */
      vec_set( inv, 1e50 );
      break;
    }

    double b = rr / rr_old ;
    vec_dmul_add( p, r, p, b );

    rr_old = rr;
  }
  //printf("CG: %d  %g %g %g %g\n",k,rr,rr_old,pMp,a);
}



void cg_propagator( double propagator[NT][NX], double source[NT][NX] )
{
  double tmp[NT][NX];
  vec_zero(tmp);
  fM_transpose( tmp, source );
  cg_MdM( propagator, tmp );
  
}























