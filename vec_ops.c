#include "Thirring.h"

extern int    ***eta;   //Staggered eta matrix
extern double m;
extern double U;
extern double mu;

extern int **field;

/* Neighbour index arrays, to be filled at the beginning
 */
extern int *tup,*xup,*tdn,*xdn;

/* Library functions
 */
void vec_neg(double **a)
{
  for (int t=0; t<NT; t++) for (int x=0; x<NX; x++) 
      a[t][x] = -a[t][x];
}
void vec_zero(double **a)
{
  for (int t=0; t<NT; t++) for (int x=0; x<NX; x++)
      a[t][x] = 0;
}
void vec_one(double **a)
{
  for (int t=0; t<NT; t++) for (int x=0; x<NX; x++)
      a[t][x] = 1;
}
void vec_set(double **a, double d)
{
  for (int t=0; t<NT; t++) for (int x=0; x<NX; x++)
      a[t][x] = d;
}
void vec_d_mul(double **a, double d)
{
  for (int t=0; t<NT; t++) for (int x=0; x<NX; x++)
      a[t][x] = a[t][x]*d;
}
void vec_assign(double **a, double **b)
{
  for (int t=0; t<NT; t++) for (int x=0; x<NX; x++) 
      a[t][x] = b[t][x];
}
void vec_add(double **a, double **b)
{
  for (int t=0; t<NT; t++) for (int x=0; x<NX; x++)  
      a[t][x] += b[t][x];
}
void vec_dmul_add(double **a, double **b, double **d, double e)
{
  for (int t=0; t<NT; t++) for (int x=0; x<NX; x++)  
      a[t][x] = b[t][x] + e * d[t][x] ;
}
double vec_dot(double **a, double **b)
{
  double s = 0 ;
  for (int t=0; t<NT; t++) for (int x=0; x<NX; x++)
      s += a[t][x]*b[t][x];
  return( s );
}

void vec_zero_occupied(double **a) {
  for (int t=0; t<NT; t++) for (int x=0; x<NX; x++)
   if( field[t][x] > 0 ) a[t][x] = 0; 
}

void vec_print_lat(double **a)
{
  for (int t=0; t<NT; t++) { 
    for (int x=0; x<NX; x++)
      printf( " %8.2f ", a[t][x] );
    printf(" \n");
  }
  printf(" \n");
}

double ** alloc_vector(){
  double ** a = malloc(NT*sizeof(double*));
  for (int t=0; t<NT; t++)
    a[t] = malloc(NX*sizeof(double));
  return a;
}

void free_vector(double ** a){
  for (int t=0; t<NT; t++)
    free(a[t]);
  free(a);
}


/* The fermion matrix
 */
#ifdef ANTISYMMETRIC
void fM(double **chi, double **psi )
{
  static double expmu, expmmu;
  static int init=1;
  if(init){
    expmu = exp(mu);
    expmmu = exp(-mu);
    init = 0;
  }
  for (int t=0; t<NT; t++) for (int x=0; x<NX; x++) {
    chi[t][x] = 0;
    if( field[t][x] == 0 ) {
      chi[t][x] = m*psi[t][x];
      int t2=tup[t];
      if( field[t2][x] == 0 ) {
        if( t2 > t ) chi[t][x] += 0.5 * eta[t][x][0] * expmu * psi[t2][x] ;
        else chi[t][x] -= 0.5 * eta[t][x][0] * expmu * psi[t2][x] ;
      }
      t2 = tdn[t];
      if( field[t2][x] == 0 ) {
        if( t2 > t ) chi[t][x] += 0.5 * eta[t][x][0] * expmmu * psi[t2][x] ;
        else chi[t][x] -= 0.5 * eta[t][x][0] * expmmu * psi[t2][x] ;
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

void fM_transpose(double **chi, double **psi)
{
  static double expmu, expmmu;
  static int init=1;
  if(init){
    expmu = exp(mu);
    expmmu = exp(-mu);
    init = 0;
  }
  for (int t=0; t<NT; t++) for (int x=0; x<NX; x++) {
    chi[t][x] = 0;
    if( field[t][x] == 0 ) {
      chi[t][x] = m*psi[t][x];
      int t2=tup[t];
      if( field[t2][x] == 0 ) {
        if( t2 > t ) chi[t][x] -= 0.5 * eta[t][x][0] * expmmu * psi[t2][x] ;
        else chi[t][x] += 0.5 * eta[t][x][0] * expmmu * psi[t2][x] ;
      }
      t2 = tdn[t];
      if( field[t2][x] == 0 ) {
        if( t2 > t ) chi[t][x] -= 0.5 * eta[t][x][0] * expmu * psi[t2][x] ;
        else chi[t][x] += 0.5 * eta[t][x][0] * expmu * psi[t2][x] ;
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
#endif

#ifdef SYMMETRIC
void fM(double **chi, double **psi )
{
  static double expmu, expmmu;
  static int init=1;
  if(init){
    expmu = exp(mu);
    expmmu = exp(-mu);
    init = 0;
  }
  for (int t=0; t<NT; t++) for (int x=0; x<NX; x++) {
    chi[t][x] = 0;
    if( field[t][x] == 0 ) {
      chi[t][x] = m*psi[t][x];
      int t2=tup[t];
      if( field[t2][x] == 0 ) {
        chi[t][x] += 0.5 * eta[t][x][0] * expmu * psi[t2][x] ;
      }
      t2 = tdn[t];
      if( field[t2][x] == 0 ) {
        chi[t][x] -= 0.5 * eta[t][x][0] * expmmu * psi[t2][x] ;
      }
      int x2 = xup[x];
      if( field[t][x2] == 0 ){
        chi[t][x] += 0.5 * eta[t][x][1] * psi[t][x2] ;
      }
      x2 = xdn[x];
      if( field[t][x2] == 0 ){
        chi[t][x] -= 0.5 * eta[t][x][1] * psi[t][x2] ;
      }
    } else {
      chi[t][x] = psi[t][x]; //An occupied site, matrix becomes identity
    }
  }
}

void fM_transpose(double **chi, double **psi)
{
  static double expmu, expmmu;
  static int init=1;
  if(init){
    expmu = exp(mu);
    expmmu = exp(-mu);
    init = 0;
  }
  for (int t=0; t<NT; t++) for (int x=0; x<NX; x++) {
    chi[t][x] = 0;
    if( field[t][x] == 0 ) {
      chi[t][x] = m*psi[t][x];
      int t2=tup[t];
      if( field[t2][x] == 0 ) {
        chi[t][x] -= 0.5 * eta[t][x][0] * expmmu * psi[t2][x] ;
      }
      t2 = tdn[t];
      if( field[t2][x] == 0 ) {
        chi[t][x] += 0.5 * eta[t][x][0] * expmu * psi[t2][x] ;
      }
      int x2 = xup[x];
      if( field[t][x2] == 0 ){
        chi[t][x] -= 0.5 * eta[t][x][1] * psi[t][x2] ;
      }
      x2 = xdn[x];
      if( field[t][x2] == 0 ){
        chi[t][x] += 0.5 * eta[t][x][1] * psi[t][x2] ;
      }
    } else {
      chi[t][x] = psi[t][x]; //An occupied site, matrix becomes identity
    }
  }
}
#endif


double * alloc_field(){
  double * a = malloc(VOLUME*sizeof(double*));
  return a;
}



/* Conjugate gradient for inverting the fermion matrix
 */
void cg_MdM( double **inv, double **source )
{
  double rr, pMp, a;
  double **r = alloc_vector();
  double **p = alloc_vector();
  double **Mp = alloc_vector();
  double **MMp = alloc_vector();
  vec_zero( inv );

  vec_assign( r,source );
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
    fM( Mp, p );
    fM_transpose( MMp, Mp );
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
  free_vector(r);
  free_vector(p);
  free_vector(Mp);
  free_vector(MMp);
  //printf("CG: %d %g %g %g %g\n",k,rr/rr_init,rr_init,pMp,a);
}



void cg_propagator( double **propagator, double **source )
{
  double **tmp = alloc_vector();
  vec_zero(tmp);

  fM_transpose( tmp, source );
  cg_MdM( propagator, tmp );
  fM( tmp, propagator ); //To test cg
  vec_dmul_add( tmp, tmp, source, -1 );
  printf(" test CG %g \n", vec_dot( tmp, tmp ));
  free_vector(tmp);
}



















/* The fermion matrix used in the bozonized version of 
 * the fermion determinant
 */

void fM_occupied(double *chi, double *psi )
{
  static double expmu, expmmu;
  static int init=1;
  if(init){
    expmu = exp(mu);
    expmmu = exp(-mu);
  }
  for (int t=0; t<NT; t++) for (int x=0; x<NX; x++) {
    chi[t*NX+x] = 0;
    if( field[t][x] == 0 ) {
      int t2=tup[t];
      if( field[t2][x] == 0 ) {
        if( t2 > t ) chi[t*NX+x] -= 0.5 * expmmu * eta[t][x][0] * psi[t2*NX+x] ;
        else chi[t*NX+x] += 0.5 * expmmu * eta[t][x][0] * psi[t2*NX+x] ;
      }
      t2 = tdn[t];
      if( field[t2][x] == 0 ) {
        if( t2 > t ) chi[t*NX+x] -= 0.5 * expmu * eta[t][x][0] * psi[t2*NX+x] ;
        else chi[t*NX+x] += 0.5 * expmu * eta[t][x][0] * psi[t2*NX+x] ;
      }
      int x2 = xup[x];
      if( field[t][x2] == 0 ){
        if( x2 > x ) chi[t*NX+x] -= 0.5 * eta[t][x][1] * psi[t*NX+x2] ;
        else chi[t*NX+x] += 0.5 * eta[t][x][1] * psi[t*NX+x2] ;
      }
      x2 = xdn[x];
      if( field[t][x2] == 0 ){
        if( x2 > x ) chi[t*NX+x] -= 0.5 * eta[t][x][1] * psi[t*NX+x2] ;
        else chi[t*NX+x] += 0.5 * eta[t][x][1] * psi[t*NX+x2] ;
      }
    } else {
      chi[t*NX+x] = psi[t*NX+x]; //An occupied site, matrix becomes identity
    }
  }
}

/*
void fM_occupied(double *chi, double *psi )
{
  static double expmu, expmmu;
  static int init=1;
  if(init){
    expmu = exp(mu);
    expmmu = exp(-mu);
  }
  for (int t=0; t<VOLUME; t++) chi[t] = psi[t];
  for (int t=0; t<NT; t++) for (int x=0; x<NX; x++) {
    chi[t*NX+x] = 0;
    int t2=tup[t];
    if( t2 > t ) chi[t*NX+x] -= 0.5 * expmmu * psi[t2*NX+x] ;
    else chi[t*NX+x] += 0.5 * expmmu * psi[t2*NX+x] ;
    
    t2 = tdn[t];
    if( t2 > t ) chi[t*NX+x] -= 0.5 * expmu * psi[t2*NX+x] ;
    else chi[t*NX+x] += 0.5 * expmu * psi[t2*NX+x] ;
    
    int x2 = xup[x];
    if( x2 > x ) chi[t*NX+x] -= 0.5 * eta[t][x][1] * psi[t*NX+x2] ;
    else chi[t*NX+x] += 0.5 * eta[t][x][1] * psi[t*NX+x2] ;
    
    x2 = xdn[x];
    if( x2 > x ) chi[t*NX+x] -= 0.5 * eta[t][x][1] * psi[t*NX+x2] ;
    else chi[t*NX+x] += 0.5 * eta[t][x][1] * psi[t*NX+x2] ;
    
    if( field[t][x] != 0 ) {   //An additional row and column to remove the
                               //determinant at (t,x)
      chi[t*NX+x] -= psi[VOLUME+t*NX+x];
      chi[VOLUME+t*NX+x] = psi[t*NX+x];
    }
  }
}
*/

void fM_occupied_sq(double *chi, double *psi )
{
  double *tmp=alloc_field();
  fM_occupied( tmp, psi );
  fM_occupied( chi, tmp );
  free(tmp);
}

double action(double *psi)
{
  double s = 0 ;
  for (int t=0; t<VOLUME; t++) s += psi[t]*psi[t];
  return( 0.5*s );
}

void vec_gaussian(double *a)
{
  for (int t=0; t<VOLUME; t++){ 
    double x1,x2;
    x1 = mersenne(), x2=mersenne();
    a[t] = sqrt( -2*log(x1) )*cos(2*M_PI*x2);
    if(t<VOLUME-1){
      a[++t] = sqrt( -2*log(x1) )*sin(2*M_PI*x2);
    }
  }
}



int cg_MdM_occupied( double *psi, double *source )
{
  double rr, pMp, a;
  static double *r,*p,*MMp,*inv;
  static int init=1;
  if(init==1){
    r = alloc_field();
    p = alloc_field(); 
    MMp = alloc_field(); 
    inv = alloc_field(); 
    init=0;
  }
  for (int t=0; t<VOLUME; t++) inv[t]=0;
  for (int t=0; t<VOLUME; t++) p[t]=r[t]=source[t];

  double rr_old = 0;
  for (int t=0; t<VOLUME; t++) rr_old+=r[t]*r[t];
  double rr_init = rr_old;
  //printf("CG: %d  %g \n",0,rr_old);

  int k;
  for( k = 1; k < CG_MAX_ITER; k++ )
  {
    fM_occupied_sq( MMp, p );
    pMp=0;
    for (int t=0; t<VOLUME; t++) pMp+=p[t]*MMp[t];
    a = rr_old / pMp ;
    for (int t=0; t<VOLUME; t++){
     inv[t]=inv[t]+a*p[t];
     r[t]=r[t]-a*MMp[t];
    }

    rr=0;
    for (int t=0; t<VOLUME; t++) rr+=r[t]*r[t];
    //printf("CG: %d  %g %g %g %g\n",k,rr,rr_old,pMp,a);
    if( rr < CG_ACCURACY ){
      fM_occupied( psi, inv );
      return(0);
    }
    if( rr/rr_init > 1e10 || isnan(rr) ) {
      /* There is an (somewhat) exact zero mode, bail */
      return(1);
    }

    double b = rr / rr_old ;
    for (int t=0; t<VOLUME; t++) p[t]=r[t]+b*p[t];

    rr_old = rr;
  }
  
  //printf("CG: %d  %g %g %g %g\n",k,rr,rr_old,pMp,a);
  return(1);
}






















