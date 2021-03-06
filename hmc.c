#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mersenne.h"
#include <complex.h>
#include <lapacke.h>


#ifndef M_PI
#define M_PI 3.14159265358979323846264338327 
#endif

/* Lattice size and dimensions */
#define NT 32
#define NX 32

#define ND 2
#define NDIRS (2*ND)

#define TUP 0
#define XUP 1
#define TDN 2
#define XDN 3

#define VOLUME (NT*NX)

//One staggered fermion
#define Nf 2 

//#define ANTISYMMETRIC //Antisymmetric boundaries
//#define SYMMETRIC     //implemented in Thirring_hop
#define OPENX       //open in space, (anti)symmetric in time

#define CG_ACCURACY 1e-30
#define CG_MAX_ITER 100000

/* Simulation parameters */
double m;
double g;
double mu;


/// Neighbour index arrays, to be filled at the beginning
int *tup,*xup,*tdn,*xdn;

/// The gauge field
double ***A;

/// The staggered phase
int    ***eta;   //Staggered eta matrix

/// Get T coordinate of neighbour in direction dir
static inline int tdir(int t, int dir){
  if( dir == TUP ) return tup[t];
  if( dir == TDN ) return tdn[t];
  return(t);
}

/// Get X coordinate of neighbour in direction dir
static inline int xdir(int x, int dir){
  if( dir == XUP ) return xup[x];
  if( dir == XDN ) return xdn[x];
  return(x);
}

/// Opposite of a direction
static inline int opp_dir(int dir){
  return ( dir + ND ) % NDIRS;
}


/// Calculate the gauge action
double calc_gauge_action(double ***A){
  double S = 0;
  for( int t=0;t<NT;t++) for( int x=0;x<NX;x++) for( int dir=0; dir<ND; dir++ ) {
    S+= 1.-cos(A[t][x][dir]);
  }
  return Nf/g*S;
}

/// Update using gauge heatbath ignoring the fermion determinant
void update_puregauge_hb(double ***A){
  for( int t=0;t<NT;t++) for( int x=0;x<NX;x++) for( int dir=0; dir<ND; dir++ ) {
    double s1 = -Nf/g*cos(A[t][x][dir]);
    double new = 2*M_PI*mersenne()-M_PI;
    double s2 = -Nf/g*cos(new);
    double edS = exp(-s2+s1);

    if( mersenne() < edS ) {
      A[t][x][dir] = new;
    }
  }
}


/// Update the gauge field according to a momentum field
double gauge_step(double ***A, double ***mom, double eps){
  for( int t=0;t<NT;t++) for( int x=0;x<NX;x++) for( int dir=0; dir<ND; dir++ ) {
    A[t][x][dir] += 2.0*eps*mom[t][x][dir];
  }
}


///Allocate a pseudofermion vector
_Complex double ** alloc_vector(){
  _Complex double ** v;
  v = (_Complex double **) malloc( NT*sizeof(double *) );
  for (int t=0; t<NT; t++){
    v[t] = (_Complex double *) malloc( NX*sizeof(_Complex double) );
  }
  return v;
}

///Free a pseudofermion vector
void free_vector(_Complex double ** v){
  for (int t=0; t<NT; t++){
    free(v[t]);
  }
  free(v);
}

/// Apply the fermion matrix to a vector
void fm_mul(_Complex double ** v_in, _Complex double ** v_out, double *** A){
  static double expmu, expmmu;
  static int init=1;
  if(init){
    expmu = exp(mu);
    expmmu = exp(-mu);
    init = 0;
  }

  for( int t=0;t<NT;t++) for( int x=0;x<NX;x++){
    // temporary variable for the output on this site
    _Complex double v;

    //Apply the mass term
    v = m*v_in[t][x];

    // Positive time direction
    double sA = sin(A[t][x][0]);
    double cA = cos(A[t][x][0]);
    int t2 = tup[t];
    if( t2 > t ){
      v += 0.5 * (cA + sA*I) * eta[t][x][0] * expmu * v_in[t2][x];
    } else {
      // Antiperiodic boundary condition
      v -= 0.5 * (cA + sA*I) * eta[t][x][0] * expmu * v_in[t2][x];
    }

    // Negative time direction
    t2 = tdn[t];
    sA = sin(A[t2][x][0]);
    cA = cos(A[t2][x][0]);
    if( t2 > t ){
      // Antiperiodic boundary condition
      v += 0.5 * (cA - sA*I) * eta[t2][x][0] * expmmu * v_in[t2][x];
    } else {
      v -= 0.5 * (cA - sA*I) * eta[t2][x][0] * expmmu * v_in[t2][x];
    }

    // Positive x direction
    int x2 = xup[x];
    sA = sin(A[t][x][1]);
    cA = cos(A[t][x][1]);
    if( x2 > x ){
      v += 0.5 * (cA + sA*I) * eta[t][x][1] * v_in[t][x2];
    } else {
      // Antiperiodic boundary condition
      v -= 0.5 * (cA + sA*I) * eta[t][x][1] * v_in[t][x2];
    }

    x2 = xdn[x];
    sA = sin(A[t][x2][1]);
    cA = cos(A[t][x2][1]);
    if( x2 > x ){
      v += 0.5 * (cA - sA*I) * eta[t][x2][1] * v_in[t][x2];
    } else {
      // Antiperiodic boundary condition
      v -= 0.5 * (cA - sA*I) * eta[t][x2][1] * v_in[t][x2];
    }

    v_out[t][x] = v;
  }
}


/// Apply the conjugate fermion matrix to a vector
void fm_conjugate_mul(_Complex double ** v_in, _Complex double ** v_out, double *** A){
  static double expmu, expmmu;
  static int init=1;
  if(init){
    expmu = exp(mu);
    expmmu = exp(-mu);
    init = 0;
  }

  for( int t=0;t<NT;t++) for( int x=0;x<NX;x++){
    // temporary variable for the output on this site
    _Complex double v;

    //Apply the mass term
    v = m*v_in[t][x];

    // Positive time direction
    double sA = sin(A[t][x][0]);
    double cA = cos(A[t][x][0]);
    int t2 = tup[t];
    if( t2 > t ){
      v += 0.5 * (cA + sA*I) * eta[t][x][0] * expmu * v_in[t2][x];
    } else {
      // Antiperiodic boundary condition
      v -= 0.5 * (cA + sA*I) * eta[t][x][0] * expmu * v_in[t2][x];
    }

    // Negative time direction
    t2 = tdn[t];
    sA = sin(A[t2][x][0]);
    cA = cos(A[t2][x][0]);
    if( t2 > t ){
      // Antiperiodic boundary condition
      v += 0.5 * (cA - sA*I) * eta[t2][x][0] * expmmu * v_in[t2][x];
    } else {
      v -= 0.5 * (cA - sA*I) * eta[t2][x][0] * expmmu * v_in[t2][x];
    }

    // Positive x direction
    int x2 = xup[x];
    sA = sin(A[t][x][1]);
    cA = cos(A[t][x][1]);
    if( x2 > x ){
      v += 0.5 * (cA + sA*I) * eta[t][x][1] * v_in[t][x2];
    } else {
      // Antiperiodic boundary condition
      v -= 0.5 * (cA + sA*I) * eta[t][x][1] * v_in[t][x2];
    }

    x2 = xdn[x];
    sA = sin(A[t][x2][1]);
    cA = cos(A[t][x2][1]);
    if( x2 > x ){
      // Antiperiodic boundary condition
      v += 0.5 * (cA - sA*I) * eta[t][x2][1] * v_in[t][x2];
    } else {
      v -= 0.5 * (cA - sA*I) * eta[t][x2][1] * v_in[t][x2];
    }

    v_out[t][x] = v;
  }
}








/// Multiply vector by the fermion matrix squared 
void fmdm_mul(_Complex double ** v_in, _Complex double ** v_out, double *** A){
  _Complex double **tmp=alloc_vector();
  fm_mul( tmp, v_in, A );
  fm_conjugate_mul( v_out, tmp, A );
  free_vector(tmp);
}



/// Construct the full fermion matrix into the array M
void fermion_matrix( _Complex double * M ){
  static double expmu, expmmu;
  static int init=1;
  if(init){
    expmu = exp(mu);
    expmmu = exp(-mu);
    init = 0;
  }

  for( int i=0;i<VOLUME*VOLUME;i++)  M[i] = 0;
    
  for( int t=0;t<NT;t++) for( int x=0;x<NX;x++){
    M[(NX*t + x) + (NX*t+x)*VOLUME] = m;
    int t2 = tup[t];
    int index = (NX*t + x) + (NX*t2+x)*VOLUME;
    double sA = sin(A[t][x][0]);
    double cA = cos(A[t][x][0]);
    if( t2 > t ) M[index] = 0.5 * (cA + sA*I) * eta[t][x][0] * expmu ;
    else M[index] = -0.5 * (cA + sA*I) * eta[t][x][0] * expmu ;
    
    t2 = tdn[t];
    sA = sin(A[t2][x][0]);
    cA = cos(A[t2][x][0]);
    index = (NX*t + x) + (NX*t2+x)*VOLUME;
    if( t2 > t ) M[index] = 0.5 * (cA - sA*I) * eta[t2][x][0] * expmmu ;
    else M[index] = -0.5 * (cA - sA*I) * eta[t2][x][0] * expmmu ;
    
    int x2 = xup[x];
    index = (NX*t + x) + (NX*t+x2)*VOLUME;
    sA = sin(A[t][x][1]);
    cA = cos(A[t][x][1]);
    if( x2 > x ) M[index] = 0.5 * (cA + sA*I) * eta[t][x][1] ;
    else M[index] = -0.5 * (cA + sA*I) * eta[t][x][1] ;
    
    x2 = xdn[x];
    index = (NX*t + x) + (NX*t+x2)*VOLUME;
    sA = sin(A[t][x2][1]);
    cA = cos(A[t][x2][1]);
    if( x2 > x ) M[index] = 0.5 * (cA - sA*I) * eta[t][x2][1] ;
    else M[index] = -0.5 * (cA - sA*I) * eta[t][x2][1] ;
  }
}


/// Calculate the determinant of the fermion matrix (constructed on the fly)
double complex determinant(){
    double complex * M; //The matrix
    M = (_Complex double *)malloc( VOLUME*VOLUME*sizeof(_Complex double) );
    fermion_matrix(M);
    int n=VOLUME;
    int info;
    int *ipiv;
    ipiv = malloc( n*sizeof(int) ); 
    double complex det = 1;
    
    LAPACK_zgetrf( &n, &n, M, &n, ipiv, &info );
    
    if( info > 0 ) printf("error %d\n", info);
    
    for(int i=0; i<n; i++) {
      det *= M[i*n+i];
      if( ipiv[i] != i+1 ) det*=-1; 
    }

    free(M);
    free(ipiv);
    return( det );
}



/// Invert the squared fermion matrix on a vector using CG
void fmdm_invert_cg(_Complex double ** v_in, _Complex double ** v_out, double *** A){
  double rr, pMp, a;
  _Complex double **r = alloc_vector();
  _Complex double **p = alloc_vector();
  _Complex double **Mp = alloc_vector();
  _Complex double **MMp = alloc_vector();
  double rr_old, rr_init;

  rr_old = 0;
  for( int t=0;t<NT;t++) for( int x=0;x<NX;x++){
    v_out[t][x] = 0;
    r[t][x] = v_in[t][x];
    p[t][x] = r[t][x];
    rr_old += creal(r[t][x])*creal(r[t][x])
            + cimag(r[t][x])*cimag(r[t][x]);
  }
  rr_init = rr_old;

  if( rr_old < CG_ACCURACY ){
    return ;
  }

  int k;
  for( k = 1; k < CG_MAX_ITER; k++ )
  {
    fm_mul(p, Mp, A);
    fm_conjugate_mul(Mp, MMp, A);
    pMp = 0;
    for( int t=0;t<NT;t++) for( int x=0;x<NX;x++)
      pMp += creal(p[t][x])*creal(MMp[t][x]) + cimag(p[t][x])*cimag(MMp[t][x]);
    a = rr_old / pMp ;
    for( int t=0;t<NT;t++) for( int x=0;x<NX;x++)
      v_out[t][x] += a*p[t][x];
    for( int t=0;t<NT;t++) for( int x=0;x<NX;x++)
      r[t][x] -= a*MMp[t][x];

    rr = 0;
    for( int t=0;t<NT;t++) for( int x=0;x<NX;x++)
      rr += creal(r[t][x])*creal(r[t][x]) + cimag(r[t][x])*cimag(r[t][x]);

    if( rr < CG_ACCURACY )
      break;
    if( rr/rr_init > 1e10 ) {
      /* There is an (somewhat) exact zero mode, bail */
      printf("Cannot invert fermion matrix\n");
      exit(1);
      break;
    }

    double b = rr / rr_old ;
    for( int t=0;t<NT;t++) for( int x=0;x<NX;x++)
      p[t][x] = r[t][x] + b*p[t][x];

    rr_old = rr;
  }




  free_vector(r);
  free_vector(p);
  free_vector(Mp);
  free_vector(MMp);
}


/// Invert the fermion matrix on a vector using CG
void fm_invert_cg(_Complex double ** v_in, _Complex double ** v_out, double *** A){
  _Complex double **tmp = alloc_vector();
  fm_conjugate_mul( v_in, tmp, A );
  fmdm_invert_cg( tmp, v_out, A);

  free_vector(tmp);
}


/// Generate random pseudofermion
double random_pseudofermion(_Complex double **v, double ***A)
{
  _Complex double **tmp = alloc_vector();
  double action = 0;
  for (int t=0; t<NT; t++) for( int x=0; x<NX; x++){ 
    double x1, x2;
    x1 = mersenne(), x2=mersenne();
    tmp[t][x] = (     sqrt( -2*log(x1) )*cos(2*M_PI*x2)
                  + I*sqrt( -2*log(x1) )*sin(2*M_PI*x2) );

    action += creal(tmp[t][x])*creal(tmp[t][x])
            + cimag(tmp[t][x])*cimag(tmp[t][x]);
  }

  fm_conjugate_mul(tmp, v, A);

  free(tmp);
  return action;
}

/// Generate a vector for stochastic estimation
void stochastic_vector(_Complex double **v)
{
  for (int t=0; t<NT; t++) for( int x=0; x<NX; x++){ 
    double x1, x2;
    x1 = mersenne(), x2=mersenne();
    v[t][x] = (     sqrt( -2*log(x1) )*cos(2*M_PI*x2)
                + I*sqrt( -2*log(x1) )*sin(2*M_PI*x2) );
  }
}


/// Calculate the fermion part of the action
double pseudofermion_action(_Complex double **v, double ***A){
  _Complex double **tmp = alloc_vector();
  double action = 0;
  fmdm_invert_cg(v, tmp, A);

  for( int t=0;t<NT;t++) for( int x=0;x<NX;x++){
    action += creal(v[t][x])*creal(tmp[t][x])
            + cimag(v[t][x])*cimag(tmp[t][x]);
  }
  
  free_vector(tmp);
  return action;
}


/// Action of the stochastic estimate of D^dagger
double stochastic_md_action(_Complex double **v, double ***A){
  _Complex double **tmp = alloc_vector();
  double action = 0;

  fm_conjugate_mul(v, tmp, A);
  for( int t=0;t<NT;t++) for( int x=0;x<NX;x++){
    action += creal(v[t][x]) * creal(tmp[t][x])
            + cimag(v[t][x]) * cimag(tmp[t][x]);
  }

  free_vector(tmp);
  return action;
}


/// Generate a random momentum field
double random_momentum(double ***mom)
{
  double action=0;
  for (int t=0; t<NT; t++){
    for (int x=0; x<NX; x++) {
      double x1, x2;
      x1 = mersenne(), x2=mersenne();
      mom[t][x][0] = sqrt( -2*log(x1) )*cos(2*M_PI*x2);
      mom[t][x][1] = sqrt( -2*log(x1) )*sin(2*M_PI*x2);

      action += mom[t][x][0]*mom[t][x][0] + mom[t][x][1]*mom[t][x][1];
    }
  }


  return action;
}


//#define CHECK_FORCE
/// Calculate the gauge action
double momentum_step(double ***mom, double ***A, _Complex double **psi, _Complex double **st, double eps){

  // Calculate the force of the gauge action and update to mom
  for( int t=0;t<NT;t++) for( int x=0;x<NX;x++) for( int dir=0; dir<ND; dir++ ) {
    //dS/dA = -Nf/g*sin(A)
    mom[t][x][dir] -= eps*Nf/g*sin(A[t][x][dir]);
  }

  // Next the pseudofermion force
  _Complex double **chi = alloc_vector();
  _Complex double **Mchi = alloc_vector();
  fmdm_invert_cg(psi, chi, A);
  fm_mul(chi, Mchi, A);
  for( int t=0;t<NT;t++) for( int x=0;x<NX;x++) {
    double force, dm;
    int t2 = tup[t];
    _Complex double cmc = conj(Mchi[t][x])*chi[t2][x];
    dm = 2 * 0.5 * eta[t][x][0] * exp(mu)
      * ( - creal(cmc) * sin(A[t][x][0]) - cimag(cmc) * cos(A[t][x][0]) );
    if( t2 > t ) force = -dm;
    else force = dm;

    cmc = conj(chi[t][x])*Mchi[t2][x];
    dm =  -2 * 0.5 * eta[t][x][0] * exp(-mu)
      * ( - creal(cmc) * sin(A[t][x][0]) - cimag(cmc) * cos(A[t][x][0]) );
    if( t2 > t ) force -= dm;
    else force += dm;
    mom[t][x][0] -= eps*force;
    
    
    // Test force calculation
    #ifdef CHECK_FORCE
    double f0 = pseudofermion_action(psi, A);
    A[t][x][0] += 0.00001;
    double f1 = pseudofermion_action(psi, A);
    A[t][x][0] -= 0.00001;
    double diff =  force - (f1-f0)/0.00001;
    if( diff*diff > 0.0001 ){
      printf("T-force at %d %d incorrect!\n",t,x);
      printf("T calculated force is %g\n", force);
      printf("T real derivative is %g\n", (f1-f0)/0.00001);
      printf("diff = %g\n", diff);
      printf("%g %g %g\n", A[t][x][0], sin(A[t][x][0]), cos(A[t][x][0]));
      cmc = conj(Mchi[t][x])*chi[t2][x];
      printf("%g %g \n", 
        2 * 0.5 * eta[t][x][0] * exp(mu) *(  creal(cmc) * sin(A[t][x][0]) ),
        2 * 0.5 * eta[t][x][0] * exp(mu) *(- cimag(cmc) * cos(A[t][x][0]) )
      );
      cmc = conj(chi[t][x])*Mchi[t2][x];
      printf("%g %g \n", 
      - 2 * 0.5 * eta[t][x][0] * exp(-mu)*(  creal(cmc) * sin(A[t][x][0]) ),
      - 2 * 0.5 * eta[t][x][0] * exp(-mu)*(- cimag(cmc) * cos(A[t][x][0]) )
      );
      printf("\n");
    }
    #endif
    


    force = 0;
    int x2 = xup[x];
    cmc = conj(Mchi[t][x])*chi[t][x2];
    dm = 2 * 0.5 * eta[t][x][1]
      * ( - creal(cmc) * sin(A[t][x][1]) - cimag(cmc) * cos(A[t][x][1]));
    if( x2 > x ) force = -dm;
    else force = dm;

    cmc = conj(chi[t][x])*Mchi[t][x2];
    dm = - 2 * 0.5 * eta[t][x][1]
      * ( - creal(cmc) * sin(A[t][x][1]) - cimag(cmc) * cos(A[t][x][1]));
    if( x2 > x ) force -= dm;
    else force += dm;
    mom[t][x][1] -= eps*force;
    
    
    // Test force calculation
    #ifdef CHECK_FORCE
    f0 = pseudofermion_action(psi, A);
    A[t][x][1] += 0.000001;
    f1 = pseudofermion_action(psi, A);
    A[t][x][1] -= 0.000001;
    diff =  force - (f1-f0)/0.000001;
    if( diff*diff > 0.0001 ){
      printf("X-force at %d %d incorrect!\n",t,x);
      printf("X calculated force is %g\n", force);
      printf("X real derivative is %g\n", (f1-f0)/0.000001);
      printf("diff = %g\n", diff);
      printf("%g %g %g\n", A[t][x][1], sin(A[t][x][1]), cos(A[t][x][1]));
      cmc = conj(Mchi[t][x])*chi[t][x2];
      printf("%g %g \n",
        2 * 0.5 * eta[t][x][1] * exp(mu) *(  creal(cmc) * sin(A[t][x][1])  ),
        2 * 0.5 * eta[t][x][1] * exp(mu) *(- cimag(cmc) * cos(A[t][x][1])  )
      );
      cmc = conj(chi[t][x])*Mchi[t][x2];
      printf("%g %g \n", 
      - 2 * 0.5 * eta[t][x][1] * exp(-mu)*(- creal(cmc) * sin(A[t][x][1])  ),
      - 2 * 0.5 * eta[t][x][1] * exp(-mu)*(- cimag(cmc) * cos(A[t][x][1])  )
      );
      printf("\n");
    }
    #endif
  }
  
  // And the force from the matrix M_conjugate
  for( int t=0;t<NT;t++) for( int x=0;x<NX;x++) {
    double force, dm;
    int t2 = tup[t];
    _Complex double cmc = conj(st[t][x])*st[t2][x];
    dm = - 0.5 * eta[t][x][0] * exp(-mu)
       * ( creal(cmc) * sin(A[t][x][0]) + cimag(cmc) * cos(A[t][x][0]) );
    if( t2 > t ) force = -dm;
    else force = dm;

    dm = 0.5 * eta[t][x][0] * exp(mu)
       * ( creal(cmc) * sin(A[t][x][0]) + cimag(cmc) * cos(A[t][x][0]) );
    if( t2 > t ) force -= dm;
    else force += dm;
    mom[t][x][0] -= eps*force;

    #ifdef CHECK_FORCE
    double f0 = stochastic_md_action(st, A);
    A[t][x][0] += 0.00001;
    double f1 = stochastic_md_action(st, A);
    A[t][x][0] -= 0.00001;
    double diff =  force - (f1-f0)/0.00001;
    if( diff*diff > 0.0001 ){
      printf("T-force at %d %d incorrect!\n",t,x);
      printf("T calculated force is %g\n", force);
      printf("T real derivative is %g\n", (f1-f0)/0.00001);
      printf("diff = %g\n", diff);
      printf("%g %g %g\n", A[t][x][0], sin(A[t][x][0]), cos(A[t][x][0]));
      printf("%g %g \n", 
        0.5 * eta[t][x][0] * exp(-mu) *(  creal(cmc) * sin(A[t][x][0]) ),
        0.5 * eta[t][x][0] * exp(-mu) *(- cimag(cmc) * cos(A[t][x][0]) )
      );
      printf("%g %g \n", 
      - 0.5 * eta[t][x][0] * exp(mu)*(  creal(cmc) * sin(A[t][x][0]) ),
      - 0.5 * eta[t][x][0] * exp(mu)*(- cimag(cmc) * cos(A[t][x][0]) )
      );
      printf("\n");
    }
    #endif

    force = 0;
    int x2 = xup[x];
    cmc = conj(st[t][x])*st[t2][x];
    dm = - 0.5 * eta[t][x][1] 
       * ( creal(cmc) * sin(A[t][x][1]) + cimag(cmc) * cos(A[t][x][1]) );
    if( x2 > x ) force = -dm;
    else force = dm;

    cmc = conj(st[t][x])*st[t2][x];
    dm = 0.5 * eta[t][x][1] 
       * ( creal(cmc) * sin(A[t][x][1]) + cimag(cmc) * cos(A[t][x][1]) );
    if( x2 > x ) force -= dm;
    else force += dm;
    mom[t][x][1] -= eps*force;
  }
  
  free_vector(chi);
  free_vector(Mchi);
}




/// Update the gauge using HMC
void update_gauge(double ***A){

  static int init = 1;
  static _Complex double **psi = NULL;
  static _Complex double **chi = NULL;
  static double ***mom = NULL;
  static double ***new_A = NULL;
  if(init){
    psi = alloc_vector();
    chi = alloc_vector();
    mom = malloc( NT*sizeof(double **) );
    new_A = malloc( NT*sizeof(double **) );
    for (int t=0; t<NT; t++){
      mom[t] = malloc( NX*sizeof(double *) );
      new_A[t] = malloc( NX*sizeof(double *) );
      for (int x=0; x<NX; x++) {
        mom[t][x] = malloc( ND*sizeof(double) );
        new_A[t][x] = malloc( ND*sizeof(double) );
      }
    }

    init = 0;
  }

  double mdm_action = random_pseudofermion(psi, A);
  double momentum_action = random_momentum(mom);
  double gauge_action = calc_gauge_action(A);
  stochastic_vector(chi);
  double md_action = stochastic_md_action(chi, A);

  printf("Start HMC: Sg %g, Smdm %g, Smd %g, Smom %g\n", gauge_action, mdm_action, md_action, momentum_action);

  for (int t=0; t<NT; t++) for (int x=0; x<NX; x++) for (int d=0; d<ND; d++)  {
    new_A[t][x][d] = A[t][x][d];
  }


  int nsteps = 10;
  double traj_length=1;
  for(int i=0; i<nsteps;i++){
    
    gauge_step(new_A,mom,traj_length*0.5/nsteps);
    momentum_step(mom, new_A, psi, chi, traj_length/nsteps);
    gauge_step(new_A,mom,traj_length*0.5/nsteps);

  }


  double mdm_action2 = pseudofermion_action(psi, new_A);

  double momentum_action2 = 0;
  for (int t=0; t<NT; t++) for (int x=0; x<NX; x++) for (int d=0; d<ND; d++) {
      momentum_action2 += mom[t][x][d]*mom[t][x][d];
  }
  double gauge_action2 = calc_gauge_action(new_A);
  double md_action2 = stochastic_md_action(chi, new_A);


  double dS = momentum_action2 - momentum_action 
            + mdm_action2      - mdm_action
            + md_action2       - md_action
            + gauge_action2    - gauge_action;


  printf("HMC End, dS %g, Sg %g, Smdm %g, Smd %g, Sm %g\n", dS, gauge_action2, mdm_action2, md_action2, momentum_action2);


  if (exp(-dS) > mersenne() ){
    printf("HMC ACCEPTED\n");
    for (int t=0; t<NT; t++) for (int x=0; x<NX; x++) for (int d=0; d<ND; d++)
      A[t][x][d] = new_A[t][x][d];
  } else {
    printf("HMC REJECTED\n");
  }

}











/// Test that fm_conjugate_mul applies the conjugate of fm_mul
void test_conjugate(double ***A){
  _Complex double **c = alloc_vector();
  _Complex double **mc = alloc_vector();
  _Complex double **mdc = alloc_vector();

  stochastic_vector(c);

  fm_mul(c, mc, A);
  fm_conjugate_mul(c, mdc, A);

  _Complex double cmdc = 0;
  _Complex double cmc  = 0;

  for( int t=0;t<NT;t++) for( int x=0;x<NX;x++){
    cmdc += conj(c[t][x])*mdc[t][x];
    cmc  += conj(mc[t][x])*c[t][x];
  }

  _Complex double cdiff = cmdc - conj(cmc);

  double diff = creal(cdiff)*creal(cdiff) + cimag(cdiff)*cimag(cdiff);
  if(diff > 0.001){
    printf("ERROR fm_conjugate_mul is not the conjugate of fm_mul\n");
    printf("Difference = %g\n", diff);
    exit(1);
  }

  free_vector(c);
  free_vector(mc);
  free_vector(mdc);
}



/// Calculate the phase of the fermion action
double fermion_phase(double ***A){
  _Complex double **c = alloc_vector();
  _Complex double **mdc = alloc_vector();

  int sources = 20;
  double imaginary_part = 0;

  for( int i=0;i<sources;i++) {
    stochastic_vector(c);
    fm_conjugate_mul(c, mdc, A);

    for( int t=0;t<NT;t++) for( int x=0;x<NX;x++){
      imaginary_part += creal(c[t][x]) * cimag(mdc[t][x])
                      - cimag(c[t][x]) * creal(mdc[t][x]);
    }
  }

  free_vector(c);
  free_vector(mdc);

  return imaginary_part/sources;
}



/// Run measurements
/// * Magnetisation: average gauge field
/// * Determinant: the fermion determinant
/// * Phase: The imaginary part of the action 
void measure(){
  static int n=0;
  static double M = 0;
  static double complex det = 0;
  static double phase = 0;

  test_conjugate(A);

  for( int t=0;t<NT;t++) for( int x=0;x<NX;x++) for( int dir=0; dir<ND; dir++) {
    M += A[t][x][dir];
  }

  phase = fermion_phase(A);
  
  n++;
  
  printf("Magnetisation %g\n",M/VOLUME);
  printf("Phase %g\n",phase);
  phase = 0; det = 0; M = 0;
}




int main(void){
    
  #ifdef DEBUG
  feenableexcept(FE_INVALID | FE_OVERFLOW);
  #endif 

  long seed=0;
  int n_loops=1, n_measure=1;
    
  /* Read in the input */
  printf(" Number of updates : ");
  scanf("%d",&n_loops);

  printf(" Updates / measurement : ");
  scanf("%d",&n_measure);
  
  printf(" m : ");
  scanf("%lf",&m);

  printf(" g : ");
  scanf("%lf",&g);

  printf(" mu : ");
  scanf("%lf",&mu);

  printf(" Random number : ");
  scanf("%ld",&seed);
  seed_mersenne( seed );
  
  /* "Warm up" the rng generator */
  for ( int i=0; i<543210; i++) mersenne();

  printf(" \n++++++++++++++++++++++++++++++++++++++++++\n");
  printf(" 2D quenched Thirring model, ( %d , %d ) lattice\n", NT, NX );
  printf(" %d updates per measurements\n", n_measure );
  printf(" m %f \n", m);
  printf(" g %f \n", g);
  printf(" mu %f \n", mu);
  printf(" Random seed %ld\n", seed );
  
  // Allocate the fields and utility arrays
  eta = malloc( NT*sizeof(int *) );
  A = malloc( NT*sizeof(double *) );
  for (int t=0; t<NT; t++){
    A[t] = malloc( (NX+1)*sizeof(double *) );
    eta[t] = malloc( (NX+1)*sizeof(int *) );
    for (int x=0; x<NX+1; x++) {
      eta[t][x] = malloc( ND*sizeof(int) );
      A[t][x] = malloc( ND*sizeof(double) );
      A[t][x][0] = A[t][x][1] = 0;
    }
  }
  xup = malloc( (NX+1)*sizeof(int) );
  xdn = malloc( (NX+1)*sizeof(int) );
  tup = malloc( NT*sizeof(int) );
  tdn = malloc( NT*sizeof(int) );
  
  /* fill up the index array and the staggered matrix */
  for ( int i=0; i<NT; i++ ) {
    tup[i] = (i+1) % NT;
    tdn[i] = (i-1+NT) % NT;
  }
  for (int i=0; i<NX+1; i++) {
    xup[i] = (i+1) % NX;
    xdn[i] = (i-1+NX) % NX;
  }

  for (int t=0; t<NT; t++) for (int x=0; x<NX; x++) {
    for( int dir=0; dir<ND; dir++ ) A[t][x][dir] = 0;
    eta[t][x][1] = 1;
    if( x%2 == 0 ){
      eta[t][x][0] = 1;
    } else {
      eta[t][x][0] = -1;
    }
#ifdef OPENX
    eta[t][NX][0] = eta[t][NX][1] = 0;
#endif
  }

  for (int i=0; i<100; i++) {
    update_puregauge_hb(A);
  }

  // Run updates and measurements
  for (int i=1; i<n_loops+1; i++) {
    update_gauge(A);
        
    if((i%n_measure)==0){
      /* Statistics */
      measure();
    }
  }
  
  
  return 0;
}















