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
#define NT 48
#define NX 6

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
#define CG_MAX_ITER 10000

/* Simulation parameters */
double m;
double g;
double mu;


/* Neighbour index arrays, to be filled at the beginning
 */
int *tup,*xup,*tdn,*xdn;

double ***A;
int    ***eta;   //Staggered eta matrix

/* Functions for fetching neighboring coordinates */
static inline int tdir(int t, int dir){
  if( dir == TUP ) return tup[t];
  if( dir == TDN ) return tdn[t];
  return(t);
}

static inline int xdir(int x, int dir){
  if( dir == XUP ) return xup[x];
  if( dir == XDN ) return xdn[x];
  return(x);
}

/* Opposite of a direction */
static inline int opp_dir(int dir){
  return ( dir + ND ) % NDIRS;
}

double gauge_action(){
    double S = 0;
    for( int t=0;t<NT;t++) for( int x=0;x<NX;x++) for( int dir=0; dir<ND; dir++ ) {
        S+= 1.-cos(A[t][x][dir]);
    }
    return Nf/g*S;
}

void update_gauge(){
    for( int t=0;t<NT;t++) for( int x=0;x<NX;x++) for( int dir=0; dir<ND; dir++ ) {
        double s1 = -Nf/g*cos(A[t][x][dir]);
        double new = 2*M_PI*mersenne()-M_PI;
        double s2 = -Nf/g*cos(new);
        double edS = exp(-s2+s1);
        //printf("New A=%f at (%d,%d), direction %d, exp(-deltaS)=%g\n", A[t][x][dir],t,x,dir, edS);

        if( mersenne() < edS ) {
            A[t][x][dir] = new;
            //printf("Accepted\n");
        }
    }
}










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



double complex determinant(){
    double complex * M; //The matrix
    M = (_Complex double *)malloc( VOLUME*VOLUME*sizeof(_Complex double) );
    fermion_matrix(M);
    int n=VOLUME;
    int info;
    int *ipiv;
    ipiv = malloc( n*sizeof(int) ); 
    double complex det = 1;
    
    /*for( int i=0;i<VOLUME;i++) { for( int j=0;j<VOLUME;j++){
        printf(" %4.2f ", creal( M[i+j*VOLUME] ) );
      }
      printf("\n");
    }
    printf("\n");
    */
    
    LAPACK_zgetrf( &n, &n, M, &n, ipiv, &info );
    
    if( info > 0 ) printf("error %d\n", info);
    
    for(int i=0; i<n; i++) {
      det *= M[i*n+i];
      if( ipiv[i] != i+1 ) det*=-1; 
    }
    
    /*
    int lwork=n*n;
    _Complex double *work;
    work = malloc( lwork*sizeof(_Complex double) );

    LAPACK_zgetri(&n, M, &n, ipiv, work, &lwork, &info);
        free(work);

    for( int i=0;i<VOLUME;i++) { for( int j=0;j<VOLUME;j++){
        printf(" %4.2f ", creal( M[i+j*VOLUME] ) );
      }
      printf("\n");
    }
    printf("\n");
    */
    
    free(M);
    free(ipiv);
    return( det );
}















int n_average;
void measure(){
  /* Measurements */
  static int n=0;
  static double M = 0;
  static double complex det = 0;
  static double sign = 1;

  double complex this_det = determinant();
  det += this_det;
  double this_sign = (creal(this_det) > 0) ? 1 : -1;
  sign += this_sign;
  
  for( int t=0;t<NT;t++) for( int x=0;x<NX;x++) for( int dir=0; dir<ND; dir++) {
    M += this_sign*A[t][x][dir];
  }
  
  n++;
  
  if(n%n_average == 0){
     printf("Magnetisation %g\n",M/(n_average*VOLUME));
     printf("Determinant Re %g Im %g\n",creal(det)/(n_average),cimag(det)/(n_average));
     printf("Sign %g\n",sign/(n_average));
     sign = 0; det = 0; M = 0;
  }

}



int main(void){
    
  #ifdef DEBUG
  feenableexcept(FE_INVALID | FE_OVERFLOW);
  #endif 

  long seed;
  int n_loops, n_measure;
    
  /* Read in the input */
  printf(" Number of updates : ");
  scanf("%d",&n_loops);

  printf(" Updates / measurement : ");
  scanf("%d",&n_measure);
  
  printf(" Average over configurations : ");
  scanf("%d",&n_average);

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
    
  eta = malloc( NT*sizeof(int *) );
  A = malloc( NT*sizeof(double *) );
  for (int t=0; t<NT; t++){
    A[t] = malloc( (NX+1)*sizeof(double *) );
    eta[t] = malloc( (NX+1)*sizeof(int *) );
    for (int x=0; x<NX+1; x++) {
      eta[t][x] = malloc( ND*sizeof(int) );
      A[t][x] = malloc( ND*sizeof(double) );
    }
  }
  xup = malloc( (NX+1)*sizeof(int) );
  xdn = malloc( (NX+1)*sizeof(int) );
  tup = malloc( NT*sizeof(int) );
  tdn = malloc( NT*sizeof(int) );
  
  /* fill up the index array */
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

  for (int i=1; i<n_loops+1; i++) {
    update_gauge();
    
    //printf("Updated gauge\n");
    //printf("Action %g\n", gauge_action());
    
    if((i%n_measure)==0){
      /* Statistics */
      measure();
    }
  }
  
  
  return 0;
}















