/****************************************************************
 * Simulate the thirring model using the fermion bag algorithm 
 * (arXiv:0910.5736). The mass term is represented as a field of 
 * monomers (occupied sites) and the four fermion term is  
 * represented as dimers (occupied links). 
 *
 ****************************************************************/
#ifdef DEBUG
#include <fenv.h>
#endif

#include "Thirring.h"

/* storage */
extern int    ***eta;   //Staggered eta matrix
extern double *Dinv;             //Storage for even to odd propagator
extern double m;
extern double mu;

/* Neighbour index arrays, to be filled at the beginning
 */
int *tup,*xup,*tdn,*xdn;


/* The fermion matrix
 */
#ifdef ANTISYMMETRIC //Antisymmetrix boundaries
double fM_index( int t1, int x1, int t2, int x2 )
{
  if(t1==t2 ){ 
#ifndef WITH_MASS_MONOMERS
    if(x1==x2){
      return m;
    } else 
#endif
    if( x2 == xup[x1] || x2 == xdn[x1] ) {
#if NX==2
      if (x2>x1) { return eta[t1][x1][1] ; }
      else { return -eta[t1][x1][1] ; }
#else
      if (x2>x1) { return 0.5*eta[t1][x1][1] ; }
      else { return -0.5*eta[t1][x1][1] ; }
#endif
    }
    else  return( 0.0 );
  }
  else if ( x1==x2  ) {
    if( t2 == tup[t1] ){
#if NT==2
      if (t2>t1) { return (exp(mu)+exp(-mu))*eta[t1][x1][0] ; }
      else { return -(exp(mu)+exp(-mu))*eta[t1][x1][0] ; }
#else
      if (t2>t1) { return 0.5*exp(mu)*eta[t1][x1][0] ; }
      else { return -0.5*exp(mu)*eta[t1][x1][0] ; }
#endif
    } else if ( t2 == tdn[t1] ) {
#if NT>2
      if (t2>t1) { return 0.5*exp(-mu)*eta[t1][x1][0] ; }
      else { return -0.5*exp(-mu)*eta[t1][x1][0] ; }
#endif
    }
    else return( 0.0 );
  } else return( 0.0 );
}
#endif

#ifdef SYMMETRIC //Symmetrix boundaries
double fM_index( int t1, int x1, int t2, int x2 )
{
  if(t1==t2 ){
#ifndef WITH_MASS_MONOMERS
    if(x1==x2){
      return m;
    } else
#endif
    if( x2 == xup[x1] ) {
#if NX==2
      return 0 ;
#else
      return 0.5*eta[t1][x1][1] ;

    } else if( x2 == xdn[x1] ) {
      return -0.5*eta[t1][x1][1] ;
#endif
    }
    else return( 0.0 );

  }
  else if ( x1==x2  ) {
  
   if( t2 == tup[t1] ){
#if NT==2
      if (t2>t1) { return (exp(mu)+exp(-mu))*eta[t1][x1][0] ; }
      else { return -(exp(mu)+exp(-mu))*eta[t1][x1][0] ; }
#else
      if (t2>t1) { return 0.5*exp(mu)*eta[t1][x1][0] ; }
      else { return -0.5*exp(mu)*eta[t1][x1][0] ; }
#endif
    } else if ( t2 == tdn[t1] ) {
#if NT>2
      if (t2>t1) { return 0.5*exp(-mu)*eta[t1][x1][0] ; }
      else { return -0.5*exp(-mu)*eta[t1][x1][0] ; }
#endif
    }
    else return( 0.0 );


  } else return( 0.0 );
}
#endif


#ifdef FULLSYMMETRIC //Symmetrix boundaries
double fM_index( int t1, int x1, int t2, int x2 )
{
  if(t1==t2 ){ 
#ifndef WITH_MASS_MONOMERS
    if(x1==x2){
      return m;
    } else 
#endif
    if( x2 == xup[x1] ) {
#if NX==2
      return 0 ;
#else
      return 0.5*eta[t1][x1][1] ;

    } else if( x2 == xdn[x1] ) {
      return -0.5*eta[t1][x1][1] ;
#endif
    }
    else return( 0.0 );
  }
  else if ( x1==x2  ) {
    if( t2 == tup[t1] ){
#if NT==2
      return (exp(mu)-exp(-mu))*eta[t1][x1][0] ; 
#else
      return 0.5*exp(mu)*eta[t1][x1][0] ;
#endif
    } else if ( t2 == tdn[t1] ) {
#if NT>2
      return -0.5*exp(-mu)*eta[t1][x1][0] ; 
#endif
    }
    else return( 0.0 );
  } else return( 0.0 );
}
#endif

#ifdef OPENX //Open in space, antisymmetric in time
double fM_index( int t1, int x1, int t2, int x2 )
{
 if( x1!=NX && x2!=NX ) {
  if(t1==t2 ){ 
#ifndef WITH_MASS_MONOMERS
    if(x1==x2){ return m; }
    else
#endif
    if( x2 == xup[x1]) { return 0.5*eta[t1][x1][1] ; }
    else if( x2 == xdn[x1] ){ return -0.5*eta[t1][x1][1] ; }
    else return( 0.0 );
  }
  else if ( x1==x2  ) {
    if( t2 == tup[t1] ){
#if NT==2
      if (t2>t1) { return (exp(mu)+exp(-mu))*eta[t1][x1][0] ; }
      else { return -(exp(mu)+exp(-mu))*eta[t1][x1][0] ; }
#else
      if (t2>t1) { return 0.5*exp(mu)*eta[t1][x1][0] ; }
      else { return -0.5*exp(mu)*eta[t1][x1][0] ; }
#endif
    } else if ( t2 == tdn[t1] ) {
#if NT>2
      if (t2>t1) { return 0.5*exp(-mu)*eta[t1][x1][0] ; }
      else { return -0.5*exp(-mu)*eta[t1][x1][0] ; }
#endif
    }
    else return( 0.0 );
  } else return( 0.0 );
 } else return( 0.0 );
}
#endif



/* Calculate the propagator matrix
 */
#ifdef WITH_MASS_MONOMERS
void calc_Dinv( )
{
  int n=VOLUME/2;
  int *ipiv;
  int info;
  int lwork=n*n;
  double *work;
  double *M;

  ipiv = malloc( n*sizeof(int) );

  double **source, **propagator;
  source = alloc_vector();
  propagator = alloc_vector();

  FILE * Dinv_file;
  char filename[100];
#ifdef ANTISYMMETRIC
  sprintf(filename, "free_propagator_eo_T%dX%d_mu%0.6f",NT,NX,mu);
#endif
#ifdef SYMMETRIC
  sprintf(filename, "free_propagator_eo_T%dX%d_mu%0.6f_sym",NT,NX,mu);
#endif
#ifdef OPENX
  sprintf(filename, "free_propagator_eo_T%dX%d_mu%0.6f_open",NT,NX,mu);
#endif
  
  Dinv_file = fopen(filename,"rb");
  if (Dinv_file){
    fread(Dinv, VOLUME*VOLUME/2, sizeof(double), Dinv_file);
    fclose(Dinv_file);
  } else {
  
    struct timeval start, end;
    gettimeofday(&start,NULL);

    // Construct the full Volume to Volume Dirac matrix
    // Odd to even here, the inverse will be even to odd
    //
    work = malloc( lwork*sizeof(double) );
    M = malloc( lwork*sizeof(double) );
    for (int t1=0; t1<NT; t1++) for (int x1=0; x1<NX; x1++)  if((t1+x1)%2==0) {
      vec_zero( source ); vec_zero( propagator );
      source[t1][x1] = 1;
      cg_propagator(propagator, source );
      //vec_zero( source ); 
      //fM(source, propagator );
      int i1 = (NX*t1 + x1)/2;
      for (int t2=0; t2<NT; t2++) for (int x2=0; x2<NX; x2++) if((t2+x2)%2==1) {
        int i2 = (NX*t2 + x2)/2;
        M[i2*n + i1] = propagator[t2][x2];
      }
    }

    // LU decompose
    /*LAPACK_dgetrf( &n, &n, M, &n, ipiv, &info );
    if( info != 0 ) {
      printf("calc_Dinv_e: sgetrf returned an error %d! \n", info);
      exit(-1);
    }

    LAPACK_dgetri(&n, M, &n, ipiv, work, &lwork, &info);
    if( info != 0 ) {
      printf("calc_Dinv_e: sgetri returned an error %d! \n", info);
      exit(-1);
    }
    */
    for(int i = 0; i < n*n; i++) Dinv[n*n+i] = M[i];
    /*for (int t1=0; t1<NT; t1++) for (int x1=0; x1<NX; x1++) if((t1+x1)%2==0) {
      vec_zero( source ); vec_zero( propagator );
      int i1 = (NX*t1 + x1)/2;
      for (int t2=0; t2<NT; t2++) for (int x2=0; x2<NX; x2++) if((t2+x2)%2==0) {
        int i2 = (NX*t2 + x2)/2;
        source[t2][x2] = M[i1*n + i2] ;
      }
      fM_transpose(propagator, source );
      for (int t2=0; t2<NT; t2++) for (int x2=0; x2<NX; x2++) if((t2+x2)%2==1) {
        int i2 = (NX*t2 + x2)/2;
        Dinv[i2*n + i1] = -propagator[t2][x2];
      }
    }
    printf("\n");
    for (int t1=0; t1<NT; t1++) for (int x1=0; x1<NX; x1++) if((t1+x1)%2==0) {
      int i1 = (NX*t1 + x1)/2;
      for (int t2=0; t2<NT; t2++) for (int x2=0; x2<NX; x2++) if((t2+x2)%2==1) {
        int i2 = (NX*t2 + x2)/2;
        printf(" %6.3f ", Dinv[i1*n + i2]);
      }
      printf("\n");
    }*/


    // Even to odd
    for (int t1=0; t1<NT; t1++) for (int x1=0; x1<NX; x1++) if((t1+x1)%2==1) {
      vec_zero( source ); vec_zero( propagator );
      source[t1][x1] = 1;
      cg_propagator(propagator, source );
      //vec_zero( source ); 
      //fM(source, propagator );
      int i1 = (NX*t1 + x1)/2;
      for (int t2=0; t2<NT; t2++) for (int x2=0; x2<NX; x2++) if((t2+x2)%2==0) {
        int i2 = (NX*t2 + x2)/2;
        M[i2*n + i1] = propagator[t2][x2];
      }
    }
    /*
    gettimeofday(&start,NULL);
    // LU decompose
    LAPACK_dgetrf( &n, &n, M, &n, ipiv, &info );
    if( info != 0 ) {
      printf("calc_Dinv: sgetrf returned an error %d! \n", info);
      exit(-1);
    }

    LAPACK_dgetri(&n, M, &n, ipiv, work, &lwork, &info);
    if( info != 0 ) {
      printf("calc_Dinv: sgetri returned an error %d! \n", info);
      exit(-1);
    }
    */
    for(int i = 0; i < n*n; i++) Dinv[i] = M[i];
    /*for (int t1=0; t1<NT; t1++) for (int x1=0; x1<NX; x1++) if((t1+x1)%2==1) {
      vec_zero( source ); vec_zero( propagator );
      int i1 = (NX*t1 + x1)/2;
      for (int t2=0; t2<NT; t2++) for (int x2=0; x2<NX; x2++) if((t2+x2)%2==1) {
        int i2 = (NX*t2 + x2)/2;
        source[t2][x2] = M[i1*n + i2] ;
      }
      fM_transpose(propagator, source );
      for (int t2=0; t2<NT; t2++) for (int x2=0; x2<NX; x2++) if((t2+x2)%2==0) {
        int i2 = (NX*t2 + x2)/2;
        Dinv[n*n + i2*n + i1] = -propagator[t2][x2];
      }
    }

    printf("\n");
    for (int t1=0; t1<NT; t1++) for (int x1=0; x1<NX; x1++) if((t1+x1)%2==0) {
      int i1 = (NX*t1 + x1)/2;
      for (int t2=0; t2<NT; t2++) for (int x2=0; x2<NX; x2++) if((t2+x2)%2==1) {
        int i2 = (NX*t2 + x2)/2;
        printf(" %6.3f ", Dinv[n*n + i1*n + i2]);
      }
      printf("\n");
    }*/


    free(work);
    free(M);

    gettimeofday(&end,NULL);
    int diff = 1e6*(end.tv_sec-start.tv_sec) + end.tv_usec-start.tv_usec;

    printf("Inverted fermion matrix in %.3g seconds\n", 1e-6*diff);
    Dinv_file = fopen(filename,"wb");
    if (Dinv_file){
      fwrite(Dinv, VOLUME*VOLUME/2, sizeof(double), Dinv_file);
      fclose(Dinv_file);
    }
  }
  free(ipiv);
  free_vector(source);
  free_vector(propagator);
}

#else

void calc_Dinv( )
{
  int n=VOLUME;
  int *ipiv;
  int info;
  int lwork=n*n;
  double *work;
  double *M;

  ipiv = malloc( n*sizeof(int) );

  FILE * Dinv_file;
  char filename[100];
#ifdef ANTISYMMETRIC
  sprintf(filename, "free_propagator_T%dX%d_m%.6f_mu%0.6f",NT,NX,m,mu);
#endif
#ifdef SYMMETRIC
  sprintf(filename, "free_propagator_T%dX%d_m%.6f_mu%0.6f_sym",NT,NX,m,mu);
#endif
#ifdef OPENX
  sprintf(filename, "free_propagator_T%dX%d_m%.6f_mu%0.6f_open",NT,NX,m,mu);
#endif
  
  Dinv_file = fopen(filename,"rb");
  if (Dinv_file){
    fread(Dinv, VOLUME*VOLUME, sizeof(double), Dinv_file);
    fclose(Dinv_file);
  } else {
  
    struct timeval start, end;
    gettimeofday(&start,NULL);

    // Construct the full Volume to Volume Dirac matrix
    work = malloc( lwork*sizeof(double) );
    M = malloc( lwork*sizeof(double) );
    for (int t1=0; t1<NT; t1++) for (int x1=0; x1<NX; x1++) {
      int i1 = (NX*t1 + x1);
      for (int t2=0; t2<NT; t2++) for (int x2=0; x2<NX; x2++) {
        int i2 = (NX*t2 + x2);
        M[i1*n + i2] = fM_index( t1, x1, t2, x2 );
      }
    }

    // LU decompose
    LAPACK_dgetrf( &n, &n, M, &n, ipiv, &info );
    if( info != 0 ) {
      printf("calc_Dinv: sgetrf returned an error %d! \n", info);
      exit(-1);
    }

    LAPACK_dgetri(&n, M, &n, ipiv, work, &lwork, &info);
    if( info != 0 ) {
      printf("calc_Dinv: sgetri returned an error %d! \n", info);
      exit(-1);
    }

    for(int i = 0; i < n*n; i++) Dinv[i] = M[i];

    free(work);
    free(M);

    gettimeofday(&end,NULL);
    int diff = 1e6*(end.tv_sec-start.tv_sec) + end.tv_usec-start.tv_usec;

    printf("Inverted fermion matrix in %.3g seconds\n", 1e-6*diff);
    Dinv_file = fopen(filename,"wb");
    if (Dinv_file){
      fwrite(Dinv, VOLUME*VOLUME/2, sizeof(double), Dinv_file);
      fclose(Dinv_file);
    }
  }
  free(ipiv);
}

#endif












