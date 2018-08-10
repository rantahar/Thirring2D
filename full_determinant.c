/****************************************************************
 * Simulate the thirring model using the fermion bag algorithm 
 * (arXiv:0910.5736). The mass term is represented as a field of 
 * monomers (occupied sites) and the four fermion term is  
 * represented as dimers (occupied links). 
 *
 ****************************************************************/
#ifdef DEBUG
#include <fenv.h>
void check_det(  );
#endif

#include "Thirring.h"

#ifndef FLUCTUATION_MATRIX

/* storage */
extern int    ***eta;   //Staggered eta matrix
extern double *Dinv;             //Storage for even to odd propagator
extern int    *evenlist, *oddlist;  //Lists of occupied sites
int    *unoccupied_evenlist, *unoccupied_oddlist;  //Lists of occupied sites
extern double m;
extern double U;
extern double mu;

void update_linklists();

/* Monomers and links
 * field stores both, 0 for empty, 1 for monomer and 2+dir for links
 */
extern int n_monomers;
extern int n_links;
extern int **field;

/* Neighbour index arrays, to be filled at the beginning
 */
extern int *tup,*xup,*tdn,*xdn;

/* Functions for fetching neighboring coordinates */
static inline int tdir(int t, int dir){
  if( dir == TUP ) return tup[t];
  else if ( dir == TDN ) return tdn[t];
  else return(t);
}

static inline int xdir(int x, int dir){
  if( dir == XUP ) return xup[x];
  else if ( dir == XDN ) return xdn[x];
  else return(x);
}

/* Opposite of a direction */
static inline int opp_dir(int dir){
  return ( dir + ND ) % NDIRS;
}

/* Turn a link on */
static inline void link_on(int t, int x, int dir){
  int t2 = tdir(t, dir);
  int x2 = xdir(x, dir);
  if ( field[t][x] == 0 && field[t2][x2] == 0 ){
    field[t][x] = dir+2;
    field[t2][x2] = opp_dir(dir)+2;
#ifdef DEBUG
    printf("Turned on link (%d,%d,%d) (%d,%d,%d)\n",t,x,dir,t2,x2,opp_dir(dir));
#endif
  } else {
    printf("Link already occupied\n");
    exit(1);
  }
}

/* Turn a monomer on at a link */
static inline void monomer_on(int t, int x, int t2, int x2){
  if ( field[t][x] == 0 && field[t2][x2] == 0 ){
    field[t][x] = MONOMER;
    field[t2][x2] = MONOMER;
#ifdef DEBUG
    printf("Turned on link (%d,%d) (%d,%d)\n",t,x,t2,x2);
#endif
  } else {
    printf("Link already occupied\n");
    exit(1);
  }
}

/* Turn a link off */
static inline void link_off(int t, int x, int dir){
  int t2 = tdir(t, dir);
  int x2 = xdir(x, dir);
  if ( field[t][x] != 0 && field[t2][x2] != 0 ){
    field[t][x] = 0;
    field[t2][x2] = 0;
#ifdef DEBUG
    printf("Turned off link (%d,%d) (%d,%d)\n",t,x,t2,x2);
#endif
  } else {
    printf("Link already off\n");
    exit(1);
  }
}


/* Check if it's legal to add a link or monomer */
static inline int is_legal(int t, int x, int nu){
  int t2,x2;
  t2 = tdir(t,nu); x2 = xdir(x,nu);
  return ( field[t][x] == 0 && field[t2][x2] == 0 );
}





/* Invert the matrix of propagators between occupied sites
 * Assing the current configuration as the background
 */
extern int current_sign;
double accepted_det = 1;
double previous_det=1;

void new_link(int t, int x, int nu){
  link_on(t,x,nu);
  accepted_det = previous_det;
  current_sign = accepted_det>0? 1:-1;
}

void new_monomer(int t1, int x1, int t2, int x2){
  n_monomers += 2;
  monomer_on(t1,x1,t2,x2);
  accepted_det = previous_det;
  current_sign = accepted_det>0? 1:-1;
}

void removed_link(int t, int x, int nu){
  n_links -= 1;
  link_off(t,x,nu);
  accepted_det = previous_det;
  current_sign = accepted_det>0? 1:-1;
}

void removed_monomer(int t1, int x1, int nu){
  n_monomers -= 2;
  link_off(t1,x1,nu);
  accepted_det = previous_det;
  current_sign = accepted_det>0? 1:-1;
}

void moved_source_monomer(int t2, int x2, int t4, int x4){
  field[t2][x2] = 0; field[t4][x4] = SOURCE_MONOMER;
  accepted_det = previous_det;
  current_sign = accepted_det>0? 1:-1;
}

void not_moved_source_monomer(int t2, int x2, int t4, int x4){
}


//NOTE: Does not work with EVENODD
/* Determinants after adding or removing links */
double determinant( ){
  int *ipiv;
  int info;
  double *M;
  double det=1;

  struct timeval start, end;
  gettimeofday(&start,NULL);
 
  ipiv = malloc( VOLUME*sizeof(int) );
  M = malloc( VOLUME*VOLUME*sizeof(double) );


  int i=0,n=0;
  //printf(" \n ");
  for (int t1=0; t1<NT; t1++) for (int x1=0; x1<NX; x1++) if(field[t1][x1]==0) {
    for (int t2=0; t2<NT; t2++) for (int x2=0; x2<NX; x2++) if(field[t2][x2]==0) {
      M[i] = fM_index( t1, x1, t2, x2 );
      //printf(" %.4f ",M[i]);
      i++;
    }
    //printf(" \n ");
    n++;
  }

  if(n>0){
    /* determinant from LU */
    LAPACK_dgetrf( &n, &n, M, &n, ipiv, &info );
    for( int a=0; a<n; a++) {
      det *= M[a*n+a];
      if(ipiv[a]!=a+1) det *= -1;
    }

    free(ipiv);
    free(M);
  }

  gettimeofday(&end,NULL);
  int diff = 1e6*(end.tv_sec-start.tv_sec) + end.tv_usec-start.tv_usec;
  //printf("Calculated determinant in %.3g seconds\n", 1e-6*diff );

  return( det );
}



void init_det(){
  static int init=1;
  if(init==1){
    accepted_det=determinant(); //Do this before the first change
    init=0;
  }
}

double det_added_link(int t, int x, int t2, int x2){
  init_det();
  field[t][x] = 1; field[t2][x2] = 1;

  previous_det = determinant();
  double detratio = fabs(previous_det/accepted_det);

  field[t][x] = 0; field[t2][x2] = 0;
  return( detratio );
}



double det_removed_link(int t, int x, int t2, int x2 ){
  init_det();
  int f1 = field[t][x], f2=field[t2][x2];
  field[t][x] = 0; field[t2][x2] = 0;

  previous_det = determinant();
  double detratio = fabs(previous_det/accepted_det);

  field[t][x] = f1; field[t2][x2] = f2;
  return( detratio );
}



double det_moved_monomer(int t, int x, int t2, int x2){
  init_det();
  int f1 = field[t][x];
  field[t][x] = 0; field[t2][x2] = 1;
  
  previous_det = determinant();
  double detratio = fabs(previous_det/accepted_det);

  field[t][x] = f1; field[t2][x2] = 0;
  return( detratio );
}






#endif




