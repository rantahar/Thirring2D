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

#ifdef FLUCTUATION_MATRIX

/* storage */
extern int    ***eta;   //Staggered eta matrix
extern double *Dinv;             //Storage for even to odd propagator
extern int    *evenlist, *oddlist;  //Lists of occupied sites
int    *unoccupied_evenlist, *unoccupied_oddlist;  //Lists of occupied sites
int *added_evensites, *added_oddsites, *removed_evenlist, *removed_oddlist;
int n_added_even=0,n_added_odd=0,n_removed_even=0,n_removed_odd=0;
extern double m;
extern double U;
extern double mu;

void update_linklists();

/* Maximum number of fluctuations from the background configuration */
extern int max_changes;

/* Monomers and links
 * field stores both, 0 for empty, 1 for monomer and 2+dir for links
 */
extern int n_monomers;
extern int n_links;
extern int **field;
int n_bc_monomers;  //Numbers of links and monomers in the background
int n_bc_links;

/* Inverse of G for the backgroud config
 */
extern double *Ginv;


/* Neighbour index arrays, to be filled at the beginning
 */
extern int *tup,*xup,*tdn,*xdn;

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





/* Matrices needed for the fluctuation matrix */
double **BAPe, **PACe, **BGCe, **BAPo, **PACo, **BGCo;
int *newsite, **new_combination;

/* Invert the matrix of propagators between occupied sites
 * Assign the current configuration as the background
 */
extern int current_sign;
double accepted_det = 1;
double fluctuation_det = 1;
int bc_sign = 1;
int fluctuation_sign=1;
#ifdef DEBUG
double bc_determinant = 1;
double current_det = 1;
#endif

void update_background( )
{

 static int init = 1;
 static double *B,*C;
 if(init==1){
   BAPe = malloc(VOLUME/2*sizeof(double*));
   PACe = malloc(VOLUME/2*sizeof(double*));
   BGCe = malloc(VOLUME/2*sizeof(double*));
   BAPo = malloc(VOLUME/2*sizeof(double*));
   PACo = malloc(VOLUME/2*sizeof(double*));
   BGCo = malloc(VOLUME/2*sizeof(double*));
   new_combination = malloc( sizeof(int*)*VOLUME/2 );
   for(int i=0; i<VOLUME/2; i++){
     BAPe[i] = malloc(VOLUME/2*sizeof(double));
     PACe[i] = malloc(VOLUME/2*sizeof(double));
     BGCe[i] = malloc(VOLUME/2*sizeof(double));
     BAPo[i] = malloc(VOLUME/2*sizeof(double));
     PACo[i] = malloc(VOLUME/2*sizeof(double));
     BGCo[i] = malloc(VOLUME/2*sizeof(double));
     new_combination[i] = malloc( sizeof(int)*VOLUME/2 );
   }
   
   B = malloc( sizeof(double)*VOLUME*VOLUME/4 );
   C = malloc( sizeof(double)*VOLUME*VOLUME/4 );
   newsite = malloc( sizeof(int)*VOLUME );
 }
 
 /* The current configuration becomes the background */
 n_bc_monomers = n_monomers;
 n_bc_links = n_links;

 double det_e=1, det_o=1; 
 int n = n_monomers/2 + n_links;
 int m = VOLUME/2 - n;
 double *Ginv_odd = Ginv+n*n;
 double *Dinv_odd = Dinv+VOLUME*VOLUME/4;

#ifdef DEBUG
 struct timeval start, end;
 gettimeofday(&start,NULL);
#endif

 /*If no occupied sites, determinant of rank 0 is 1 */
 if( n > 0 ){
  int *ipiv,*ipiv_o;
  int info;
  
  ipiv = malloc( n*sizeof(int) ); ipiv_o = malloc( n*sizeof(int) );
  
  /* Find occupied and unoccupied sites, construct lists */
  int i=0,j=0,a=0,b=0;
  for (int t=0; t<NT; t++) for (int x=0; x<NX; x++){
    if(field[t][x]>0) {
      if((t+x)%2==0) {
        evenlist[i] = (NX*t + x)/2;
        i++;
      } else {
        oddlist[j] = (NX*t + x)/2;
        j++;
      }
    } else {
      if((t+x)%2==0) {
        unoccupied_evenlist[a] = (NX*t + x)/2;
        a++;
      } else {
        unoccupied_oddlist[b] = (NX*t + x)/2;
        b++;
      }
    }
  }
  if(i!=j || j!=n){
    printf("Number on occupied sites doesn't match, i=%d, j=%d, n=%d, n_links=%d, n_monomers=%d\n",i,j,n,n_links,n_monomers);
    //print_config();
    exit(1);
  }

  /* Construct the inverse of the occupied to occupied propagator */
  /* get propagator */
  for(int i=0; i<n; i++){
    for(int j=0; j<n; j++){
      Ginv[i*n+j] = Dinv[evenlist[i]*VOLUME/2+oddlist[j]];
      Ginv_odd[i*n+j] = Dinv_odd[oddlist[i]*VOLUME/2+evenlist[j]];
    }
  }

  /* LU */
  LAPACK_dgetrf( &n, &n, Ginv, &n, ipiv, &info );
  if( info != 0 ) {
    printf("sgetrf returned error %d (zero determinant has been accepted)! \n", info);
#ifdef DEBUG
    printf(" Incorrect determinant, accepted det %g , accepted factor %g \n", current_det, accepted_det);
#endif
    exit(-1);
  }
  LAPACK_dgetrf( &n, &n, Ginv_odd, &n, ipiv_o, &info );
  if( info != 0 ) {
    printf("sgetrf returned error %d (zero determinant has been accepted)! \n", info);
#ifdef DEBUG
    printf(" Incorrect determinant, accepted det %g , accepted factor %g \n", current_det, accepted_det);
#endif
    exit(-1);
  }

#ifdef DEBUG
  /* determinant from LU ( for debugging ) */
  for(int i=0; i<n; i++) {
    if(ipiv[i]==i+1) { det_e *= Ginv[i*n+i];}
    else { det_e *= -Ginv[i*n+i]; }
    if(ipiv_o[i]==i+1) { det_o *= -Ginv_odd[i*n+i];}
    else { det_o *= Ginv_odd[i*n+i]; }
  }
#endif

  /* The inversion */
  int lwork=n*n;
  double *work;
  work = malloc( lwork*sizeof(double) );
  if (NULL == work) {
    printf("failed to allocate work matrix in update_background (size n=%d) \n",n);
    exit(-1);
  }
  
  LAPACK_dgetri(&n, Ginv, &n, ipiv, work, &lwork, &info);
  if( info != 0 ) {
    printf("sgetri returned error %d! \n", info);
    //exit(-1);
  }
  LAPACK_dgetri(&n, Ginv_odd, &n, ipiv_o, work, &lwork, &info);
  if( info != 0 ) {
    printf("sgetri returned error %d! \n", info);
    //exit(-1);
  }
  
  free(work);
  free(ipiv); free(ipiv_o);
 }


#ifdef DEBUG
  current_det = bc_determinant*accepted_det;
  double det_diff =  fabs(det_e*det_o - current_det);
  printf("EXACT %g %g det %g  accepted %g  diff %g, accepted factor %g \n", det_e, det_o, det_e*det_o, current_det, det_diff, accepted_det);
  if(init == 0 && det_diff/fabs(det_e*det_o)>0.0000001){
    printf(" Incorrect determinant, old bc det %g, fluctuation %g \n", bc_determinant,accepted_det);
    exit(1);
  }
  if(current_sign*det_e*det_o < 0){
    printf("Wrong sign\n");
    exit(1);
  }

  bc_determinant = det_e*det_o;
  current_sign = (det_e*det_o)<0 ? -1 : 1;
#endif

  accepted_det = fluctuation_det = 1;
  fluctuation_sign = 1;
  bc_sign = current_sign ;

#ifdef DEBUG
  gettimeofday(&end,NULL);
  double diff = 1e6*(end.tv_sec-start.tv_sec) + end.tv_usec-start.tv_usec;
  printf("Inverted propagator matrices in %.3g seconds\n", 1e-6*diff);
#endif 


  for( int a=0; a<VOLUME; a++ ) newsite[a] = 1;
  for( int a=0; a<VOLUME/2; a++ ) for( int b=0; b<VOLUME/2; b++ ) 
    new_combination[a][b] = 1;

  n_added_even = n_added_odd = 0;
  n_removed_even = n_removed_odd = 0;

  if(init==1) init = 0;
}



void new_link(int t, int x, int nu){
  link_on(t,x,nu);
#ifdef DEBUG
  current_det = bc_determinant*fluctuation_det;
#endif
  accepted_det = fluctuation_det;
  current_sign = bc_sign*fluctuation_sign;

  int t2 = tdir(t,nu), x2=xdir(x,nu);
  int s1 = (t*NX + x)/2, s2 = (t2*NX+x2)/2;
  
  if( (t+x)%2 == 1 ) {
    int s=s2; s2=s1; s1=s; 
  }
  
  /* Check for the site in added monomers */ 
  int new_site = 1;
  for( int k=0; k<n_removed_even; k++ ) if( evenlist[removed_evenlist[k]] == s1 ) {
    n_removed_even--;
    new_site = 0;
  }
  if( new_site ){
    n_added_even++; /* Already added to lists in det_added_link() */
  }

  new_site = 1;
  for( int k=0; k<n_removed_odd; k++ ) if( oddlist[removed_oddlist[k]] == s2 ) {
    n_removed_odd--;
    new_site = 0;
  }
  if( new_site ){
    n_added_odd++; /* Already added to lists in det_added_link() */
  }
  
  update_linklists();
  if( n_added_even == max_changes || n_added_odd == max_changes  ){
    update_background();
    n_removed_even=n_removed_odd=0;
    n_added_even=n_added_odd=0;
  }
}
void new_monomer(int t1, int x1, int t2, int x2){
  n_monomers += 2;
  monomer_on(t1,x1,t2,x2);

#ifdef DEBUG
  current_det = bc_determinant*fluctuation_det;
#endif
  accepted_det = fluctuation_det;
  current_sign = bc_sign*fluctuation_sign;

  int s1 = (t1*NX + x1)/2, s2 = (t2*NX+x2)/2;
  
  if( (t1+x1)%2 == 1 ) {
    int s=s2; s2=s1; s1=s; 
  }
  
  /* Check for the site in added monomers */ 
  int new_site = 1;
  for( int k=0; k<n_removed_even; k++ ) if( evenlist[removed_evenlist[k]] == s1 ) {
    n_removed_even--;
    new_site = 0;
  }
  if( new_site ){
    n_added_even++; /* Already added to lists in det_added_link() */
  }

  new_site = 1;
  for( int k=0; k<n_removed_odd; k++ ) if( oddlist[removed_oddlist[k]] == s2 ) {
    n_removed_odd--;
    new_site = 0;
  }
  if( new_site ){
    n_added_odd++; /* Already added to lists in det_added_link() */
  }

  update_linklists();
  if( n_added_even == max_changes || n_added_odd == max_changes  ){
    update_background();
    n_removed_even=n_removed_odd=0;
    n_added_even=n_added_odd=0;
  }
}
void removed_link(int t, int x, int nu){
  n_links -= 1;
  link_off(t,x,nu);

#ifdef DEBUG
  current_det = bc_determinant*fluctuation_det;
#endif
  accepted_det = fluctuation_det;
  current_sign = bc_sign*fluctuation_sign;

  int t2 = tdir(t,nu), x2=xdir(x,nu);
  int s1 = (t*NX + x)/2, s2 = (t2*NX+x2)/2;
  
  if( (t+x)%2 == 1 ) {
    int s=s2; s2=s1; s1=s; 
  }
  
  /* Check for the site in added monomers */ 
  int new_site = 1;
  for( int i=0; i<n_added_even; i++ ) if( added_evensites[i] == s1 ) {
    n_added_even--;
    new_site = 0;
  }
  if( new_site ){
    n_removed_even++; /* Already added to lists in det_added_link() */
  }

  new_site = 1;
  for( int i=0; i<n_added_odd; i++ ) if( added_oddsites[i] == s2 ) {
    n_added_odd--;
    new_site = 0;
  }
  if( new_site ){
    n_removed_odd++; /* Already added to lists in det_added_link() */
  }
 
#ifdef DEBUG
  printf("Removed link at (%d,%d,%d) %g \n",t,x,nu,current_det);
  printf(" indexes %d %d %d %d \n",n_added_even,n_added_odd,n_removed_even,n_removed_odd);
#endif
  update_linklists();
  if( n_removed_even == max_changes || n_removed_odd == max_changes ){
    update_background();
    n_removed_even=n_removed_odd=0;
    n_added_even=n_added_odd=0;
  }
}
void removed_monomer(int t1, int x1, int nu){
  n_monomers -= 2;
  link_off(t1,x1,nu);
  int t2 = tdir(t1,nu), x2=xdir(x1,nu);

#ifdef DEBUG
  current_det = bc_determinant*fluctuation_det;
#endif
  accepted_det = fluctuation_det;
  current_sign = bc_sign*fluctuation_sign;

  if( (t1+x1)%2 == 1 ) { int t=t2, x=x2; t2=t1;x2=x1; t1=t;x1=x;}

  /* Check for the site in added monomers */ 
  int s = (t1*NX + x1)/2;
  int new_site = 1;
  for( int i=0; i<n_added_even; i++ ) if( added_evensites[i] == s ) {
    n_added_even--;
    new_site = 0;
  } 
  if( new_site ){
    n_removed_even++; /* Already added to lists in det_added_link() */
  }

  s = (t2*NX + x2)/2; 
  new_site = 1;
  for( int i=0; i<n_added_odd; i++ ) if( added_oddsites[i] == s ) {
    n_added_odd--;
    new_site = 0;
  } 
  if( new_site ){
    n_removed_odd++; /* Already added to lists in det_added_link() */
  }

#ifdef DEBUG
  printf("Removed monomer at (%d,%d,%d) %g \n",t1,x1,nu,current_det);
  printf(" indexes %d %d %d %d \n",n_added_even,n_added_odd,n_removed_even,n_removed_odd);
#endif
  update_linklists();
  if(n_removed_even == max_changes || n_removed_odd == max_changes ){
    update_background();
    n_removed_even=n_removed_odd=0;
    n_added_even=n_added_odd=0;
  }
}


int moved_new_site=1;
int moved_old_site=1;
void moved_source_monomer(int t2, int x2, int t4, int x4){
  int s = (t2*NT+x2)/2;
  if( moved_old_site ) {
    if( (t2+x2)%2==0 ) { /* The lists are already ordered properly, so just adjust lengths */
      n_added_even++;
      n_removed_even++;
    } else {
      n_removed_odd++;
      n_added_odd++;
    }
  }
  if( !moved_new_site ) {   // Found in the list of removed sites
    if( (t2+x2)%2==0 ) {
      n_added_even--;
      n_removed_even--;
    } else {
      n_removed_odd--;
      n_added_odd--;
    }
  }

  field[t2][x2] = 0; field[t4][x4] = SOURCE_MONOMER;
#ifdef DEBUG
  current_det = bc_determinant*fluctuation_det;
#endif
  accepted_det = fluctuation_det;
  current_sign = bc_sign*fluctuation_sign;

  #ifdef DEBUG
  print_config();
  check_det();
  #endif
  if(n_added_even == max_changes || n_added_odd == max_changes ||
     n_removed_even == max_changes || n_removed_odd == max_changes ){
  //if(n_added_even == 2 || n_added_odd == 2 ||
  //  n_removed_even == 2 || n_removed_odd == 2 ){
    update_background();
    n_removed_even=n_removed_odd=0;
    n_added_even=n_added_odd=0;
  }
}

void not_moved_source_monomer(int t2, int x2, int t4, int x4){
  int n = n_bc_monomers/2 + n_bc_links;
  int s = (t4*NX+x4)/2;
  if(!moved_old_site){
    if( (t2+x2)%2==0 ) {
      added_evensites[n_added_even-1] = (t2*NX+x2)/2;
    } else {
      added_oddsites[n_added_odd-1] = (t2*NX+x2)/2;
    }
  }
  if( !moved_new_site ) {
    if( (t2+x2)%2==0 ) {
      for( int k=0; k<n; k++ ) if( evenlist[k] == s )
        removed_evenlist[n_removed_even-1] = k;
    } else {
      for( int k=0; k<n; k++ ) if( oddlist[k] == s )
        removed_oddlist[n_removed_odd-1] = k;
    }
  }
}



/* Determinants after adding or removing links */
double extended_determinant( int me, int mo, int re, int ro ){
  int n = n_bc_monomers/2 + n_bc_links;
  int mr = me + ro;

#ifdef DEBUG
  printf("extended_determinant, me %d mo %d re %d ro %d\n",me,mo,re,ro);
  printf("extended_determinant, n_added_even %d n_added_odd %d n_removed_even %d n_removed_odd %d, n %d\n",n_added_even,n_added_odd,n_removed_even,n_removed_odd,n);
#endif

  if(mr==0){
    fluctuation_sign = 1;
    fluctuation_det = 1;
    return 1;
  }

  /* Construct new parts of the fluctuation matrix */
  double *Dinv_odd = Dinv+VOLUME*VOLUME/4;
  double *Ginv_odd = Ginv+n*n;
  for( int b=0; b<me; b++) if( newsite[added_evensites[b]] == 1 ){
    for( int a=0; a<n; a++)
    { 
      double sum = 0;
      for( int l=0; l<n; l++) sum += Ginv_odd[a*n+l]*Dinv_odd[oddlist[l]*VOLUME/2+added_evensites[b]];
      PACo[added_evensites[b]][a] = -sum;
      
      sum=0;
      for( int l=0; l<n; l++) sum += Dinv[added_evensites[b]*VOLUME/2+oddlist[l]]*Ginv[l*n+a];
      BAPe[a][added_evensites[b]] = sum;
    }
    newsite[added_evensites[b]] = 0;
  }

  for( int a=0; a<mo; a++) if( newsite[VOLUME/2+added_oddsites[a]] == 1 ){
    for( int b=0; b<n; b++)  { 
      double sum = 0;
      for( int k=0; k<n; k++) sum += Ginv[b*n+k]*Dinv[evenlist[k]*VOLUME/2+added_oddsites[a]];
      PACe[added_oddsites[a]][b] = -sum;

      sum = 0;
      for( int l=0; l<n; l++)  sum += Dinv_odd[added_oddsites[a]*VOLUME/2+evenlist[l]]*Ginv_odd[l*n+b];
      BAPo[b][added_oddsites[a]] = sum;
    }
    newsite[VOLUME/2+added_oddsites[a]] = 0;
  }


  for( int b=0; b<me; b++) for( int a=0; a<mo; a++) 
  if( new_combination[added_evensites[b]][added_oddsites[a]] == 1 )
  {
    double sum = 0;
    for( int l=0; l<n; l++)
      sum += Dinv[added_evensites[b]*VOLUME/2+oddlist[l]]*PACe[added_oddsites[a]][l];
    BGCe[added_evensites[b]][added_oddsites[a]] = Dinv[added_evensites[b]*VOLUME/2+added_oddsites[a]] + sum;

    sum=0;
    for( int l=0; l<n; l++) 
      sum += Dinv_odd[added_oddsites[a]*VOLUME/2+evenlist[l]]*PACo[added_evensites[b]][l];
    BGCo[added_oddsites[a]][added_evensites[b]] = Dinv_odd[added_oddsites[a]*VOLUME/2+added_evensites[b]] + sum;

    new_combination[added_evensites[b]][added_oddsites[a]] = 0;
  }


  /* Even to odd */
  double *F;
  F = malloc(mr*mr*sizeof(double));
  for( int b=0; b<ro; b++) for( int a=0; a<re; a++)
    F[a*mr+b] = Ginv[removed_oddlist[b]*n+removed_evenlist[a]];
  for( int a=0; a<mo; a++) for( int b=0; b<ro; b++)
    F[(re+a)*mr+b] = PACe[added_oddsites[a]][removed_oddlist[b]];
  for( int a=0; a<re; a++) for( int b=0; b<me; b++)
    F[a*mr+(ro+b)] = BAPe[removed_evenlist[a]][added_evensites[b]];
  for( int b=0; b<me; b++) for( int a=0; a<mo; a++)
    F[(re+a)*mr+(ro+b)] = BGCe[added_evensites[b]][added_oddsites[a]];


  /* determinant from LU */
  double det=1;
  int ipiv[mr], info;
  LAPACK_dgetrf( &mr, &mr, F, &mr, ipiv, &info );
  for( int a=0; a<mr; a++) {
    det *= F[a*mr+a];
    if(ipiv[a]!=a+1) det *= -1;
  }

  /* Odd to even */
  for( int b=0; b<re; b++) for( int a=0; a<ro; a++)
    F[a*mr+b] = Ginv[n*n+removed_evenlist[b]*n+removed_oddlist[a]];
  for( int a=0; a<me; a++) for( int b=0; b<re; b++)
    F[(a+ro)*mr+b] = PACo[added_evensites[a]][removed_evenlist[b]];
  for( int a=0; a<ro; a++) for( int b=0; b<mo; b++)
    F[a*mr+(b+re)] = BAPo[removed_oddlist[a]][added_oddsites[b]];
  for( int b=0; b<mo; b++) for( int a=0; a<me; a++)
    F[(ro+a)*mr+(re+b)] = BGCo[added_oddsites[b]][added_evensites[a]];

  /* determinant from LU */
  LAPACK_dgetrf( &mr, &mr, F, &mr, ipiv, &info );
  for( int a=0; a<mr; a++) {
    det *= -F[a*mr+a];
    if(ipiv[a]!=a+1) det *= -1;
  }

  free(F);

  fluctuation_sign = det<0 ? -1: 1;
  fluctuation_det = det;

  return( det );
}





double det_added_link(int t, int x, int t2, int x2){

  int n = n_bc_monomers/2 + n_bc_links;
  int me=n_added_even, mo=n_added_odd, re=n_removed_even, ro=n_removed_odd;

  /* Make sure that (t,x) is even and (t2,x2) is odd */
  if( (t+x)%2 == 1 ) { int t1=t2, x1=x2; t2=t;x2=x; t=t1;x=x1; }

  int new_even_site = 1;// even_site_in_removed_list((NX*t + x)/2, n_removed_even);
  for( int k=0; k<n_removed_even; k++ ) if( evenlist[removed_evenlist[k]] == (NX*t + x)/2 ) {
    re--;
    int tmp = removed_evenlist[k];
    for( int l=k; l<n_removed_even-1; l++) removed_evenlist[l] = removed_evenlist[l+1];
    removed_evenlist[n_removed_even-1] = tmp;
    new_even_site = 0;
    break;
  } 
  if( new_even_site ) {
    me++;
    added_evensites[n_added_even] = (NX*t + x)/2 ;
  }

  int new_odd_site = 1;
  for( int k=0; k<n_removed_odd; k++ ) if( oddlist[removed_oddlist[k]] == (NX*t2 + x2)/2 ) {
    ro--;
    int tmp = removed_oddlist[k];
    for( int l=k; l<n_removed_odd-1; l++) removed_oddlist[l] = removed_oddlist[l+1];
    removed_oddlist[n_removed_odd-1] = tmp;
    new_odd_site = 0;
    break;
  } 
  if( new_odd_site ){
    mo++;
    added_oddsites[n_added_odd] = (NX*t2 + x2)/2;
  }

  double det = extended_determinant(me,mo,re,ro);
  double detratio = fabs(det/accepted_det);
#ifdef DEBUG
  printf("Adding at (%d,%d) (%d) and (%d,%d) (%d)\n",t,x,(t*NX+t)/2,t2,x2,(t2*NX+x2)/2);
  printf(" new det %g  %g %g\n",current_det*detratio, det, accepted_det );
#endif

  return( detratio );
}






double det_removed_link(int t, int x, int t2, int x2 ){
  int n = n_bc_monomers/2 + n_bc_links;

  //Make sure that (t,x) is even and (t2,x2) is odd
  if( (t+x)%2 == 1 ) { int t1=t2, x1=x2; t2=t;x2=x; t=t1;x=x1; }

  // Check for the removed links in added, if found, reorder, otherwise ad to removed list
  int me=n_added_even, mo=n_added_odd, re=n_removed_even, ro=n_removed_odd;

  int new_even_site=1;
  for( int k=0; k<n_added_even; k++ ) if( added_evensites[k] == (NX*t+x)/2 ) {
    me--;
    for( int l=k; l<n_added_even-1; l++) added_evensites[l] = added_evensites[l+1];
    added_evensites[n_added_even-1] = (NX*t+x)/2;
    new_even_site = 0;
    break;
  } 
  if( new_even_site ){
    re++;
    for( int k=0; k<n; k++) if( evenlist[k] == (NX*t+x)/2 )
      removed_evenlist[n_removed_even] = k;
  }

  int new_odd_site=1;
  for( int k=0; k<n_added_odd; k++ ) if( added_oddsites[k] == (NX*t2+x2)/2 ) {
    mo--;
    for( int l=k; l<n_added_odd-1; l++) added_oddsites[l] = added_oddsites[l+1];
    added_oddsites[n_added_odd-1] = (NX*t2+x2)/2;
    new_odd_site = 0;
    break;
  } 
  if( new_odd_site ){
    ro++;
    for( int k=0; k<n; k++) if( oddlist[k] == (NX*t2+x2)/2 )
      removed_oddlist[n_removed_odd] = k;
  }

  double det = extended_determinant(me,mo,re,ro);
  double detratio = fabs(det/accepted_det);
#ifdef DEBUG
  printf("Removing at (%d,%d) and (%d,%d), values %d and %d\n",t,x,t2,x2,field[t][x],field[t2][x2]);
  printf(" new det %g  %g %g\n",current_det*detratio, det, accepted_det );
#endif

  return( detratio );
}



double det_moved_monomer(int t, int x, int t2, int x2){
  int n = n_bc_monomers/2 + n_bc_links;
  double detratio;

  // Check for the chosen sites in added and removed lists */
  moved_old_site = 1;
  int s = (t*NX+x)/2;
  if( (t+x)%2==0 ) for( int k=0; k<n_added_even; k++ ) if( added_evensites[k] == s )
    moved_old_site = 0;
  if( (t+x)%2==1 ) for( int k=0; k<n_added_odd; k++ ) if( added_oddsites[k] == s )
    moved_old_site = 0;

  moved_new_site = 1;
  s = (t2*NX+x2)/2;
  if( (t+x)%2==0 ) for( int k=0; k<n_removed_even; k++ ) if( evenlist[removed_evenlist[k]] == s )
    moved_new_site = 0;
  if( (t+x)%2==1 ) for( int k=0; k<n_removed_odd; k++ ) if( oddlist[removed_oddlist[k]] == s )
    moved_new_site = 0;


  /* Make sure that (t,x) is even and (t2,x2) is odd */
  if( (t+x)%2 != (t2+x2)%2) { 
    //printf("Attempting to move to opposite parity\n");
    //exit(1);
    return 0;
  }
  if( (t+x)%2 == 0 ) {

    int me = n_added_even, re = n_removed_even;

    /* Remove if found in list, otherwise add a to list of removed */
    int new_even_site = 1;
    for( int k=0; k<n_removed_even; k++ ) if( evenlist[removed_evenlist[k]] == (NX*t2 + x2)/2 ) {
      int tmp = removed_evenlist[k];
      for( int l=k; l<n_removed_even-1; l++) removed_evenlist[l] = removed_evenlist[l+1];
      removed_evenlist[n_removed_even-1] = tmp;
      new_even_site = 0;
      re--;
      break;
    } 
    if( new_even_site ) {
      me++;
      added_evensites[n_added_even] = (NX*t2 + x2)/2 ;
    }

    new_even_site=1;
    for( int k=0; k<me; k++ ) if( added_evensites[k] == (NX*t+x)/2 ) {
      for( int l=k; l<me-1; l++) added_evensites[l] = added_evensites[l+1];
      added_evensites[me-1] = (NX*t+x)/2;
      new_even_site = 0;
      me--;
      break;
    } 
    if( new_even_site ){
      for( int k=0; k<n; k++) if( evenlist[k] == (NX*t+x)/2 )
        removed_evenlist[re] = k;
      re++;
    }

    double det = extended_determinant(me,n_added_odd,re,n_removed_odd);
    detratio = fabs(det/accepted_det);
#ifdef DEBUG
    printf("Moving (%d,%d) (%d) to (%d,%d) (%d)\n",t,x,(NX*t+x)/2,t2,x2,(NX*t2+x2)/2);
    printf(" new det %g  %g %g\n",current_det*detratio, det, accepted_det );
#endif
    
  } else {

    int mo = n_added_odd, ro = n_removed_odd;

    int new_odd_site = 1;
    for( int k=0; k<n_removed_odd; k++ ) if( oddlist[removed_oddlist[k]] == (NX*t2 + x2)/2 ) {
      int tmp = removed_oddlist[k];
      for( int l=k; l<n_removed_odd-1; l++) removed_oddlist[l] = removed_oddlist[l+1];
      removed_oddlist[n_removed_odd-1] = tmp;
      new_odd_site = 0;
      ro--;
      break;
    } 
    if( new_odd_site ) {
      mo++;
      added_oddsites[n_added_odd] = (NX*t2 + x2)/2 ;
    }

    new_odd_site=1;
    for( int k=0; k<mo; k++ ) if( added_oddsites[k] == (NX*t+x)/2 ) {
      for( int l=k; l<mo-1; l++) added_oddsites[l] = added_oddsites[l+1];
      added_oddsites[mo-1] = (NX*t+x)/2;
      new_odd_site = 0;
      mo--;
      break;
    } 
    if( new_odd_site ){
      for( int k=0; k<n; k++) if( oddlist[k] == (NX*t+x)/2 )
        removed_oddlist[ro] = k;
      ro++;
    }


    double det = extended_determinant(n_added_even,mo,n_removed_even,ro);
    detratio = fabs(det/accepted_det);
#ifdef DEBUG
    printf("Moving (%d,%d) (%d) to (%d,%d) (%d)\n",t,x,(NX*t+x)/2,t2,x2,(NX*t2+x2)/2);
    printf(" new det %g  %g %g\n",current_det*detratio, det, accepted_det );
#endif
  }
  return( detratio );
}





#ifdef DEBUG
/* Check the determinant
 */
void check_det(  )
{
 double det_e=1,det_o=1; 
 int n = n_monomers/2 + n_links;
 int * check_oddlist, * check_evenlist;
 double * check_Ginv, * check_Ginv_odd;

 check_Ginv = malloc( VOLUME*VOLUME*sizeof(double)/4 );
 check_Ginv_odd = malloc( VOLUME*VOLUME*sizeof(double)/4 );
 check_evenlist = malloc( VOLUME*sizeof(int)/2 );
 check_oddlist  = malloc( VOLUME*sizeof(int)/2 );

 /*If no occupied sites, determinant of rank 0 is 1 */
 if( n > 0 ){
  int ipiv[n],ipiv_o[n];
  int info;
  
  /* Find occupied and unoccupied sites, construct lists */
  int i=0,j=0,a=0,b=0;
  for (int t=0; t<NT; t++) for (int x=0; x<NX; x++){
    if(field[t][x]>0) {
      if((t+x)%2==0) {
        check_evenlist[i] = (NX*t + x)/2;
        i++;
      } else {
        check_oddlist[j] = (NX*t + x)/2;
        j++;
      }
    }
  }

  if(i!=j || j!=n){
    printf("Number of occupied sites doesn't match, i=%d, j=%d, n=%d, n_links=%d, n_monomers=%d\n",i,j,n,n_links,n_monomers);
    //print_config();
    exit(1);
  }

  /* Construct the inverse of the occupied to occupied propagator */
  /* get propagator */
  for(int i=0; i<n; i++){
    for(int j=0; j<n; j++){
      check_Ginv[i*n+j] = Dinv[check_evenlist[i]*VOLUME/2+check_oddlist[j]];
      check_Ginv_odd[i*n+j] = Dinv[VOLUME/2*VOLUME/2 + check_oddlist[i]*VOLUME/2+check_evenlist[j]];
    }
  }

  /* LU */
  LAPACK_dgetrf( &n, &n, check_Ginv, &n, ipiv, &info );
  if( info != 0 ) {
    printf("sgetrf returned error %d (zero determinant has been accepted)! \n", info);
#ifdef DEBUG
    printf(" Incorrect determinant, accepted det %g , accepted factor %g \n", current_det, accepted_det); 
#endif
    exit(-1);
  }
  LAPACK_dgetrf( &n, &n, check_Ginv_odd, &n, ipiv_o, &info );
  if( info != 0 ) {
    printf("sgetrf returned error %d (zero determinant has been accepted)! \n", info);
#ifdef DEBUG
    printf(" Incorrect determinant, accepted det %g , accepted factor %g \n", current_det, accepted_det); 
#endif
    exit(-1);
  }

  /* determinant from LU ( for debugging ) */
  for(int i=0; i<n; i++) {
    if(ipiv[i]==i+1) { det_e *= check_Ginv[i*n+i];}
    else { det_e *= -check_Ginv[i*n+i]; }
    if(ipiv_o[i]==i+1) { det_o *= -check_Ginv_odd[i*n+i];}
    else { det_o *= check_Ginv_odd[i*n+i]; }
  }
  
 }

  double current_det = bc_determinant*accepted_det;
  double det_diff =  fabs(det_e*det_o - current_det);
  printf("CHECK %g %g det %g  accepted %g  diff %g, accepted factor %g \n", det_e, det_o, det_e*det_o, current_det, det_diff, accepted_det);
  if( det_diff/fabs(det_e*det_o)>0.0000001){
    printf(" Incorrect determinant, old bc det %g, fluctuation %g \n", bc_determinant,accepted_det);
    exit(1);
  }
  if(current_sign*det_e*det_o < 0){
    printf("Wrong sign\n");
    exit(1);
  }

  free(check_Ginv);
  free(check_Ginv_odd);
  free(check_evenlist);
  free(check_oddlist);
}

/* LAPACK LU loses to cg when L~200 */
void calc_Dinv_cg( )
{
  int n=VOLUME/2;
  double ** source = alloc_vector();
  double ** propagator = alloc_vector();

  struct timeval start, end;
  gettimeofday(&start,NULL);

  /* Construct the full Volume to Volume Dirac matrix
   * Odd to even here, the inverse will be even to odd
   */
  for (int t1=0; t1<NT; t1++) for (int x1=0; x1<NX; x1++) if((t1+x1)%2==1) {
    int i1 = (NX*t1 + x1)/2;
    vec_zero(source);
    source[t1][x1] = 1;
    cg_propagator( propagator, source );
    for (int t2=0; t2<NT; t2++) for (int x2=0; x2<NX; x2++) if((t2+x2)%2==0) {
      int i2 = (NX*t2 + x2)/2;
      Dinv[i2*n + i1] = propagator[t2][x2];
    }
  }
 
  free_vector(source);
  free_vector(propagator);

  gettimeofday(&end,NULL);
  int diff = 1e6*(end.tv_sec-start.tv_sec) + end.tv_usec-start.tv_usec;
  printf("Inverted fermion matrix in %.3g seconds\n", 1e-6*diff);
}


#endif

#endif




