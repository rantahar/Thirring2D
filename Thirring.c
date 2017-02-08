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

/* storage */
int    ***eta;   //Staggered eta matrix
double *Dinv;             //Storage for even to odd propagator
int    *evenlist, *oddlist;  //Lists of occupied sites
int    *unoccupied_evenlist, *unoccupied_oddlist;  //Lists of occupied sites
int *added_evensites, *added_oddsites, *removed_evenlist, *removed_oddlist;
int n_added_even=0,n_added_odd=0,n_removed_even=0,n_removed_odd=0;
double m;
double U;
double mu;

/* Maximum number of fluctuations from the background configuration */
int max_changes;

/* Monomers and links
 * field stores both, 0 for empty, 1 for monomer and 2+dir for links
 */
int n_monomers=0;
int n_links=0;
int **field;
int n_bc_monomers;  //Numbers of links and monomers in the background
int n_bc_links;

/* Inverse of G for the backgroud config
 */
double *Ginv;


/* Neighbour index arrays, to be filled at the beginning
 */
int *tup,*xup,*tdn,*xdn;

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




/* The fermion matrix
 */
#ifdef ANTISYMMETRIC //Antisymmetrix boundaries
double fM_index( int t1, int x1, int t2, int x2 )
{
  if(t1==t2 ){ 
    if( x2 == xup[x1] || x2 == xdn[x1] ) {
#if NX==2
      if (x2>x1) { return eta[t1][x1][1] ; }
      else { return -eta[t1][x1][1] ; }
#else
      if (x2>x1) { return 0.5*eta[t1][x1][1] ; }
      else { return -0.5*eta[t1][x1][1] ; }
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

#ifdef OPENX //Open in space, antisymmetric in time
double fM_index( int t1, int x1, int t2, int x2 )
{
 if( x1!=NX && x2!=NX ) {
  if(t1==t2 ){ 
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
void calc_Dinv( )
{
  int n=VOLUME/2;
  int *ipiv;
  int info;
  int lwork=n*n;
  double *work;
  double *M;

  ipiv = malloc( n*sizeof(int) );

  FILE * Dinv_file;
  char filename[100];
#ifdef ANTISYMMETRIC
  sprintf(filename, "free_propagator_T%dX%d_mu%0.6f",NT,NX,mu);
#endif
#ifdef OPENX
  sprintf(filename, "free_propagator_T%dX%d_mu%0.6f_open",NT,NX,mu);
#endif
  
  Dinv_file = fopen(filename,"rb");
  if (Dinv_file){
    fread(Dinv, VOLUME*VOLUME/2, sizeof(double), Dinv_file);
    fclose(Dinv_file);
  } else {
  
    struct timeval start, end;
    gettimeofday(&start,NULL);

    /* Construct the full Volume to Volume Dirac matrix
     * Odd to even here, the inverse will be even to odd
     */
    work = malloc( lwork*sizeof(double) );
    M = malloc( lwork*sizeof(double) );
    for (int t1=0; t1<NT; t1++) for (int x1=0; x1<NX; x1++) if((t1+x1)%2==0) {
      int i1 = (NX*t1 + x1)/2;
      for (int t2=0; t2<NT; t2++) for (int x2=0; x2<NX; x2++) if((t2+x2)%2==1) {
        int i2 = (NX*t2 + x2)/2;
        M[i2*n + i1] = fM_index( t1, x1, t2, x2 );
      }
    }

    // LU decompose
    LAPACK_dgetrf( &n, &n, M, &n, ipiv, &info );
    if( info != 0 ) {
      printf("calc_Dinv: sgetrf returned an error %d! \n", info);
      exit(-1);
    }
    double det = 0;
    for(int i=0; i<n; i++) {
      det += log(fabs(M[i*n+i]));
    }

    LAPACK_dgetri(&n, M, &n, ipiv, work, &lwork, &info);
    if( info != 0 ) {
      printf("calc_Dinv: sgetri returned an error %d! \n", info);
      exit(-1);
    }

    for(int i = 0; i < n*n; i++) Dinv[i] = M[i];

    /* Even to odd */
    for (int t1=0; t1<NT; t1++) for (int x1=0; x1<NX; x1++) if((t1+x1)%2==1) {
      int i1 = (NX*t1 + x1)/2;
      for (int t2=0; t2<NT; t2++) for (int x2=0; x2<NX; x2++) if((t2+x2)%2==0) {
        int i2 = (NX*t2 + x2)/2;
        M[i2*n + i1] = fM_index( t1, x1, t2, x2 );
      }
    }
   
    gettimeofday(&start,NULL);
    // LU decompose
    LAPACK_dgetrf( &n, &n, M, &n, ipiv, &info );
    if( info != 0 ) {
      printf("calc_Dinv: sgetrf returned an error %d! \n", info);
      exit(-1);
    }

    for(int i=0; i<n; i++) {
      det += log(fabs(M[i*n+i]));
    } 

    LAPACK_dgetri(&n, M, &n, ipiv, work, &lwork, &info);
    if( info != 0 ) {
      printf("calc_Dinv: sgetri returned an error %d! \n", info);
      exit(-1);
    }
  
    for(int i = 0; i < n*n; i++) Dinv[n*n+i] = M[i];

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

void write_config(){
  FILE * config_file;
  char filename[100];
  sprintf(filename, "config_checkpoint_T%dX%d_U%.6gm%.6gmu%.6g",NT,NX,U,m,mu);

  int * buffer = malloc(NX*NT*sizeof(int));
  for (int t=0; t<NT; t++) for (int x=0; x<NX; x++) 
    buffer[t*NX+x] = field[t][x];
  
  config_file = fopen(filename,"wb");
  if (config_file){
    fwrite(buffer, VOLUME, sizeof(int), config_file);
    fclose(config_file);
  } else {
    printf("Could not write configuration\n");
    exit(1);
  }
  free(buffer);
}


void read_config(char * filename){
  FILE * config_file;
  
  config_file = fopen(filename,"rb");
  int * buffer = malloc(NX*NT*sizeof(int));
  if (config_file){
    fread(buffer, VOLUME, sizeof(int), config_file);
    fclose(config_file);
  } else {
    printf("Could not read configuration\n");
    exit(1);
  }
  for (int t=0; t<NT; t++) for (int x=0; x<NX; x++) 
    field[t][x] = buffer[t*NX+x];
  free(buffer);
}




/* Matrices needed for the fluctuation matrix */
double **BAPe, **PACe, **BGCe, **BAPo, **PACo, **BGCo;
int *newsite, **new_combination;

/* Invert the matrix of propagators between occupied sites
 * Assing the current configuration as the background
 */
double previous_accepted_det = 1;
double previous_det = 1;
double det_save = 1;
int det_sign=1;
int previous_sign=1;
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

 struct timeval start, end;
 gettimeofday(&start,NULL);

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
    printf(" Incorrect determinant, accepted det %g , accepted factor %g \n", det_save, previous_accepted_det);
#endif
    exit(-1);
  }
  LAPACK_dgetrf( &n, &n, Ginv_odd, &n, ipiv_o, &info );
  if( info != 0 ) {
    printf("sgetrf returned error %d (zero determinant has been accepted)! \n", info);
#ifdef DEBUG
    printf(" Incorrect determinant, accepted det %g , accepted factor %g \n", det_save, previous_accepted_det);
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
    //exit(-1);
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

 //printf( "even %4.3g odd %4.3g  n=%d\n", det_e, det_o, n );
#ifdef DEBUG
  double det_diff =  fabs(fabs(det_e*det_o) - det_save);
  printf("EXACT %g %g det %g  accepted %g  diff %g, accepted factor %g \n", det_e, det_o, det_e*det_o, det_save, det_diff, previous_accepted_det);
  if(init == 0 && det_diff/(det_e*det_o)>0.0000001){
    printf(" Incorrect determinant, det %g  accepted %g  diff %g, accepted factor %g \n",fabs(det_e*det_o), det_save, det_diff, previous_accepted_det);
    exit(1);
  }
#endif
  previous_accepted_det = previous_det = 1;
#ifdef DEBUG
  previous_sign = det_sign = (det_e*det_o)<0 ? -1 : 1;
  det_save = fabs(det_e*det_o);
#endif

  gettimeofday(&end,NULL);
  double diff = 1e6*(end.tv_sec-start.tv_sec) + end.tv_usec-start.tv_usec;
  printf("Inverted propagator matrices in %.3g seconds\n", 1e-6*diff);

  /* Mark fluctutation matrix at all sites as not calculated */
  for( int a=0; a<VOLUME; a++ ) newsite[a] = 1;
  for( int a=0; a<VOLUME/2; a++ ) for( int b=0; b<VOLUME/2; b++ ) 
    new_combination[a][b] = 1;

  if(init==1) init = 0;
}


/* Update the links and monomers
 */
int *legalemptylinks;
int *legalmonomers;
int *links;   /* all links are legal */
int n_legalemptylinks=0;
int n_legalmonomers=0;

void update_linklists(){
  static int init=1;
  if(init==1){
    legalemptylinks = malloc(ND*VOLUME*sizeof(int));
    legalmonomers = malloc(ND*VOLUME*sizeof(int));
    links = malloc(ND*VOLUME*sizeof(int));
    init = 0;
  }

  n_links=0; n_legalemptylinks=0; n_legalmonomers=0;
  for (int t=0; t<NT; t++) for (int x=0; x<NX; x++)
    if( field[t][x] >= LINK_TUP && field[t][x] < LINK_TUP+ND  ){    //  Count each link once
      int nu = field[t][x] - LINK_TUP;
      links[n_links] = nu*NX*NT + t*NX + x ;
      n_links ++;
  }

  for (int t=0; t<NT; t++) for (int x=0; x<NX; x++) for (int nu=0; nu<ND; nu++)
    if( (field[t][x] == 0) && (field[tdir(t,nu)][xdir(x,nu)] == 0) ){
      legalemptylinks[n_legalemptylinks] = nu*NX*NT + t*NX + x;
      n_legalemptylinks ++;
  }

  for (int t=0; t<NT; t++) for (int x=0; x<NX; x++) for (int nu=0; nu<ND; nu++)
    if( (field[t][x] == MONOMER) && (field[tdir(t,nu)][xdir(x,nu)] == MONOMER) ){
      legalmonomers[n_legalmonomers] = nu*NX*NT + t*NX + x ;
      n_legalmonomers ++;
  }
#ifdef DEBUG
  printf("update_linklists: n_links = %d, n_legalmonomers = %d, n_legalemptylinks = %d\n",n_links,n_legalmonomers,n_legalemptylinks);
#endif
}

int n_legal_monomers_after_adding(int t1, int x1, int t2, int x2){
  int n_new=1;
  for (int nu=0; nu<NDIRS; nu++){
    if( field[tdir(t1,nu)][xdir(x1,nu)] == MONOMER )
      n_new ++;
    if( field[tdir(t2,nu)][xdir(x2,nu)] == MONOMER )
      n_new ++;
  }
  return n_legalmonomers+n_new;
}

int n_legal_links_after_removing(int t1, int x1, int t2, int x2){
  int n_new=1;
  for (int nu=0; nu<NDIRS; nu++)
  {
    if( field[tdir(t1,nu)][xdir(x1,nu)] == 0 )
      n_new ++;
    if( field[tdir(t2,nu)][xdir(x2,nu)] == 0 )
      n_new ++;
  }
#if NT==2
  if(t2==tdir(t1,0)) n_new += 1;
#endif
#if NX==2
  if(x2==xdir(x1,1)) n_new += 1;
#endif
  //printf("nl %d nn %d\n",n_legalemptylinks,n_new);
  return n_legalemptylinks+n_new;
}



void new_link(int t, int x, int nu){
  link_on(t,x,nu);
#ifdef DEBUG
  det_save = det_save*previous_det/previous_accepted_det;
#endif
  previous_accepted_det = previous_det;
  det_sign = previous_sign;

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
  det_save = det_save*previous_det/previous_accepted_det;
#endif
  previous_accepted_det = previous_det;
  det_sign = previous_sign;

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
  det_save = det_save*previous_det/previous_accepted_det;
#endif
  previous_accepted_det = previous_det;
  det_sign = previous_sign;

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
  printf("Removed link at (%d,%d,%d) %g \n",t,x,nu,det_save);
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
  det_save = det_save*previous_det/previous_accepted_det;
#endif
  previous_accepted_det = previous_det;
  det_sign = previous_sign;

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
  printf("Removed monomer at (%d,%d,%d) %g \n",t1,x1,nu,det_save);
  printf(" indexes %d %d %d %d \n",n_added_even,n_added_odd,n_removed_even,n_removed_odd);
#endif
  update_linklists();
  if(n_removed_even == max_changes || n_removed_odd == max_changes ){
    update_background();
    n_removed_even=n_removed_odd=0;
    n_added_even=n_added_odd=0;
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
    previous_sign = 1;
    previous_det = 1;
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

  previous_sign = det<0 ? -1: 1;
  previous_det = fabs(det);

  return( fabs(det) );
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
  double detratio = det/previous_accepted_det;
#ifdef DEBUG
  printf("Adding at (%d,%d) (%d) and (%d,%d) (%d)\n",t,x,(t*NX+t)/2,t2,x2,(t2*NX+x2)/2);
  printf(" new det %g  %g %g\n",det_save*detratio, det, previous_accepted_det );
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
  double detratio = det/previous_accepted_det;
#ifdef DEBUG
  printf("Removing at (%d,%d) and (%d,%d), values %d and %d\n",t,x,t2,x2,field[t][x],field[t2][x2]);
  printf(" new det %g  %g %g\n",det_save*detratio, det, previous_accepted_det );
#endif

  return( detratio );
}




double det_moved_monomer(int t, int x, int t2, int x2){
  int n = n_bc_monomers/2 + n_bc_links;
  double detratio;

  /* Make sure that (t,x) is even and (t2,x2) is odd */
  if( (t+x)%2 != (t2+x2)%2) { 
    //printf("Attempting to move to opposite parity\n");
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
    detratio = det/previous_accepted_det;
#ifdef DEBUG
    printf("Moving (%d,%d) (%d) to (%d,%d) (%d)\n",t,x,(NX*t+x)/2,t2,x2,(NX*t2+x2)/2);
    printf(" new det %g  %g %g\n",det_save*detratio, det, previous_accepted_det );
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
    detratio = det/previous_accepted_det;
#ifdef DEBUG
    printf("Moving (%d,%d) (%d) to (%d,%d) (%d)\n",t,x,(NX*t+x)/2,t2,x2,(NX*t2+x2)/2);
    printf(" new det %g  %g %g\n",det_save*detratio, det, previous_accepted_det );
#endif
  }
  return( detratio );
}








/* Move monomer along a set of links  */
int move_monomer(){
  int success = 0;
  if( n_monomers > 0 ){

    /* Pick monomer  */
    int i=0,t,x,m = mersenne()*n_monomers + 1 ;
    for (t=0; i<m; t++) for (x=0; x<NX && i<m ; x++) {
      if( field[t][x] == MONOMER ) i++;
    }
    t--;x--;

    for( int iter = 0; mersenne()>flip_exit_propability; iter++ ){
      /* pick direction, check for links and flip */
      int nu = mersenne()*NDIRS;
      int t2, x2, linkdir;
      t2 = tdir(t,nu); x2=xdir(x,nu);

      if( field[t2][x2] >= LINK_TUP ) {
        linkdir = field[t2][x2] - LINK_TUP;
        int t3 = tdir(t2,linkdir), x3 = xdir(x2,linkdir);

        field[t][x] = 0;
        link_off(t2,x2,linkdir); link_on(t,x,nu);
        field[t3][x3] = MONOMER;
      
        t = t3; x = x3;
        success = 1;
      }
    }
    
  }
  update_linklists();
  return success;
}


/* Suggest adding a link
 */
int add_link()
{
  int success = 0;
  if(n_legalemptylinks > 0){

    /* Draw a legal site */
    int s = legalemptylinks[ (int) (mersenne()*n_legalemptylinks) ];
    int t, x, nu;
    x = s % NX; t = (s/NX)%NT; nu = s/(NX*NT);
    int t2,x2;
    t2 = tdir(t,nu); x2 = xdir(x,nu);

    /* Calculate the determinant fraction */
    double d = det_added_link(t,x,t2,x2);

    if( mersenne() < 0.5 ) {
      double p=( n_legalemptylinks/(double)n_legal_monomers_after_adding(t,x,t2,x2) );
      p*= m*m*d;
      if( mersenne() < p ) {
        new_monomer(t,x,t2,x2);
        success = 1;
      }
    } else {
      double p=( n_legalemptylinks/((double)n_links+1) ) * U*d;
      //printf(" Adding, p %g le %d ll %d\n", p, n_legalemptylinks, n_links+1);
      if( mersenne() < p ) {
        new_link(t,x,nu);
        success = 1;
      }
    }
  }
  
  return success;
}

int remove_link()
{
  int success = 0;
  /* Try removing either a link or a monomer */
  if( mersenne() < 0.5 ){
    if(n_legalmonomers > 0){
      int s = legalmonomers[ (int)(mersenne()*n_legalmonomers) ];
      int t, x, nu;
      x = s % NX; t = (s/NX)%NT; nu = s/(NX*NT);
      int t2,x2;
      t2 = tdir(t,nu); x2 = xdir(x,nu);

      double M = det_removed_link(t,x,t2,x2);
      double p = ( n_legalmonomers/((double)n_legal_links_after_removing(t,x,t2,x2)) );
      p *= M/(m*m);
      if( mersenne() < p ){
        removed_monomer(t,x,nu);
        success = 1;
      }
    }
  } else {
    if( n_links > 0){
      int s = links[ (int)(mersenne()*n_links) ];
      int t, x, nu;
      x = s % NX; t = (s/NX)%NT; nu = s/(NX*NT);
      int t2,x2;
      t2 = tdir(t,nu); x2 = xdir(x,nu);

      double M = det_removed_link(t,x,t2,x2);
      double p = ( n_links/((double) n_legal_links_after_removing(t,x,t2,x2) ) ) * M/U;
      //printf(" Removing, p %g det %g ll %d le %d\n",p,M*M,n_links,n_legal_links_after_removing(t,x,t2,x2));
      if( mersenne() < p ){
        removed_link(t,x,nu);
        success = 1;
      }
    }
  }
  return success;
}



int update()
{
  int changes=0;
  if(mersenne()>0.5){
    changes+=remove_link();
  } else {
    changes+=add_link();
  }
  #ifdef DEBUG
  //check_det( );
  //print_config();
  #endif
  move_monomer();
  return changes;
}




/* Functions for measurements
 */

/* Propagator
 */
void measure_propagator(){ 
  double prop[NT];
  double j[NT];
  for( int t1=0; t1<NT; t1++){
    prop[t1]=0; j[t1] = 0;
  }

  double **source, **propagator;
  source = alloc_vector();
  propagator = alloc_vector();
  //vec_zero( source );

  //int t1=1;
  for( int t1=0;t1<NT;t1++)  for( int x1=0;x1<NX;x1++) if( field[t1][x1] == 0 )  {
      vec_zero( source );
      //for( int x1=0;x1<NX;x1++) if( field[t1][x1] == 0 ) source[t1][x1] = (mersenne()>0.5) ? -1:1 ;
      source[t1][x1] = 1;

      vec_zero( propagator );
      cg_propagator(propagator,source);

      /*if( t1==0 && x1==0 ){
        printf("Dinv %g\n",propagator[0][1]);
        printf("Dinv %g\n",propagator[1][0]);
        printf("Dinv %g\n",propagator[3][0]);
        //printf("Dinv %g\n",propagator[5][0]);
      }*/

      //for( int x1=0;x1<NX;x1++)
      for( int t2=0; t2<NT; t2++)
        prop[(t2-t1+NT)%NT] += propagator[t2][x1];

      //for( int x1=0;x1<NX;x1++)
      if( field[t1][x1] == 0 )
        j[t1] += exp(mu)*eta[t1][x1][0]*propagator[tup[t1]][x1];
      //for( int x1=0;x1<NX;x1++)
      if( field[t1][x1] == 0 ) 
        j[tdn[t1]] += exp(-mu)* eta[tdn[t1]][x1][0]*propagator[tdn[t1]][x1];
      //source[t1][x1] = 0;
  }
  j[NT-1] = -j[NT-1]; //Fix boundary condition
  
  for( int t2=0; t2<NT; t2++) printf("Propagator %d %g\n", t2, det_sign*prop[t2]/(VOLUME) );
  for( int t1=0;t1<NT;t1++)
    printf("Charge %d %g\n", t1, det_sign*j[t1]/2 );

  free_vector(source);
  free_vector(propagator);
}



/* Measure the susceptibility using a worm algorithm. Introduce 2 source monomers
 * with the weigth J (=U*NDIRS/V). Allow one to move around (produce configurations)
 * until it contacts with the the other one, remove with the appropriate weight.
 * The number of intermediate configurations counts the Z_J/Z_0, which is the 
 * susceptibility , Z_J/Z_0 = U*NDIRS/V * dZ_J/dJ |_J=0.
 */
void measure_susceptibility(){
 int n = n_bc_monomers/2 + n_bc_links;
 int steps = 0;
 int n_attempts=100;

 /* Do multiple attemps, these are cheap and the result is usually 0 */
 update_linklists();
 if(n_links > 0) for( int attempt=0; attempt<n_attempts; attempt++ ){
   /* Pick a site with a link */
   int s = links[ (int)(mersenne()*n_links) ];
   int t1, x1, nu, t2,x2;
   x1 = s % NX; t1 = (s/NX)%NT; nu = s/(NX*NT);
   if(mersenne()>0.5) {  //pick one of the two sites
     t1 = tdir(t1,nu); x1 = xdir(x1,nu);
   }

   /* Chose a link, switch it to a pair of source monomers */
   int dir = field[t1][x1] - LINK_TUP;
   link_off(t1,x1,dir);
    
   t2 = tdir(t1,dir), x2 = xdir(x1,dir);
   field[t1][x1] = SOURCE_MONOMER ; field[t2][x2] = SOURCE_MONOMER ;
   
   steps+=det_sign;

   #ifdef DEBUG
   print_config();
   check_det();
   #endif
   
   for(;;){
     /* Now we are at (t2,x2), and the link is off. Try to move. */
     int dir = NDIRS*mersenne();
     int t3 = tdir(t2,dir), x3 = xdir(x2,dir);

     if( mersenne() > 0.5 ){
      /* Try to exchange with an occupied neighbor */
      if( t3==t1 && x3==x1 ) {
        /* Back at the original site, turn into a link */
        field[t1][x1]=0; field[t2][x2]=0;
        link_on(t2,x2,dir);
        break;

      } else if( field[t3][x3] >= LINK_TUP ) {
        /* found a link, flip it */
        int linkdir = field[t3][x3] - LINK_TUP;
        int t4 = tdir(t3,linkdir), x4 = xdir(x3,linkdir);

        field[t2][x2] = 0;
        link_off(t3,x3,linkdir); link_on(t2,x2,dir);
        field[t4][x4] = SOURCE_MONOMER;

        t2 = t4; x2 = x4;

      } else if( field[t3][x3] == MONOMER ) {
        /* Monomer found, change places */
        field[t2][x2] = MONOMER; field[t3][x3] = SOURCE_MONOMER;

        t2 = t3; x2 = x3;
      }
      #ifdef DEBUG
      print_config();
      #endif
     } else {
      /*  Try to hop over */
      int dir2; 
      do dir2 = NDIRS*mersenne();
      while( dir2==opp_dir(dir));
      int t4 = tdir(t3,dir2), x4 = xdir(x3,dir2);

      if( field[t4][x4] == 0 ) {
        /* Check for the chosen sites in added and removed lists */
        int s = (t2*NX+x2)/2;
        int old_site=1;
        if( (t2+x2)%2==0 ) for( int k=0; k<n_added_even; k++ ) if( added_evensites[k] == s )
           old_site = 0;
        if( (t2+x2)%2==1 ) for( int k=0; k<n_added_odd; k++ ) if( added_oddsites[k] == s )
           old_site = 0;

        s = (t4*NX+x4)/2;
        int new_site=1;
        if( (t2+x2)%2==0 ) for( int k=0; k<n_removed_even; k++ ) if( evenlist[removed_evenlist[k]] == s )
           new_site = 0;
        if( (t2+x2)%2==1 ) for( int k=0; k<n_removed_odd; k++ ) if( oddlist[removed_oddlist[k]] == s )
           new_site = 0;

        /* Get the difference in the determinant */
        double det;
        det = det_moved_monomer( t2, x2, t4, x4 );

        if( mersenne() < det ){
          /* Accepted */
          int s = (t2*NT+x2)/2;
          if( old_site ) {
            if( (t2+x2)%2==0 ) { /* The lists are already ordered properly, so just adjust lengths */
              n_added_even++;
              n_removed_even++;
            } else {
              n_removed_odd++;
              n_added_odd++;
            }
          }
          if( !new_site ) {   // Found in the list of removed sites
            if( (t2+x2)%2==0 ) {
              n_added_even--;
              n_removed_even--;
            } else {
              n_removed_odd--;
              n_added_odd--;
            }
          }

          field[t2][x2] = 0; field[t4][x4] = SOURCE_MONOMER;
          t2 = t4; x2 = x4;
#ifdef DEBUG
          det_save = det_save*previous_det/previous_accepted_det;
#endif
          previous_accepted_det = previous_det;
          det_sign = previous_sign;

          #ifdef DEBUG
          print_config();
          check_det();
          #endif
          //if(n_added_even == max_changes || n_added_odd == max_changes ||
          //   n_removed_even == max_changes || n_removed_odd == max_changes ){
          if(n_added_even == 2 || n_added_odd == 2 ||
             n_removed_even == 2 || n_removed_odd == 2 ){
            update_background();
            n_removed_even=n_removed_odd=0;
            n_added_even=n_added_odd=0;
          }
        } else {
          /* Rejected */
          if(!old_site){
            if( (t2+x2)%2==0 ) {
              added_evensites[n_added_even-1] = (t2*NX+x2)/2;
            } else {
              added_oddsites[n_added_odd-1] = (t2*NX+x2)/2;
            }
          }
          if( !new_site ) {
            if( (t2+x2)%2==0 ) {
              for( int k=0; k<n; k++ ) if( evenlist[k] == s )
                removed_evenlist[n_removed_even-1] = k;
            } else {
              for( int k=0; k<n; k++ ) if( oddlist[k] == s )
                removed_oddlist[n_removed_odd-1] = k;
            }
          }
        }
      }
    } //neighbouring site (t3,x3)

    steps+=det_sign;
   } //steps
   update_linklists();
 } //attempts


 printf("Susceptibility %g \n",(double)steps*n_links/(U*2*NDIRS*n_attempts*VOLUME));
  
}




void print_config()
{
  for (int t=0; t<NT; t++) {
    for (int x=0; x<NX; x++){
      int empty = 1;
      if(field[t][x]==MONOMER)  { empty = 0; printf(" o "); }
      if(field[t][x]==LINK_TUP) { empty = 0; printf(" | "); }
      if(field[t][x]==LINK_XUP) { empty = 0; printf(" --"); }
      if(field[t][x]==LINK_TDN) { empty = 0; printf(" | "); }
      if(field[t][x]==LINK_XDN) { empty = 0; printf("-- "); }
      if(field[t][x]==SOURCE_MONOMER) { empty = 0; printf(" x "); } //A source monomer
      if(empty==1) { printf(" . "); }
    }
    printf(" \n");
  }
  printf(" \n");
  //usleep(100000);
}

static int measurement = 0;
void measure()
{ 
#ifdef DEBUG
  //print_config();
#endif

  measure_propagator();
  measure_susceptibility();

  measurement++;
}







/* Main function
 */
int main(int argc, char* argv[])
{
  #ifdef DEBUG
  feenableexcept(FE_INVALID | FE_OVERFLOW);
  #endif 

  int i,n_loops,n_measure,n_thermalize;
  long seed;

  /* Read in the input */
  printf(" Number of updates : ");
  scanf("%d",&n_loops);

  printf(" Updates / measurement : ");
  scanf("%d",&n_measure);
  
  printf(" Thermalization updates : ");
  scanf("%d",&n_thermalize);
  printf("%d\n",n_thermalize);

  printf(" m : ");
  scanf("%lf",&m);

  printf(" U : ");
  scanf("%lf",&U);

  printf(" mu : ");
  scanf("%lf",&mu);

  printf(" Random number : ");
  scanf("%ld",&seed);
  seed_mersenne( seed );

  printf(" maximum size of fluctuation matrix : ");
  scanf("%d",&max_changes);

  char start_config[100];
  printf(" Start configuration : ");
  scanf("%s",&start_config);


  /* "Warm up" the rng generator */
  for (i=0; i<543210; i++) mersenne();

  printf(" \n++++++++++++++++++++++++++++++++++++++++++\n");
  printf(" 4D free fermion, ( %d , %d ) lattice\n", NT, NX );
  printf(" %d updates per measurements\n", n_measure );
  printf(" m %f \n", m);
  printf(" U %f \n", U);
  printf(" mu %f \n", mu);
  printf(" Size of fluctuation matrix %d\n", max_changes );
  printf(" Random seed %ld\n", seed );
  printf(" Start configuration %s\n", start_config );

  field = malloc( NT*sizeof(int *) );
  eta = malloc( NT*sizeof(int *) );
  for (int t=0; t<NT; t++){
    field[t] = malloc( (NX+1)*sizeof(int) );
    eta[t] = malloc( (NX+1)*sizeof(int *) );
    for (int x=0; x<NX+1; x++) {
     eta[t][x] = malloc( 2*sizeof(int) );
    }
  }
  xup = malloc( (NX+1)*sizeof(int) );
  xdn = malloc( (NX+1)*sizeof(int) );
  tup = malloc( NT*sizeof(int) );
  tdn = malloc( NT*sizeof(int) );

  if(strcmp(start_config,"cold")==0){
    printf(" Starting from a cold configuration\n" );
    for (int t=0; t<NT; t++) for (int x=0; x<NX; x++) {
      field[t][x] = 0;
    }
  } else {
    printf(" Reading configuration file\n" );
    read_config(start_config);
  }
#ifdef OPENX
  for (int t=0; t<NT; t++) {
    field[t][NX] = EMPTY; //Site doesn't exist, no links or monomers, but not free either
  }
#endif


  /* allocate propagator and lists */
  Dinv = malloc( VOLUME*VOLUME*sizeof(double)/2 );
  Ginv = malloc( VOLUME*VOLUME*sizeof(double)/2 );
  evenlist = malloc( VOLUME*sizeof(int)/2 );
  oddlist  = malloc( VOLUME*sizeof(int)/2 );
  unoccupied_evenlist = malloc( VOLUME*sizeof(int)/2 );
  unoccupied_oddlist  = malloc( VOLUME*sizeof(int)/2 );

  added_evensites = malloc(VOLUME/2*sizeof(int));
  added_oddsites = malloc(VOLUME/2*sizeof(int));
  removed_evenlist = malloc(VOLUME/2*sizeof(int));
  removed_oddlist = malloc(VOLUME/2*sizeof(int));
  
  if (NULL == Dinv) {
    fprintf(stderr, "failed to allocate Dinv\n");
    return(-1);
  }
  if (NULL == Ginv) {
    fprintf(stderr, "failed to allocate Ginv\n");
    return(-1);
  }
  if (NULL == evenlist) {
    fprintf(stderr, "failed to allocate evenlist\n");
    return(-1);
  }
  if (NULL == oddlist) {
    fprintf(stderr, "failed to allocate oddlist\n");
    return(-1);
  }

  /* fill up the index array */
  for (i=0; i<NT; i++) {
    tup[i] = (i+1) % NT;
    tdn[i] = (i-1+NT) % NT;
  }
  for (i=0; i<NX+1; i++) {
    xup[i] = (i+1) % NX;
    xdn[i] = (i-1+NX) % NX;
  }
#ifdef OPENX  //Never match boundaries to actual sites
  xdn[0] = NX;
  xup[NX-1] = NX;
  xup[NX] = NX;
#endif

  /* fill the staggered eta array */
  for (int t=0; t<NT; t++) for (int x=0; x<NX; x++) {
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

  
  /* calculate propagators */
  calc_Dinv( );
  int changes = 0;
  struct timeval start, end;
  gettimeofday(&start,NULL);

  /* Background and fluctuation matrix */
  update_linklists();
  update_background( );

/*
  double detratio= det_added_link( 0, 0, 0, 1);

  printf(" %g ",detratio);

  new_link(0 , 0, 1);

  detratio= det_added_link( 1, 2, 1, 3);

  printf(" %g ",detratio);
*/

  /*print_config();
  int t2 = 63, x2 = 18;
  for(int mu=0; mu<4; mu++) for(int nu=0; nu<=mu; nu++){
    printf("  mu %d nu %d\n",mu,nu);
    int t3,x3,t4,x4;
    t3 = tdir(t2,mu); x3 = xdir(x2,mu);
    t4 = tdir(t3,nu); x4 = xdir(x3,nu);
    det_moved_monomer( t2, x2, t4, x4 );
  }
  exit(0);  
*/
  /* and the update/measure loop */
  for (i=1; i<n_loops+1; i++) {

    /* Update */
    changes += update();

    if((i%n_measure)==0){

      /* Update the background and count links and monomers */
      update_background();
      n_added_even = n_added_odd = 0;
      n_removed_even = n_removed_odd = 0;
      update_linklists();
      
      /* Time and report */
      gettimeofday(&end,NULL);
      int diff = 1e6*(end.tv_sec-start.tv_sec) + end.tv_usec-start.tv_usec;
      printf("\n%d, %d updates in %.3g seconds, %d successfull changes, %g changes/update\n", i, n_measure, 1e-6*diff,changes,changes/((double) i));
      
      print_config();

      /* Statistics */
      printf("MONOMERS %d \n", det_sign*n_monomers);
      printf("LINKS %d \n", det_sign*n_links);
      #ifdef DEBUG
      printf("Determinant %g \n", det_sign*det_save);
      #endif DEBUG
      printf("Sign %d \n", det_sign);

      /* Write configuration */
      //write_config();

      gettimeofday(&start,NULL);

      if(i>n_thermalize) {
        /* Do measurements */
        measure();

        /* Time measurements */
        gettimeofday(&end,NULL);
        diff = 1e6*(end.tv_sec-start.tv_sec) + end.tv_usec-start.tv_usec;
        printf("Measurement %d done in %.3g seconds\n", i/n_measure, 1e-6*diff);
        gettimeofday(&start,NULL);
      }
    }
  }

  printf(" ** simulation done\n");

  return(1);
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
    printf(" Incorrect determinant, accepted det %g , accepted factor %g \n", det_save, previous_accepted_det); 
#endif
    exit(-1);
  }
  LAPACK_dgetrf( &n, &n, check_Ginv_odd, &n, ipiv_o, &info );
  if( info != 0 ) {
    printf("sgetrf returned error %d (zero determinant has been accepted)! \n", info);
#ifdef DEBUG
    printf(" Incorrect determinant, accepted det %g , accepted factor %g \n", det_save, previous_accepted_det); 
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

 printf( "even %4.3g odd %4.3g  n=%d\n", det_e, det_o, n );
 double diff =  fabs(fabs(det_e*det_o) - det_save);
#ifdef DEBUG
 printf("CHECK det %g  accepted %g  diff %g, accepted factor %g \n", det_e*det_o, det_save, diff, previous_accepted_det);
 if(diff/fabs(det_e*det_o)>0.00001){
    printf(" Incorrect determinant, accepted det %g , accepted factor %g \n", det_save, previous_accepted_det); 
   exit(1);
 }
#endif

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














