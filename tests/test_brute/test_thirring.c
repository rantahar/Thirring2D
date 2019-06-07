/****************************************************************
 * A brute force integration of the partition function for the
 * 2D Thirring model. Usefull for testing the correctness
 * of a simulation algorithm.
 * 
 * Go through all possible configurations and calculate the
 * determinant and measurables for each one.
 *
 ****************************************************************/



#include "test_thirring.h"

/* storage */
int    ***eta;   //Staggered eta matrix
double *Dinv;             //Storage for even to odd propagator
int    *evenlist, *oddlist;  //Lists of occupied sites
int    *unoccupied_evenlist, *unoccupied_oddlist;  //Lists of occupied sites
double m;
double U;
double mu;


/* Monomers and links
 * field stores both, 0 for empty, 1 for monomer and 2+dir for links
 */
int n_monomers=0;
int n_links=0;
int **field;
int n_bc_monomers;  //Numbers of links and monomers in the background
int n_bc_links;

/* Inverse of G for the background config
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

/* Turn two monomer on around a link */
static inline void monomer_on(int t, int x, int t2, int x2){
  if ( field[t][x] == 0 && field[t2][x2] == 0 ){
    field[t][x] = MONOMER;
    field[t2][x2] = MONOMER;
#ifdef DEBUG
    printf("Turned on monomers (%d,%d) and (%d,%d)\n",t,x,t2,x2);
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


/* Check if we're allowed to add a link or monomer */
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
  
  free(ipiv);
}


/* Invert the matrix of propagators between occupied sites
 * Assign the current configuration as the background
 */
double previous_accepted_det = 1;
double previous_det = 1;
double det_save = 1;
void update_background( )
{

 static int init = 1;

 /* The current configuration becomes the background */
 n_bc_monomers = n_monomers;
 n_bc_links = n_links;

 double det_e=1, det_o=1;
 det_save = 1; 
 int n = n_monomers/2 + n_links;
 double *Ginv_odd = Ginv+n*n;
 double *Dinv_odd = Dinv+VOLUME*VOLUME/4;

 struct timeval start, end;
 gettimeofday(&start,NULL);

 /*If no occupied sites, determinant of rank 0 is 1 */
 if( n_monomers + n_links > 0 ){
  int *ipiv,*ipiv_o;
  int info;
  
  ipiv = malloc( n*sizeof(int) ); ipiv_o = malloc( n*sizeof(int) );
  
  /* Find occupied and unoccupied sites, construct lists */
  int i=0,j=0;
  for (int t=0; t<NT; t++) for (int x=0; x<NX; x++){
    if(field[t][x]>0) {
      if((t+x)%2==0) {
        evenlist[i] = (NX*t + x)/2;
        i++;
      } else {
        oddlist[j] = (NX*t + x)/2;
        j++;
      }
    }
  }

  if( n_monomers%2 == 0 && i==j ){
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
    LAPACK_dgetrf( &n, &n, Ginv_odd, &n, ipiv_o, &info );

    /* determinant from LU */
    for(int i=0; i<n; i++) {
      if(ipiv[i]==i+1) { det_e *= Ginv[i*n+i];}
      else { det_e *= -Ginv[i*n+i]; }
      if(ipiv_o[i]==i+1) { det_o *= -Ginv_odd[i*n+i];}
      else { det_o *= Ginv_odd[i*n+i]; }
    }
  
    det_save = det_e*det_o;
  } else {
    det_save = 0;
  }
 }

 gettimeofday(&end,NULL);
 double diff = 1e6*(end.tv_sec-start.tv_sec) + end.tv_usec-start.tv_usec;
 //printf("Inverted propagator matrices in %.3g seconds, det=%g\n", 1e-6*diff,det_save);
 //if(det_save==0) sleep(1);
 if(init==1) init = 0;
}


/* Update the links and monomers
 */
int legalemptylinks[ND*VOLUME];
int legalmonomers[ND*VOLUME];
int links[ND*VOLUME];   /* all links are legal */
int n_legalemptylinks=0;
int n_legalmonomers=0;

void update_linklists(){
  n_links=0; n_monomers=0;
  for (int t=0; t<NT; t++) for (int x=0; x<NX; x++)
    if( field[t][x] >= LINK_TUP && field[t][x] < LINK_TUP+ND  ){    //  Count each link once
      int nu = field[t][x] - LINK_TUP;
      n_links ++;
  }

  for (int t=0; t<NT; t++) for (int x=0; x<NX; x++)
    if(field[t][x] == MONOMER) {
      n_monomers ++;
  }
#ifdef DEBUG
  printf("update_linklists: n_links = %d, n_legalmonomers = %d, n_legalemptylinks = %d\n",n_links,n_legalmonomers,n_legalemptylinks);
#endif
}








/* Print a configuration */
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

/* Measure the fermion and chiral currents (j and cj) and the chiral charge squared (cs)
 */
void measure_charge( double * j, double * cj, double * cs){ 

  double **source, **propagator;
  source = alloc_vector();
  propagator = alloc_vector();

  double cj1=0;

  for( int x1=0;x1<NX;x1++)  {
    int t1=0;
    if(field[t1][x1]==0){
      vec_zero( source );
      source[t1][x1] = 1;
      vec_zero( propagator );
      cg_propagator(propagator,source);

      *j += exp(mu)*eta[t1][x1][0]*propagator[tup[t1]][x1];
      cj1 += ((x1+t1)%2 ==0 ? 1:-1 )* exp(mu)*eta[t1][x1][0]*propagator[tup[t1]][x1];
    }
    if( field[t1][x1]== LINK_TUP ) {
      double cl = ((t1+x1)%2==0 ? 1:-1 )*2;
      cj1 += 2*cl;
    }
   
    t1=1;
    if(field[t1][x1]==0){
      vec_zero( source );
      source[t1][x1] = 1;
      vec_zero( propagator );
      cg_propagator(propagator,source);

      *j += exp(-mu)* eta[tdn[t1]][x1][0]*propagator[tdn[t1]][x1];
      cj1 -= ((x1+tdn[t1])%2 ==0 ? 1:-1 )* exp(-mu)* eta[tdn[t1]][x1][0]*propagator[tdn[t1]][x1];
    }
  }

  *j/=2;
  *cj = cj1/2;

  *cs=cj1*cj1/4;

  free_vector(source);
  free_vector(propagator);
}


static double ZL[VOLUME/2+1][VOLUME+1];
static double ZLabs[VOLUME/2+1][VOLUME+1];
static double j[VOLUME/2+1][VOLUME+1];
static double cj[VOLUME/2+1][VOLUME+1];
static double cs[VOLUME/2+1][VOLUME+1];


/* Go through all possible configurations */
void update( int i0 )
{
  
  //printf(" Links %d, Monomers %d, Determinant %g \n", n_links, n_monomers, det_save*det_save );
  ZL[n_links][n_monomers] += det_save;
  ZLabs[n_links][n_monomers] += fabs(det_save);

  double jm=0, cjm=0, csm=0;
  measure_charge(&jm,&cjm,&csm);
  j[n_links][n_monomers] += det_save*jm;
  cj[n_links][n_monomers] += det_save*cjm;
  cs[n_links][n_monomers] += det_save*csm;

  /* Starting from top left, set each site to each possible state */
  for( int i=i0; i< VOLUME ; i++ ) {
    int t = (i%VOLUME)/NT, x = i%NT;
    if(field[t][x]==0){
      for(int nu=0; nu<2; nu++) {
        int t2 = tdir(t,nu), x2 = xdir(x,nu);
        if(field[t2][x2]==0) {
          link_on(t,x,nu);
          update_linklists();
          update_background();

          /* Recursively update the rest of the lattice */
          update(i);
          link_off(t,x,nu);
        }
      }
    
      field[t][x] = 1;
      update_linklists();
      update_background();

      update(i);
      field[t][x]=0;
    }
  }
}






/* Main function
 */
int main(int argc, char* argv[])
{
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

  printf(" \n++++++++++++++++++++++++++++++++++++++++++\n");
  printf(" 4D free fermion, ( %d , %d ) lattice\n", NT, NX );
  printf(" %d updates per measurements\n", n_measure );
  printf(" m %f \n", m);
  printf(" U %f \n", U);
  printf(" mu %f \n", mu);
  printf(" Random seed %ld\n", seed );


  /* allocate propagator and lists */
  Dinv = malloc( VOLUME*VOLUME*sizeof(double)/2 );
  Ginv = malloc( VOLUME*VOLUME*sizeof(double)/2 );
  evenlist = malloc( VOLUME*sizeof(int)/2 );
  oddlist  = malloc( VOLUME*sizeof(int)/2 );
  unoccupied_evenlist = malloc( VOLUME*sizeof(int)/2 );
  unoccupied_oddlist  = malloc( VOLUME*sizeof(int)/2 );
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

  /* fill up the index array */
  for (i=0; i<NT; i++) {
    tup[i] = (i+1) % NT;
    tdn[i] = (i-1+NT) % NT;
  }

  for (i=0; i<NX; i++) {
    xup[i] = (i+1) % NX;
    xdn[i] = (i-1+NX) % NX;
  }
#ifdef OPENX  //Never match boundaries to actual sites
  xdn[0] = NX;
  xup[NX-1] = NX;
  xup[NX] = NX;
#endif


  /* Fill the staggered eta array */
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

  
  /* Fill monomers and links */
  for (int t=0; t<NT; t++) for (int x=0; x<NX; x++) {
    field[t][x] = 0;
  }
#ifdef OPENX
  for (int t=0; t<NT; t++) {
    field[t][NX] = EMPTY; //Site doesn't exist, no links or monomers, but not free either
  }
#endif

  update_linklists();

  /* Calculate propagators */
  calc_Dinv( );
  det_save = 1;
  
  /* Set the start configuration as the background configuration */
  update_background( );

  /* Initialize measurements */
  for(int nl=0; nl<VOLUME/2+1; nl++ ) for(int nm=0; nm<VOLUME+1; nm++ ){
    ZL[nl][nm] = 0;
    ZLabs[nl][nm] = 0;
    j[nl][nm] = 0;
    cj[nl][nm] = 0;
    cs[nl][nm] = 0;
  }

  /* And finally run the integral (starting from site 0) */
  update(0);
  
  for(int nl=0; nl<VOLUME/2+1; nl++ ) for(int nm=0; nm<VOLUME+1; nm++ )
    printf(" Links %d Monomers %d det %g abs(det) %g det*charge %g det*chiralcharge %g det*chiralsusc %g\n",nl,nm,ZL[nl][nm],ZLabs[nl][nm],j[nl][nm],cj[nl][nm],cs[nl][nm]);

  printf(" ** simulation done\n");

  return(1);
}









