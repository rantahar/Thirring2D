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
int    eta[NT][NX][ND];   //Staggered eta matrix
int    *evenlist, *oddlist;  //Lists of occupied sites
double m;
double U;
double mu;


/* Monomers and links
 * field stores both, 0 for empty, 1 for monomer and 2+dir for links
 */
int n_monomers=0;
int n_links=0;
int field[NT][NX];
int n_bc_monomers;  //Numbers of links and monomers in the background
int n_bc_links;

/* Pseudofermion fields */
double *psi,*chi;

/* Neighbour index arrays, to be filled at the beginning
 */
int tup[NT],xup[NX],tdn[NT],xdn[NX];

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

    double S1,S2, edS=0;
    vec_gaussian(psi);
    fM_occupied( chi, psi ); //Now chi is from correct distribution
    S1 = action(psi);

    monomer_on(t,x,t2,x2); //make suggested change
    if( cg_MdM_occupied( psi, chi ) ==0 ){
      S2 = action(psi);
    } else {
      S2 = 1e200; //Unable to invert, zero mode
    }
    edS = exp(S1-S2);
    //printf(" action %g %g %g\n",S1,S2,edS);
    link_off(t,x,nu); //change back
    
    /* Calculate the determinant fraction */
    double d = edS;
    if( mersenne() < 0.5 ) {
      double p=( n_legalemptylinks/(double)n_legal_monomers_after_adding(t,x,t2,x2) );
      p*= m*m*d;
      if( mersenne() < p ) {
        n_monomers += 2;
        monomer_on(t,x,t2,x2);
        update_linklists();
        success = 1;
      }
    } else {
      double p=( n_legalemptylinks/((double)n_links+1) ) * U*d;
      //printf(" Adding, p %g le %d ll %d\n", p, n_legalemptylinks, n_links+1);
      if( mersenne() < p ) {
        link_on(t,x,nu);
        update_linklists();
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

      double S1,S2,edS=0;
      vec_gaussian(psi);
      fM_occupied( chi, psi );//Now chi is from correct distribution
      S1 = action(psi);

      int f1=field[t][x], f2=field[t2][x2];
      field[t][x] = 0; field[t2][x2] = 0; //make suggested change
      if( cg_MdM_occupied( psi, chi ) ==0 ){
        S2 = action(psi);
      } else {
        S2 = 1e200; //Unable to invert, zero mode
      }
      edS = exp(S1-S2);
      field[t][x] = f1; field[t2][x2] = f2; //back

      double p = ( n_legalmonomers/((double)n_legal_links_after_removing(t,x,t2,x2)) );
      p *= edS/(m*m);
      if( mersenne() < p ){
        n_monomers -= 2;
        link_off(t,x,nu);
        update_linklists();
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

      double S1,S2,edS=0;
      vec_gaussian(psi);
      fM_occupied( chi, psi );//Now chi is from correct distribution
      S1 = action(psi);

      int f1=field[t][x], f2=field[t2][x2];
      field[t][x] = 0; field[t2][x2] = 0; //make suggested change
      if( cg_MdM_occupied( psi, chi ) ==0 ){
        S2 = action(psi);
      } else {
        S2 = 1e200; //Unable to invert, zero mode
      }
      edS = exp(S1-S2);
      field[t][x] = f1; field[t2][x2] = f2; //back

      double p = ( n_links/((double) n_legal_links_after_removing(t,x,t2,x2) ) ) * edS/U;
      //printf(" Removing, p %g det %g ll %d le %d\n",p,M*M,n_links,n_legal_links_after_removing(t,x,t2,x2));
      if( mersenne() < p ){
        n_links -= 1;
        link_off(t,x,nu);
        update_linklists();
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
  check_det( );
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

  double source[NX][NT], propagator[NT][NX];
  vec_zero( source );

  for( int t1=0;t1<NT;t1++) for( int x1=0;x1<NX;x1++) if( field[t1][x1] == 0 )
  {
    source[t1][x1] = 1;
    cg_propagator(propagator,source);
    
    for( int t2=0; t2<NT; t2++)
      prop[(t2-t1+NT)%NT] += propagator[t2][x1];

    j[t1] += propagator[tup[t1]][x1];
    j[tdn[t1]] += propagator[tdn[t1]][x1];

    source[t1][x1] = 0;
  }
  
  for( int t2=0; t2<NT; t2++) printf("Propagator %d %g\n", t2, prop[t2]/(VOLUME) );
  for( int t1=0;t1<NT;t1++)  printf("Current %d %g\n", t1, j[t1]/(2*NX) );
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
 for( int attempt=0; attempt<n_attempts; attempt++ ){
  /* Pick a site  */
  int s1 = mersenne()*VOLUME ;
  int t1 = s1/NX, x1=s1%NX, t2,x2;

  //printf("measure_susceptibility: Initial site %d, (%d,%d), field %d\n",s1,t1,x1,field[t1][x1]);
  if(field[t1][x1] >= LINK_TUP ){
   /* Chose a link, switch it to a pair of source monomers */
   int dir = field[t1][x1] - LINK_TUP;
   link_off(t1,x1,dir); 
   
   t2 = tdir(t1,dir), x2 = xdir(x1,dir);
   field[t1][x1] = SOURCE_MONOMER ; field[t2][x2] = SOURCE_MONOMER ;
   
   steps++;
   
   for(;; steps++ ){
     //print_config();
     //printf("measure_susceptibility: At site (%d,%d), field %d\n",t2,x2,field[t2][x2]);
     /* Now we are at (t2,x2), and the link is off. Try to move. */
     int dir = NDIRS*mersenne();
     int t3 = tdir(t2,dir), x3 = xdir(x2,dir);

     //printf("measure_susceptibility: Trying site (%d,%d), field %d\n",t3,x3,field[t3][x3]);

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
     } else {
      /*  Try to hop over */
      int dir2 = NDIRS*mersenne();
      int t4 = tdir(t3,dir2), x4 = xdir(x3,dir2);

      //printf("measure_susceptibility: Hop over, new site (%d,%d), field %d\n",t4,x4,field[t4][x4]);

      if( field[t4][x4] == 0 ) {
        /* Check for the chosen sites in added and removed lists */
        int s = (t2*NX+x2)/2;

        double S1,S2, edS;
        vec_gaussian(psi);
        fM_occupied( chi, psi ); //Now chi is from correct distribution
        S1 = action(psi);

        field[t2][x2] = 0; field[t4][x4] = SOURCE_MONOMER; //make suggested change
        if( cg_MdM_occupied( psi, chi ) == 0 ){
          S2 = action(psi);
        } else {
          S2 = 1e200; //Unable to invert, zero mode
        }
        edS = exp(S1-S2);
        field[t2][x2] = SOURCE_MONOMER; field[t4][x4] = 0; //back

        /* Get the difference in the determinant */
        double det = edS;
        if( mersenne() < det ){
          /* Accepted */
          field[t2][x2] = 0; field[t4][x4] = SOURCE_MONOMER;
          t2 = t4; x2 = x4;
        }
      }

    } //neighbouring site (t3,x3)
   } //steps 
  } //First site site (t2,x2)
 } //attempts

 update_linklists();

 printf("Susceptibility %g \n",(double)steps/(U*4*NDIRS*n_attempts));
  
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

  //measure_propagator();
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

  /* "Warm up" the rng generator */
  for (i=0; i<543210; i++) mersenne();

  printf(" \n++++++++++++++++++++++++++++++++++++++++++\n");
  printf(" 4D free fermion, ( %d , %d ) lattice\n", NT, NX );
  printf(" %d updates per measurements\n", n_measure );
  printf(" m %f \n", m);
  printf(" U %f \n", U);
  printf(" mu %f \n", mu);
  printf(" Random seed %ld\n", seed );


  /* allocate propagator and lists */
  evenlist = malloc( VOLUME*sizeof(int)/2 );
  oddlist  = malloc( VOLUME*sizeof(int)/2 );

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

  for (i=0; i<NX; i++) {
    xup[i] = (i+1) % NX;
    xdn[i] = (i-1+NX) % NX;
  }

  /* fill the staggered eta array */
  for (int t=0; t<NT; t++) for (int x=0; x<NX; x++) {
    eta[t][x][0] = 1;
    if( t%2 == 0 ){
      eta[t][x][1] = 1;
    } else {
      eta[t][x][1] = -1;
    }
  }

  psi =  alloc_field();
  chi =  alloc_field();

  
  /* fill monomers and links */
  for (int t=0; t<NT; t++) for (int x=0; x<NX; x++) {
    field[t][x] = 0;
  }
  update_linklists();

  
  /* and the update/measure loop */
  int changes = 0;
  struct timeval start, end;
  gettimeofday(&start,NULL);
  for (i=1; i<n_loops+1; i++) {

    /* Update */
    changes += update();

    if((i%n_measure)==0){

      /* Time and report */
      gettimeofday(&end,NULL);
      int diff = 1e6*(end.tv_sec-start.tv_sec) + end.tv_usec-start.tv_usec;
      printf("\n%d, %d updates in %.3g seconds, %d successfull changes, %g changes/update\n", i, n_measure, 1e-6*diff,changes,changes/((double) i));
      
      /* Statistics */
      printf("MONOMERS %d \n", n_monomers);
      printf("LINKS %d \n", n_links);

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
















