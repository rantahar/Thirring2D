/****************************************************************
 * Simulate the 1+1D Thirring model using the fermion bag algorithm 
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
extern int    *unoccupied_evenlist, *unoccupied_oddlist;  //Lists of occupied sites
extern int *added_evensites, *added_oddsites, *removed_evenlist, *removed_oddlist;
extern int n_added_even,n_added_odd,n_removed_even,n_removed_odd;
double m, m_update;
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

/* Inverse of the propagator matrix for the background config
 */
double *Ginv;


/* Neighbour index arrays, to be filled at the beginning
 */
int *tup,*xup,*tdn,*xdn;

int current_sign=1;

double det_added_link(int t, int x, int t2, int x2);
double det_removed_link(int t, int x, int t2, int x2 );
double det_moved_monomer(int t, int x, int t2, int x2);
void new_link(int t, int x, int nu);
void new_monomer(int t1, int x1, int t2, int x2);
void removed_link(int t, int x, int nu);
void removed_monomer(int t1, int x1, int nu);
void moved_source_monomer(int t2, int x2, int t4, int x4);
void not_moved_source_monomer(int t2, int x2, int t4, int x4);
#ifdef FLUCTUATION_MATRIX
void update_background( );
#endif
double determinant();
double determinant_mu( double mu);



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




/* Update the lists oflinks and monomers
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

    if( mersenne() < 0.5 ) {
      if( m_update > 0 ){
        double d = det_added_link(t,x,t2,x2);
        double p=( n_legalemptylinks/(double)n_legal_monomers_after_adding(t,x,t2,x2) );
        p*= m_update*m_update*d;
        if( mersenne() < p ) {
          new_monomer(t,x,t2,x2);
          success = 1;
        }
      }
    } else {
      double d = det_added_link(t,x,t2,x2);
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
      p *= M/(m_update*m_update);
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



/* A full update consists of a number of add and remove attempts */
int update()
{
  int changes=0;
  update_linklists();
  if(mersenne()>0.5){
    changes+=remove_link();
  } else {
    changes+=add_link();
  }
  
  #ifdef DEBUG
  //check_det( );
  //print_config();
  #endif
  //move_monomer();
  return changes;
}




/* Functions for measurements
 */

/* Propagator
 */
void measure_propagator(){ 
  double prop[NT];
  double j[NT];
  double c[NT];
  for( int t1=0; t1<NT; t1++){
    prop[t1]=0; j[t1] = 0; c[t1]=0;
  }

  double **source, **propagator;
  source = alloc_vector();
  propagator = alloc_vector();

  for( int t1=0;t1<2;t1++) 
  for( int x1=0;x1<NX;x1++) {
    int bc_dn = 1, bc_up = 1;
    if(t1 == 0) bc_dn = -1;  if(t1 == NT-1) bc_up = -1;

    if( field[t1][x1] == 0 ) {
      vec_zero( source );
      source[t1][x1] = 1;

      vec_zero( propagator );
      cg_propagator(propagator,source);

      j[t1] += bc_up*exp(mu)*eta[t1][x1][0]*propagator[tup[t1]][x1];
      c[t1] += bc_up*((x1+t1)%2 ==0 ? 1:-1 )* exp(mu)*eta[t1][x1][0]*propagator[tup[t1]][x1];

      j[tdn[t1]] += bc_dn*exp(-mu)* eta[tdn[t1]][x1][0]*propagator[tdn[t1]][x1];
      c[tdn[t1]] -= bc_dn*((x1+tdn[t1])%2 ==0 ? 1:-1 )* exp(-mu)* eta[tdn[t1]][x1][0]*propagator[tdn[t1]][x1];
      
    }

    if( field[t1][x1]== LINK_TUP ) { 
      double cl = ((t1+x1)%2==0 ? 1:-1 )*2;
      c[t1] += 2*cl;
    }
  }

  
  for( int t2=0; t2<NT; t2++) printf("Propagator %d %g\n", t2, current_sign*prop[t2]/(VOLUME) );
  for( int t1=0;t1<1;t1++) printf("Charge %d %g\n", t1, current_sign*j[t1]/2 );
  for( int t1=0;t1<1;t1++)  printf("Qchi %d %g\n", t1, current_sign*c[t1]/2 );
  printf("Qchi2  %g\n", current_sign*c[0]*c[0]/4 );

  free_vector(source);
  free_vector(propagator);
}


extern double accepted_det;
extern double fluctuation_det;
extern int bc_sign;
extern int fluctuation_sign;
extern int moved_old_site, moved_new_site;

/* Measure the susceptibility using a worm algorithm. Introduce 2 source monomers
 * with the weigth J (=U*NDIRS/V). Allow one to move around (produce configurations)
 * until it contacts with the the other one, remove with the appropriate weight.
 * The number of intermediate configurations counts the Z_J/Z_0, which is the 
 * susceptibility , Z_J/Z_0 = U*NDIRS/V * dZ_J/dJ |_J=0.
 */
void measure_susceptibility(){
 int n = n_bc_monomers/2 + n_bc_links;
 int steps = 0;
 int n_attempts=20;
 double scale_factor = (double)n_links/(NDIRS*n_attempts*VOLUME);

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
   
   n_monomers += 2;
   n_links -= 1;
   steps+=current_sign;

   #ifdef DEBUG
   //print_config();
   check_det();
   #endif
   
   for(;;){
     /* Now we are at (t2,x2), and the link is off. Try to move. */
     int dir = NDIRS*mersenne();
     int t3 = tdir(t2,dir), x3 = xdir(x2,dir);

     if( mersenne() > 0.5 ) {
       /* Try to exchange with an occupied neighbor */
       if( t3==t1 && x3==x1 ) {
         /* Back at the original site, turn into a link */
         field[t1][x1]=0; field[t2][x2]=0;
         link_on(t2,x2,dir);
         n_monomers -= 2;
         n_links += 1;
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
     } else if( field[t3][x3]!=EMPTY ) {
#ifdef WITH_MASS_MONOMERS
       /*  Try to hop over */
       int dir2; 
       do dir2 = NDIRS*mersenne();
       while( dir2==opp_dir(dir));
       int t4 = tdir(t3,dir2), x4 = xdir(x3,dir2);

       if( mersenne() < 0 ){
         if( field[t4][x4] == 0 ) {
           /* Get the difference in the determinant */
           double det = det_moved_monomer( t2, x2, t4, x4 );
           if( mersenne() < det ){
             /* Accepted */
             moved_source_monomer(t2,x2,t4,x4);
             t2 = t4; x2 = x4;
           } else {
             /* Rejected */
             not_moved_source_monomer(t2,x2,t4,x4);
           }
         }
       } else {
         if( field[t3][x3]==0 && field[t4][x4]==0 ) {
           double p = U*det_added_link(t3,x3,t4,x4);
           if( mersenne() < p ){
             new_link(t3,x3,dir2);
             field[t2][x2] = 0;
             link_off(t3,x3,dir2); link_on(t2,x2,dir);
             field[t4][x4] = SOURCE_MONOMER;

             t2 = t4; x2 = x4;
           }
         } else if(field[t3][x3] == (LINK_TUP + dir2) ) {
           double p = det_removed_link(t2,x2,t3,x3)/U;
           if( mersenne() < p ){
             
             field[t2][x2] = 0;
             link_off(t3,x3,dir2); link_on(t2,x2,dir);
             field[t4][x4] = SOURCE_MONOMER;
             removed_link(t2,x2,dir);

             t2 = t4; x2 = x4;
           }
         }
       }
#else
      /* Hop to the neighbour */
      if( field[t3][x3] == 0 ) {
       /* Get the difference in the determinant */
        double det;
        det = fabs(det_moved_monomer( t2, x2, t3, x3 ));

        if( mersenne() < det ){
          /* Accepted */
          moved_source_monomer(t2,x2,t3,x3);
          t2 = t3; x2 = x3;
        } else {
          /* Rejected */
          not_moved_source_monomer(t2,x2,t3,x3);
        }
      }
#endif
    } //neighbouring site (t3,x3)

    steps+=current_sign;
   } //steps
   update_linklists();
 } //attempts


 printf("Susceptibility %g \n",(double)steps*scale_factor);
  
}




void print_config()
{
  for (int t=0; t<NT; t++) {
    for (int x=0; x<NX; x++){
      int empty = 1;
      if(field[t][x]==MONOMER)  { empty = 0; printf(" o "); }
      if(field[t][x]==LINK_TUP) { empty = 0; printf(" v "); }
      if(field[t][x]==LINK_XUP) { empty = 0; printf("  >"); }
      if(field[t][x]==LINK_TDN) { empty = 0; printf(" ^ "); }
      if(field[t][x]==LINK_XDN) { empty = 0; printf("<  "); }
      if(field[t][x]==SOURCE_MONOMER) { empty = 0; printf(" x "); } //A source monomer
      if(empty==1) { printf(" . "); }
    }
    printf(" \n");
  }
  printf(" \n");
}


static int measurement = 0;
void measure()
{ 
#ifdef DEBUG
  //print_config();
#endif

  measure_propagator( current_sign );
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
  printf(" 2D free fermion, ( %d , %d ) lattice\n", NT, NX );
  printf(" %d updates per measurements\n", n_measure );
  printf(" m %f \n", m);
  printf(" U %f \n", U);
  printf(" mu %f \n", mu);
  printf(" Size of fluctuation matrix %d\n", max_changes );
  printf(" Random seed %ld\n", seed );
  printf(" Start configuration %s\n", start_config );

  /* Set monomer and link fields and counters */
#ifdef WITH_MASS_MONOMERS
  m_update = m;
  m=0;
#else
  m_update = 0;
#endif

  field = malloc( NT*sizeof(int *) );
  eta = malloc( NT*sizeof(int *) );
  for (int t=0; t<NT; t++){
    field[t] = malloc( (NX+1)*sizeof(int) );
    eta[t] = malloc( (NX+1)*sizeof(int *) );
    for (int x=0; x<NX+1; x++) {
     eta[t][x] = malloc( 2*sizeof(int) );
    }
  }

  // Set the neighbour arrays
  xup = malloc( (NX+1)*sizeof(int) );
  xdn = malloc( (NX+1)*sizeof(int) );
  tup = malloc( NT*sizeof(int) );
  tdn = malloc( NT*sizeof(int) );

  for (int t=0; t<NT; t++) for (int x=0; x<NX; x++) {
    field[t][x] = 0;
  }
#ifdef OPENX
  for (int t=0; t<NT; t++) {
    field[t][NX] = EMPTY; //Site doesn't exist, no links or monomers, but not free either
  }
#endif


  /* allocate propagator and lists */
#ifdef FLUCTUATION_MATRIX
  Dinv = malloc( VOLUME*VOLUME*sizeof(double) );
  Ginv = malloc( VOLUME*VOLUME*sizeof(double) );
  evenlist = malloc( VOLUME*sizeof(int) );
  oddlist  = malloc( VOLUME*sizeof(int) );
  unoccupied_evenlist = malloc( VOLUME*sizeof(int) );
  unoccupied_oddlist  = malloc( VOLUME*sizeof(int) );

  added_evensites = malloc(VOLUME*sizeof(int));
  added_oddsites = malloc(VOLUME*sizeof(int));
  removed_evenlist = malloc(VOLUME*sizeof(int));
  removed_oddlist = malloc(VOLUME*sizeof(int));
  
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
#endif

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
#ifdef PROPAGATOR_MATRIX
  calc_Dinv( );
#endif

  if( strcmp(start_config,"cold")!=0 ) {
    printf(" Reading configuration file\n" );
    read_config(start_config);
  } else {
    printf(" Starting from a cold configuration\n" );
  }

  /* Background and fluctuation matrix */
  update_linklists();
#ifdef FLUCTUATION_MATRIX
  update_background( );
#endif

  /* and the update/measure loop */
  int changes = 0;
  struct timeval start, end;
  gettimeofday(&start,NULL);
  for (i=1; i<n_loops+1; i++) {

    /* Update */
    changes += update();
    
    if((i%n_measure)==0){

      /* Update the background and count links and monomers */
      update_linklists();
      
      /* Time and report */
      gettimeofday(&end,NULL);
      int diff = 1e6*(end.tv_sec-start.tv_sec) + end.tv_usec-start.tv_usec;
      printf("\n%d, %d updates in %.3g seconds, %d successfull changes, %g changes/update\n", i, n_measure, 1e-6*diff,changes,changes/((double) i));

      /* Statistics */
      printf("MONOMERS %d \n", current_sign*n_monomers);
      printf("LINKS %d \n", current_sign*n_links);
      printf("Determinant base %g \n", accepted_det);
      #ifdef DEBUG
      print_config();
      #endif DEBUG
      printf("Sign %d \n", current_sign);

      /* Write configuration */
      write_config();

      gettimeofday(&start,NULL);

      if(i>n_thermalize) {
        /* Do measurements */
        double tdet = determinant_mu(0), tdet2;
        printf("Determinant at mu=0 %g \n", tdet);
        printf("Determinant ratio %g \n", accepted_det/tdet);
        int sector = 0;
        for(double tmu=0;tmu<=mu;tmu+=SECTOR_STEP){
          tdet2 = determinant_mu(tmu);
          if( tdet*tdet2 < 0 ){
            sector++;
          }
          tdet = tdet2;
        }
        printf("Sector %d\n",sector);

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



























