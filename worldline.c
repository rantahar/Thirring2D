/******************************************************************
 * Simulate the 1+1d Thirring model using the fermion bag algorithm 
 * (arXiv:0910.5736) using a worldline representation. The mass
 * term is represented as a field of monomers (occupied sites) 
 * and the four fermion term is represented as dimers 
 * (occupied links). 
 *
 ******************************************************************/
#ifdef DEBUG
#include <fenv.h>
#endif

#include "Thirring.h"

/* storage */
int    ***eta;   //Staggered eta matrix
double m;
double U;
double mu;

/* Monomers and links
 * field stores both, 0 for empty, 1 for monomer and 2+dir for links
 */
int n_monomers=0;
int n_links=0;
int **field;

/* Fermion worldline links */
int **diraclink;

/* Size of the fluctuation matrix, used in measurements
 */
int max_changes;

/* Pseudofermion fields */
double *psi,*chi;

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





/* Try to add or remove a link at a given site */
int update_link_at(int s)
{
  int success = 0;
  int t = s%NT, x = s/NT;

  if( field[t][x] >= LINK_TUP ) {
    /* Remove link at t,x */
    //Note the factor of 4 from the dirac operator,
    //each dirac link has 0.5
    int dir = field[t][x]-2;
    int t2 = tdir(t,dir), x2 = xdir(x,dir);
    if( mersenne() < 1/(4*U) ) {
      link_off(t,x,dir);
      /* Replace with 2 opposing arrows */
      diraclink[t][x] = dir;
      diraclink[t2][x2] = opp_dir(dir);
      n_links--;
      success = 1;
    }
  } else {
    /* No link, add if possible */
    int dir = diraclink[t][x];
    int t2 = tdir(t,dir), x2 = xdir(x,dir);
    int nu = diraclink[t2][x2];

    if( dir<NDIRS && dir == opp_dir(nu) ){
      //Two opposing arrows, easy to add
      //Note the factor of 4 from the dirac operator,
      //each dirac link has 0.5
      if( mersenne() < 4*U ) {
        link_on(t,x,dir);
        diraclink[t][x] = NDIRS;
        diraclink[t2][x2] = NDIRS;
        n_links++;
        success = 1;
      }
    }
  }
  return success;
}


int update_link()
{
  return update_link_at( (int) (mersenne()*VOLUME) );
}


/* Try to add a link at given coordinates */
int add_link_at(int t, int x)
{
  int success = 0;

  if( field[t][x] == 0 ) {

    int dir = diraclink[t][x];
    int t2 = tdir(t,dir), x2 = xdir(x,dir);
    int nu = diraclink[t2][x2];

    if( dir<NDIRS && nu == opp_dir(dir) ){
      //Two opposing arrows, easy to add
      //Note the factor of 4 from the dirac operator,
      //each dirac link has 0.5
      if( mersenne() < 4*U ) {
        link_on(t,x,dir);
        diraclink[t][x] = NDIRS;
        diraclink[t2][x2] = NDIRS;
        n_links++;
        success = 1;
      }
    }
  }
  return success;
}


/* Try to remove a link at a given site */
int remove_link_at(int t, int x)
{
  int success = 0;

  if( field[t][x] >= LINK_TUP && field[t][x] <= LINK_XDN ) {
    /* Remove link at t,x */
    //Note the factor of 4 from the dirac operator,
    //each dirac link has 0.5
    int dir = field[t][x]-2;
    int t2 = tdir(t,dir), x2 = xdir(x,dir);
    if( mersenne() < 1./(4*U) ) {
      link_off(t,x,dir);
      /* Replace with 2 opposing arrows */
      diraclink[t][x] = dir;
      diraclink[t2][x2] = opp_dir(dir);
      n_links--;
      success = 1;
    }
  }
  return success;
}


/* Get the sign of a given link */
int linksign( int t, int x, int dir ){
    int e = eta[t][x][dir%ND];
    if( dir == TDN ) e *= -1;
    if( dir == XDN ) e *= -1;
    if( t == NT-1 && dir == TUP ) e*= -1;
    if( t == 0 && dir == TDN ) e*= -1;
#ifdef ANTISYMMETRIC
    if( x == NX-1 && dir == XUP ) e*= -1;
    if( x == 0 && dir == XDN ) e*= -1;
#endif
    return e;
}

/* Calculate the sign of the propagator between the head and tail of a worm,
   the fermion propagator. Just follow the worm and count the signs. */
int find_sign(int t0, int x0, int t, int x){
  int dir = diraclink[t0][x0];
  int sign = linksign(t0,x0,dir);
  int t1 = tdir(t0, dir), x1= xdir(x0, dir);
  while( t1 != t || x1 != x ){
    dir = diraclink[t1][x1];
    sign *= linksign(t1,x1,dir);
    t1 = tdir(t1, dir), x1 = xdir(x1, dir);
  }
  return( sign );
}




//Update the Dirac background using a worm update
//Can also measure the fermion propagator
int update_dirac_background(){
  int propagator[NT];
  for(int t=0;t<NT;t++) propagator[t]=0;

  //Pick a site
  int t= mersenne()*NT, x=mersenne()*NX;

  //We break a link to create a fermion correlator.
  //The point the original link points to is the
  //starting point of the correlator
  int dir = diraclink[t][x];
  
  double p = 2;    //There is always a factor 0.5 for each link
  if( dir == TUP ) p *= exp(-mu);
  if( dir == TDN ) p *= exp(mu);
  
  if(field[t][x] == 0 && mersenne() < p ){
    int t0 = tdir(t, dir), x0 = xdir(x, dir);
    int done = 0;
    diraclink[t][x] = 10;
    
    for(int i=0;;i++) {
      //Also calculate the propagator between the two defects
      //This is the fermion propagator
      propagator[ (t-t0+NT)%NT ]+=find_sign(t0,x0,t,x);

      //Pick a random direction to create a link
      int dir = mersenne()*NDIRS;

      double p=0.5;
      if( dir == TUP ) p *= exp(mu);
      if( dir == TDN ) p *= exp(-mu);

      // Check if the worm closes
      int t2 = tdir(t, dir), x2 = xdir(x, dir);
      if( t2 == t0 && x2 == x0 ) {
        if( mersenne() < p ){
          diraclink[t][x] = dir;
          break;
        }
      } else {
        
        // Branch between adding/removing a link and propagating the worm
        if(mersenne()<0.5){
        
         // If the middle site is empty, we can point a link to it
         // There will then be a second link pointing at it, which we need
         // to remove. That will become the head of the worm
         if(field[t2][x2]==0) {
          int removeddir;
          int t3,x3;
          for( int dir2 = 0; dir2<NDIRS; dir2++ ) {
            t3 = tdir(t2, dir2), x3 = xdir(x2, dir2);
            if( diraclink[t3][x3] == opp_dir(dir2) ){
              //Found the site pointing to t2,x2
              removeddir = diraclink[t3][x3];
              break;
            }
          }

          //Calculate the propability and flip the link
          p*=2;
          if( removeddir == TUP ) p *= exp(-mu);
          if( removeddir == TDN ) p *= exp(mu);
          if( mersenne() < p ){
            //flip the link
            diraclink[t][x] = dir;
            diraclink[t3][x3] = 10;
            t=t3; x=x3;
          }
         }

        } else {
    
          // Check if we can add or remove link at the site
          // This is the only place in the worm where dimers
          // are actually updated
          if( field[t2][x2] == 0 ){
            add_link_at( t2, x2 );
          } else {
            remove_link_at( t2, x2 ); 
          }
       }
      }
    }
  }

  //for(int t=0;t<NT;t++) printf("PROPAGATOR %d %d\n",t,propagator[t]);
  return 1;
}



/* A full update function. A single worm update followed by a number of random
   link updates */
int update()
{
  int changes=0;

  update_dirac_background();

  int n_update = NX;  
  for(int i=0;i<n_update;i++) update_link();
  return changes;
}











void print_config()
{
  for (int t=NT; t--; ) {
    for (int x=0; x<NX; x++){
      int empty = 1;
      if(field[t][x]==MONOMER)  { empty = 0; printf(" o "); }
      if(field[t][x]==LINK_TUP) { empty = 0; printf(" | "); }
      if(field[t][x]==LINK_XUP) { empty = 0; printf(" --"); }
      if(field[t][x]==LINK_TDN) { empty = 0; printf(" | "); }
      if(field[t][x]==LINK_XDN) { empty = 0; printf("-- "); }
      if(field[t][x]==SOURCE_MONOMER) { empty = 0; printf(" x "); } //A source monomer
      if(diraclink[t][x]==10) printf(" . ");
      if(diraclink[t][x]==TUP) printf(" ^ ");
      if(diraclink[t][x]==TDN) printf(" v ");
      if(diraclink[t][x]==XUP) printf(" > ");
      if(diraclink[t][x]==XDN) printf(" < ");
    }
    printf(" \n");
  }
  printf(" \n");
}


/* Measure the fermion and chiral charges */
void measure_charge(int *c, int *q){
  int _c = 0;
  int _q = 0;
  for (int x=0; x<NX; x++){
    if( diraclink[0][x] == TUP ){
      _c+= 1;
      _q+= ((x)%2==1?1:-1); 
    }
    if( diraclink[1][x] == TDN ){
      _c-= 1;
      _q+= ((x)%2==1?1:-1);
    }
    if( field[0][x] == LINK_TUP ){
      _q+= 2*((x)%2==1?1:-1);
    }
  }
  *c = _c;
  *q = _q;
}








//In the precense of (source) monomers, the sign
//can be negative. Assuming two monomers, find 
//the sign of the configuration
//Essentially count the surfaces crossed
//on a path between the monomers
int sign_with_monomers(int t0, int x0, int t1, int x1){
  int charge = 0;
  int ds, t = t0, x = x0;
  //First travel in t-direction
  if(t0<t1) {
    for(;t<t1-1;){
      t++;
      if(diraclink[t][x]==XUP) charge++;
      if(diraclink[t][xup[x]]==XDN) charge++;
    }
  }
  else {
    for(;t>t1;t--){
      if(diraclink[t][x]==XUP) charge++;
      if(diraclink[t][xup[x]]==XDN) charge++;
    }
  }

  //Then in the x direction
  if(x0<x1) {
    for(;x<x1-1;){
      x++;
      if(diraclink[t][x]==TUP) charge++;
      if(diraclink[tup[t]][x]==TDN) charge++;
    }
  }
  else {
    for(;x>x1;x--){
      if(diraclink[t][x]==TUP) charge++;
      if(diraclink[tup[t]][x]==TDN) charge++;
    }
  }
  return charge%2 == 0 ? 1:-1;
}


/* Find a link pointing at a given site */
int find_link_pointing_at( int t, int x ){
  int t2,x2, dir;
  for( dir = 0; dir<NDIRS; dir++ ) {
    t2 = tdir(t, dir), x2 = xdir(x, dir);
    if( diraclink[t2][x2] == opp_dir(dir) ){
      //Found a site pointing to t,2
      break;
    }
  }
  return dir;
}


/* A worm for measuring the susceptibility in a wordline background
 * does not take the signs into account */

/* Propagate the worm by moving one of the sources to an unoccupied site.
 * This will always trigger a worldline worm update.
 * For detailed balance the worm can start in two ways:
 *   Moving the source and creating a defect in the background
 *   Just creating starting a worm update next to the source
 * Equivalently it can close in two way
 *   Hitting the source and moving it to close the worm
 *   Closing by hitting the defect
 */
int move_source_monomer(int *ts0, int *xs0, int dir){

  int t = *ts0, x = *xs0;
  int t2 = tdir(t, dir), x2 = xdir(x, dir);
  int t0,x0;  //point where the propagator starts
  int ts,xs;  //the location of the other end of the worm (source)
  
  // Check correct use
  if( field[t][x] != SOURCE_MONOMER ){
    exit(1);
  }
  
  // This algorithm is only used when the neighbouring site
  // is unoccupied
  if( field[t2][x2] != 0 ){
    return 0;
  }
  
  // Start by finding the neighboring sites
  int done=1;
  int dir2 = find_link_pointing_at(t2,x2);
  int t3 = tdir(t2, dir2), x3 = xdir(x2, dir2);
  int dir3 = find_link_pointing_at(t3,x3);
  int t4 = tdir(t3, dir3), x4 = xdir(x3, dir3);


  // Branch between moving the source by two sites and keeping it at place
  // (the second one is required for detailed balance)
  if( mersenne() > 0.5) {

    double p = 2;
    if( dir == TUP ) p *= exp(mu);
    if( dir == TDN ) p *= exp(-mu);
    if( diraclink[t3][x3] == TUP ) p *= exp(-mu);
    if( diraclink[t3][x3] == TDN ) p *= exp(mu);
    if( diraclink[t4][x4] == TUP ) p *= exp(-mu);
    if( diraclink[t4][x4] == TDN ) p *= exp(mu);

    if( mersenne() < p ){
      // Remove the source monomer at the original point and move it to point 3
      field[t][x] = 0;
      field[t3][x3] = SOURCE_MONOMER;

      // This also means cutting the Worldline link at the site
      diraclink[t4][x4] = 10;
      diraclink[t3][x3] = NDIRS;
      diraclink[t][x] = dir; //this now points to site 2
      
      // Remember the new sites
      t0 = t; x0 = x;
      ts = t3; xs = x3;
      t=t4; x=x4;
      done = 0;
    }
    
  } else {

    // Don't move the source, just start a worm
    double p = 2;
    if( diraclink[t4][x4] == TUP ) p *= exp(-mu);
    if( diraclink[t4][x4] == TDN ) p *= exp(mu);

    if( mersenne() < p ){
      t0 = t3, x0 = x3;
      diraclink[t4][x4] = 10;
      ts = t; xs = x;
      t = t4; x = x4;
      done = 0;
    }
    
  }
  
  for(;done==0;){
      // Now keep updating the worm until it closes
      int dir = mersenne()*NDIRS;
      t2 = tdir(t, dir), x2 = xdir(x, dir);

      //If pointing at the source, move it back and close the worm
      if( t2==ts && x2 == xs ) {

        // Find the relevant directions
        int t3,x3;
        int dir2 = diraclink[t0][x0];
        t3 = tdir(t0, dir2), x3 = xdir(x0, dir2);
        for(dir2 = 0; dir2<NDIRS; dir2++){
            if( tdir(t2,dir2) == t3 && xdir(x2,dir2) == x3 ){
                break;
            }
        }
      
       // Calculate the propability of moving the source
       double p = 0.5;
       int dir0 = diraclink[t0][x0];
       if( dir == TUP ) p *= exp(mu);
       if( dir == TDN ) p *= exp(-mu);
       if( dir2 == TUP ) p *= exp(mu);
       if( dir2 == TDN ) p *= exp(-mu);
       if( diraclink[t0][x0] == TUP ) p *= exp(-mu);
       if( diraclink[t0][x0] == TDN ) p *= exp(mu);

       // Move and close the worm
       if( mersenne() < p) {
          field[t0][x0] = SOURCE_MONOMER;
          field[t2][x2] = 0;
          diraclink[t][x] = dir;
          diraclink[t0][x0] = NDIRS;
          diraclink[t2][x2] = dir2;
          *ts0 = t0; *xs0 = x0; 
          done = 1;
        }

      } else if(t2 == t0 && x2 == x0) {

        //If pointing at the defect, try to close the worm
        double p = 0.5;
        if( dir == TUP ) p *= exp(mu);
        if( dir == TDN ) p *= exp(-mu);
        if( mersenne() < p ){
          diraclink[t][x] = dir;
          *ts0 = ts; *xs0 = xs; 
          done = 1;
        }

      } else if(field[t2][x2]==0) {

        // This is the standard worm step
        // Find the directions and the propability of the move
        int removeddir;
        int t3,x3; 
        for( int dir2 = 0; dir2<NDIRS; dir2++ ) {
          t3 = tdir(t2, dir2), x3 = xdir(x2, dir2);
          if( diraclink[t3][x3] == opp_dir(dir2) ){
            break;
          }
        }
        double p=1;
        if( dir == TUP ) p *= exp(mu);
        if( dir == TDN ) p *= exp(-mu);
        if( diraclink[t3][x3] == TUP ) p *= exp(-mu);
        if( diraclink[t3][x3] == TDN ) p *= exp(mu);
        
        //flip and follow to a new site with propability p
        if( mersenne() < p ){
          //flip the link
          diraclink[t][x] = dir;
          diraclink[t3][x3] = 10;
          t=t3; x=x3;
        }
      }
  }
}



/* The actual worm for susceptibility. Does not take signs into account */
double measure_susceptibility_with_background( ){
  #ifdef OPENX
  int nsteps = 0;
  int n_attempts = 1;
  double scale_factor = n_links/(2.*VOLUME*n_attempts);
  
  for(int i=0; i<n_attempts;i++){
    //Pick a random link
    int t=0,x=0;
    if(n_links>0) do{
      t= mersenne()*NT, x=mersenne()*NX;
    } while(field[t][x] < LINK_TUP);
  
    if( field[t][x] >= LINK_TUP){
      //Replace the link with two source monomers
      int dir = field[t][x]-2;
      int t2 = tdir(t, dir), x2 = xdir(x, dir);

      field[t][x] = SOURCE_MONOMER;
      field[t2][x2] = SOURCE_MONOMER;
      n_links--;
    
      //Propagate the sources untill the worm closes
      for(;;){
        nsteps+=sign_with_monomers(t,x,t2,x2);
        int dir = mersenne()*NDIRS;
        int t3 = tdir(t2, dir), x3 = xdir(x2, dir);
      
        if( t3==t && x3==x ) {
           // Back at the original site, turn into a link 
           field[t][x]=0; field[t2][x2]=0;
           link_on(t2,x2,dir);
           n_links++;
           break;

         } else if( field[t3][x3] >= LINK_TUP ) {
           // found a link, flip it  
           int linkdir = field[t3][x3] - LINK_TUP;
           int t4 = tdir(t3,linkdir), x4 = xdir(x3,linkdir);

           field[t2][x2] = 0;
           link_off(t3,x3,linkdir); link_on(t2,x2,dir);
           field[t4][x4] = SOURCE_MONOMER;

           t2 = t4; x2 = x4;

         } else {
           // No link at found, run a worm that may move the
           // Source by two sites by updating the background field
           move_source_monomer(&t2,&x2,dir);
         }
      }
    }
  }

  return (double)nsteps*scale_factor;
  #else
  return 0;
  #endif
}









/* Figure out the sign of the current
   configuration. In the most general
   case we go over all fermion loops. */
int find_config_sign(){
  int sign = 1;

  // Mark the checked sites to avoid double counting
  int checked[NT][NX];
  for(int t=0; t<NT; t++) for(int x=0;x<NX;x++)
    checked[t][x] = 0;

  // Cycle trough all the sites and follow around each loop
  for(int t=0; t<NT; t++) for(int x=0;x<NX;x++) if(field[t][x]==0) if(checked[t][x] == 0){
    int dir = diraclink[t][x];
    int loop_sign = linksign(t,x,dir);
    checked[t][x] = 1;
    int t1 = tdir(t, dir), x1= xdir(x, dir);
    while( t1 != t || x1 != x ){
      dir = diraclink[t1][x1];
      loop_sign *= linksign(t1,x1,dir);
      checked[t1][x1] = 1;
      t1 = tdir(t1, dir), x1 = xdir(x1, dir);
    }
    //loop closed
    sign *=-loop_sign;
  }
  return sign;
}





/* Main function
 */
int main(int argc, char* argv[])
{
  #ifdef DEBUG
  feenableexcept(FE_INVALID | FE_OVERFLOW);
  #endif 

  int i,n_loops,n_measure,n_average;
  long seed;

  /* Read in the input */
  printf(" Number of updates : ");
  scanf("%d",&n_loops);

  printf(" Updates / measurement : ");
  scanf("%d",&n_measure);
  
  printf(" Average over : ");
  scanf("%d",&n_average);
  printf("%d\n",n_average);

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

  /* "Warm up" the rng generator */
  for (i=0; i<543210; i++) mersenne();

  printf(" \n++++++++++++++++++++++++++++++++++++++++++\n");
  printf(" 2D Thirring model, ( %d , %d ) lattice\n", NT, NX );
  printf(" %d updates per measurements\n", n_measure );
  printf(" m %f \n", m);
  printf(" U %f \n", U);
  printf(" mu %f \n", mu);
  printf(" Size of fluctuation matrix %d\n", max_changes );
  printf(" Random seed %ld\n", seed );

  /* Allocate location and field arrays */
  field = malloc( NT*sizeof(int *) );
  diraclink = malloc( NT*sizeof(int *) );
  eta = malloc( NT*sizeof(int *) );
  for (int t=0; t<NT; t++){
    field[t] = malloc( (NX+1)*sizeof(int) );
    diraclink[t] = malloc( (NX+1)*sizeof(int) );
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

  /* fill monomers and links */
  for (int t=0; t<NT; t++) for (int x=0; x<NX; x++) {
    field[t][x] = 0;
    diraclink[t][x] = XUP + 2*(x%2) ;

  }

  int ** field_copy = malloc( NT*sizeof(int *) );
  for (int t=0; t<NT; t++) field_copy[t] = malloc( (NX+1)*sizeof(int) );
  
#ifdef OPENX
  for (int t=0; t<NT; t++) {
    field[t][NX] = EMPTY; //Site doesn't exist, no links or monomers, but not free either
    field_copy[t][NX] = EMPTY;
    diraclink[t][NX] = EMPTY;
  }
#endif
  
  /* and the update/measure loop */
  int sum_monomers = 0;
  int sum_links = 0;
  int sum_charge = 0;
  int sum_c2 = 0;
  int sum_q = 0;
  int sum_q2 = 0;
  double sum_susc = 0;
  double sum_susc_wb = 0;
  int sum_sign=0;

  struct timeval start, end;
  double updatetime=0, measuretime = 0;
  gettimeofday(&start,NULL);
  for (i=1; i<n_loops+1; i++) {

    /* Update */
    update();

    if((i%n_measure)==0){

      /* Time and report */
      gettimeofday(&end,NULL);
      updatetime += 1e6*(end.tv_sec-start.tv_sec) + end.tv_usec-start.tv_usec;

      gettimeofday(&start,NULL);

      int sign = find_config_sign();
      sum_sign += sign;
      //measure_propagator(); //This includes an invertion and therefore takes time

      sum_monomers += sign*n_monomers;
      sum_links += sign*n_links;

      int c, q;
      measure_charge(&c, &q);
      sum_charge += sign*c;
      sum_c2 += sign*c*c;
      sum_q += sign*q;
      sum_q2 += sign*q*q; 

      sum_susc_wb += sign*measure_susceptibility_with_background();
      gettimeofday(&end,NULL);
      measuretime += 1e6*(end.tv_sec-start.tv_sec) + end.tv_usec-start.tv_usec;

      if((i%(n_average*n_measure))==0){
        printf("\n%d, %d updates in %.3g seconds\n", i, n_average*n_measure, 1e-6*updatetime);
        printf("%d, %d measurements in %.3g seconds\n", i, n_average, 1e-6*measuretime);
        updatetime = 0; measuretime = 0;

        printf("MONOMERS %g \n", (double)sum_monomers/n_average);
        printf("LINKS %g \n", (double)sum_links/n_average);
        printf("CHARGE %g %g \n", (double)sum_charge/n_average, (double)sum_c2/n_average);
        printf("QCHI %g %g \n", (double)sum_q/n_average, (double)sum_q2/n_average);
        //printf("Susc %g \n", (double)sum_susc/n_average);
        printf("SUSCEPTIBILITY %g \n", (double)sum_susc_wb/n_average);
        printf("SIGN %g\n", (double)sum_sign/n_average);

        sum_monomers = 0; sum_links = 0; sum_charge = 0;
        sum_c2 = 0; sum_q = 0; sum_q2 = 0; sum_susc = 0;
        sum_susc_wb = 0; sum_sign = 0;
      }
      
      gettimeofday(&start,NULL);
    }
  }

  printf(" ** simulation done\n");

  return(1);
}














