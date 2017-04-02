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





/* Suggest adding a link
 */
int update_link_at(int s)
{
  int success = 0;
  int t = s/NT, x = s%NT;

  //printf("Testing (%d,%d) %d \n",t,x,field[t][x]);

  if( field[t][x] >= LINK_TUP ) {
    /* Remove link at t,x */
    //Note the factor of 4 from the dirac operator,
    //each dirac link has 0.5
    int dir = field[t][x]-2;
    int t2 = tdir(t,dir), x2 = xdir(x,dir);
    //printf("Removing link at (%d,%d) %d, (%d,%d) \n",t,x,dir,t2,x2);
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

    //printf("Adding link at (%d,%d) %d, (%d,%d) \n",t,x,dir,t2,x2);

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


int add_link_at(int t, int x)
{
  int success = 0;
  //printf("Trying to add at (%d,%d) %d \n",t,x,field[t][x]-2);

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


int remove_link_at(int t, int x)
{
  int success = 0;
  //printf("Trying to remove at (%d,%d) %d \n",t,x,field[t][x]-2);

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
  }
  return success;
}


int linksign(int t,int x,int dir){
    int e = eta[t][x][dir%ND];
    if( dir == TDN ) e *= -1;
    if( dir == XDN ) e *= -1;
    if( x == NX-1 && dir == XUP ) e*= -1;
    if( t == NT-1 && dir == TUP ) e*= -1;
    if( x == 0 && dir == XDN ) e*= -1;
    if( t == 0 && dir == TDN ) e*= -1;
    return e;
}

int find_sign(int t0, int x0, int t, int x){
  int dir = diraclink[t0][x0];
  int sign = linksign(t0,x0,dir);
  //printf("sign (%d,%d) %d %d\n",t0,x0,dir,sign);
  int t1 = tdir(t0, dir), x1= xdir(x0, dir);
  while( t1 != t || x1 != x ){
    dir = diraclink[t1][x1];
    sign *= linksign(t1,x1,dir);
    //printf("sign (%d,%d) %d %d\n",t1,x1,dir,sign);
    t1 = tdir(t1, dir), x1 = xdir(x1, dir);
  }
  return( sign );
}




//Update the Dirac background using a worm update
//Measures the fermion propagator
int update_dirac_background(){
  int propagator[NT];
  for(int t=0;t<NT;t++) propagator[t]=0;

  //Pick a site
  int t= mersenne()*NT, x=mersenne()*NX;
  
  //We break a link to create a fermion correlator.
  //The point the original link points to is the
  //starting point of the correlator
  int dir = diraclink[t][x];
  int t0 = tdir(t, dir), x0 = xdir(x, dir);
  
  double p = 2;    //There is always a factor 0.5 for each link
  if( dir == TUP ) p *= exp(-mu);
  if( dir == TDN ) p *= exp(mu);
  
  if(field[t][x] == 0 && mersenne() < p ){
    //printf("STARTING worm at (%d,%d) %g \n",t,x,p);
    int done = 0;
    diraclink[t][x] = 10;
    
    for(int i=0;;i++) {
      //printf("WORM at (%d,%d) \n",t,x);
      //printf("Sign %d\n",find_sign(t0,x0,t,x));
      //print_config();
      propagator[ (t-t0+NT)%NT ]+=find_sign(t0,x0,t,x);

      //Pick a random direction to create a link
      int dir = mersenne()*NDIRS;

      double p=0.5;
      if( dir == TUP ) p *= exp(mu);
      if( dir == TDN ) p *= exp(-mu);

      //printf(" adding (%d,%d) mu %d, %g \n",t,x,dir,p);

      int t2 = tdir(t, dir), x2 = xdir(x, dir);
      if( t2 == t0 && x2 == x0 ) {
        if( mersenne() < p ){
          diraclink[t][x] = dir;
          //printf("WORM closed\n");
          break;
        }
      } else {
      
       if(mersenne()<0.5){
        int removeddir; //If we accept, this one will be removed
        int t3=t,x3=x;  //New site
        for( int dir2 = 0; dir2<NDIRS; dir2++ ) {
          t3 = tdir(t2, dir2), x3 = xdir(x2, dir2);
          if( diraclink[t3][x3] == opp_dir(dir2) ){
            //Found a site pointing to t2,x2
            removeddir = diraclink[t3][x3];
            break;
          }
        }
        p*=2;
        if( removeddir == TUP ) p *= exp(-mu);
        if( removeddir == TDN ) p *= exp(mu);
      
        //printf(" adding (%d,%d) mu %d, removing (%d,%d) %d %g \n",t,x,dir, t3,x3,removeddir,p);
        
        //flip and follow to a new site with propability p
         if(field[t2][x2]==0) if( mersenne() < p ){
          //flip the link
          diraclink[t][x] = dir;
          diraclink[t3][x3] = 10;
          t=t3; x=x3;
        }
       } else {
    
        // Check if we can add or remove link over the site 
        if( field[t2][x2] == 0 ){
          add_link_at( t2, x2 );
        } else {
          remove_link_at( t2, x2 ); 
        }
       }
      }
    }
  }
  
  for(int t=0;t<NT;t++) printf("PROPAGATOR %d %d\n",t,propagator[t]);


  //printf("Background updated\n");
  return 1;
}




int update()
{
  int changes=0;

  update_dirac_background();
  
  //for(int i=0;i<10;i++) changes+=update_link();
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
      //if(empty==1) {
        if(diraclink[t][x]==10) printf(" . ");
        if(diraclink[t][x]==TUP) printf(" ^ ");
        if(diraclink[t][x]==TDN) printf(" v ");
        if(diraclink[t][x]==XUP) printf(" > ");
        if(diraclink[t][x]==XDN) printf(" < ");
      //}
    }
    printf(" \n");
  }
  printf(" \n");
  //usleep(100000);
}


void measure_charge(){
  int c = 0;
  for (int x=0; x<NX; x++){
    if( diraclink[0][x] == TUP ) c++;
    if( diraclink[1][x] == TDN ) c--;
  }
  printf("CHARGE %d\n",c);
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
    diraclink[t][x] = TUP ;

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
      printf("\n%d, %d updates in %.3g seconds\n", i, n_measure, 1e-6*diff);

      print_config();
      
      /* Statistics */
      printf("MONOMERS %d \n", n_monomers);
      printf("LINKS %d \n", n_links);

      // The measurements can change the configuration
      //without taking the dirac background into account
      for (int t=0; t<NT; t++) for (int x=0; x<NX; x++) field_copy[t][x] = field[t][x];
      
      measure_charge();
      measure_propagator();
      measure_susceptibility();

      n_links = 0;
      for (int t=0; t<NT; t++) for (int x=0; x<NX; x++){
          field[t][x] = field_copy[t][x];
          if(field[t][x] >= LINK_TUP ) n_links++;
      }
      n_links/=2;
      
      gettimeofday(&start,NULL);
    }
  }

  printf(" ** simulation done\n");

  return(1);
}









/*
 This does not take the sign into account
 
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

//Update the Dirac background using a worm update
//Measures the fermion propagator
void move_source_monomer(int *ts0, int *xs0, int dir){

  int t = *ts0, x = *xs0;
  int t2 = tdir(t, dir), x2 = xdir(x, dir);
  int t0,x0;  //point where the propagator starts
  int ts,xs;  //the location of the source
  
  //printf("START move monomer (%d,%d) %d\n",t,x,dir);
  if( field[t][x] != SOURCE_MONOMER ){
    //printf("No source monomer here\n");
    exit(1);
  }
  
  if( field[t2][x2] != 0 ){
    return ;
  }
  
  int done=1;
  
  if( mersenne() > 0.5) {
    //printf("Starting by moving \n");
    int dir2 = find_link_pointing_at(t2,x2);
    int t3 = tdir(t2, dir2), x3 = xdir(x2, dir2);
    int dir3 = find_link_pointing_at(t3,x3);
    int t4 = tdir(t3, dir3), x4 = xdir(x3, dir3);

    double p = 1;
    if( mersenne() < p ){

      field[t][x] = 0;
      field[t3][x3] = SOURCE_MONOMER;
      diraclink[t4][x4] = 10;
      diraclink[t3][x3] = NDIRS;
      diraclink[t][x] = dir;
      
      t0 = t; x0 = x;
      ts = t3; xs = x3;
      t=t4; x=x4;
      done = 0;
    }
    
  } else {
    //printf("Starting by removing a link \n");
    int dir2 = find_link_pointing_at(t2,x2);
    int t3 = tdir(t2, dir2), x3 = xdir(x2, dir2);
    int dir3 = find_link_pointing_at(t3,x3);
    int t4 = tdir(t3, dir3), x4 = xdir(x3, dir3);
    
    t0 = t3, x0 = x3;
    diraclink[t4][x4] = 10;
    ts = t; xs = x;
    t = t4; x = x4;
    done = 0;
    
  }
  
  for(;done==0;){
      //print_config();
      // Now keep updating the worm until it closes
      int dir = mersenne()*NDIRS;
      //printf("move (%d,%d) to dir %d\n",t,x,dir);
      t2 = tdir(t, dir), x2 = xdir(x, dir);
      if( t2==ts && x2 == xs ) {
        //Pointing at the source, move it back and close the worm
        int t3,x3;
        int dir2 = diraclink[t0][x0];
        t3 = tdir(t0, dir2), x3 = xdir(x0, dir2);
        for(dir2 = 0; dir2<NDIRS; dir2++){
            if( tdir(t2,dir2) == t3 && xdir(x2,dir2) == x3 ){
                break;
            }
        }
        
        double p = 1;
        if( mersenne() < p ){
          field[t0][x0] = SOURCE_MONOMER;
          field[t2][x2] = 0;
          diraclink[t][x] = dir;
          diraclink[t0][x0] = NDIRS;
          diraclink[t2][x2] = dir2;
          *ts0 = t0; *xs0 = x0; 
          done = 1;
        }
      } else if(t2 == t0 && x2 == x0) {
        //Ending up at the beginning of the correlator
        //close the worm
        double p = 1;
        if( mersenne() < 0.5 && mersenne() < p ){
          diraclink[t][x] = dir;
          //printf("WORM closed\n");
          *ts0 = ts; *xs0 = xs; 
          done = 1;
        }
      } else {
        //This is the standard worm step
        int removeddir; //If we accept, this one will be removed
        int t3,x3;  //New site
        for( int dir2 = 0; dir2<NDIRS; dir2++ ) {
          t3 = tdir(t2, dir2), x3 = xdir(x2, dir2);
          if( diraclink[t3][x3] == opp_dir(dir2) ){
            //Found a site pointing to t2,x2
            removeddir = diraclink[t3][x3];
            break;
          }
        }
        double p=1;
        if( removeddir == TUP ) p *= exp(-mu);
        if( removeddir == TDN ) p *= exp(mu);
        
        //printf(" adding (%d,%d) mu %d, removing (%d,%d) %d %g \n",t,x,dir, t3,x3,removeddir,p);
        
        //flip and follow to a new site with propability p
         if(field[t2][x2]==0) if( mersenne() < p ){
          //flip the link
          diraclink[t][x] = dir;
          diraclink[t3][x3] = 10;
          t=t3; x=x3;
        }
      }
  }


}

void measure_susceptibility(){
  int nsteps = 0;
  //Pick a site
  int t= mersenne()*NT, x=mersenne()*NX;
  printf("STARTING susceptibility at (%d,%d)\n",t,x);
  if( field[t][x] >= LINK_TUP){
    print_config();
    //Replace the link with two source monomers
    int dir = field[t][x]-2;
    int t2 = tdir(t, dir), x2 = xdir(x, dir);

    field[t][x] = SOURCE_MONOMER;
    field[t2][x2] = SOURCE_MONOMER;
    n_links--;
    
    nsteps++;
    for(;;){
      print_config();
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
         move_source_monomer(&t2,&x2,dir);
       }
    }  
  }

  print_config();
  printf("Susceptibility %g \n",(double)2*nsteps);

}
*/






