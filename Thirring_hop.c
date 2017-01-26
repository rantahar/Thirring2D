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
double m;
double U;
double mu;

/* Monomers and links
 * field stores both, 0 for empty, 1 for monomer and 2+dir for links
 */
int n_monomers=0;
int n_links=0;
int field[NT][NX];
int diraclink[NT][NX];

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





/* Suggest adding a link
 */
int update_link()
{
  int success = 0;

  /* Draw a legal site */
  int s = (int) (mersenne()*VOLUME);
  int t = s%NX, x = s/NX;
  
  if( field[t][x] > MONOMER ) {
    /* Remove link at t,x */
    //Note the factor of 4 from the dirac operator,
    //each dirac link has 0.5
    if( mersenne() < 1/(4*U) ) {
      int mu = field[t][x]-2;
      int t2 = tdir(t,mu), x2 = xdir(x,mu);
      link_off(t,x,mu);
      /* Replace with 2 opposing arrows */
      diraclink[t][x] = mu;
      diraclink[t2][x2] = opp_dir(mu);
      n_links--;
      success = 1;
    }
  } else {
    /* No link, add if possible */
    int mu = diraclink[t][x];
    int t2 = tdir(t,mu), x2 = xdir(x,mu);
    int nu = diraclink[t2][x2];

    //printf("Adding link at %d %d %d\n",t,x,mu);

    if( mu<NDIRS && mu == opp_dir(nu) ){  //Two opposing arrows, easy to add
      //Note the factor of 4 from the dirac operator,
      //each dirac link has 0.5
      if( mersenne() < 4*U ) {
        link_on(t,x,mu);
        diraclink[t][x] = NDIRS;
        diraclink[t2][x2] = NDIRS;
        n_links++;
        success = 1;
        //printf("Accepted");
      }
    }
  }
  
  return success;
}




//Update the Dirac background using a worm update
int update_dirac_background(){

  //Pick a site
  int t= mersenne()*NT, x=mersenne()*NT;
  int done = 0;

  if(field[t][x] == 0 ) do {
    //Pick a random direction to turn the link
    int mu = mersenne()*NDIRS;

    //printf(" t %d x %d mu0 %d mu %d field %d \n",t,x,diraclink[t][x],mu, field[t][x]);

    //follow to a new site
    int t2 = tdir(t, mu), x2 = xdir(x, mu);
    if( field[t2][x2] == 0 ) { //Cannot point at a monomer

      //flip the link
      diraclink[t][x] = mu;

      //See if we now have two sites pointing to the same place
      done = 1;
      int t3,x3; //New site
      for( int nu = 0; nu<NDIRS; nu++ ) {
        t3 = tdir(t2, nu), x3 = xdir(x2, nu);
        if( diraclink[t3][x3] == opp_dir(nu) ){
          //Found a site pointing to t2,x2
          if( !( t3==t && x3==x) ){
            //And it's a new one, flip it next
            done = 0;
            break;
          }
        }
      }

      //move to the new site
      t = t3; x = x3;
      //printf(" site found at %d %d \n",t3,x3);
    }
  } while( done == 0 );  //Quit if no new site is found
  

  //printf("Background updated\n");
  return 1;
}




int update()
{
  int changes=0;

  update_dirac_background();

  changes += update_link();
  
  return changes;
}

















/* Measure the susceptibility using a worm algorithm. Introduce 2 source monomers
 * with the weigth J (=U*NDIRS/V). Allow one to move around (produce configurations)
 * until it contacts with the the other one, remove with the appropriate weight.
 * The number of intermediate configurations counts the Z_J/Z_0, which is the 
 * susceptibility , Z_J/Z_0 = U*NDIRS/V * dZ_J/dJ |_J=0.
 */
void measure_susceptibility(){
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
    // printf("measure_susceptibility: At site (%d,%d), field %d\n",t2,x2,field[t2][x2]);
     /*if(steps>100000){
       printf("Step number %d\n",steps);
       print_config();
       if(steps>100010) exit(1);
     }*/
     /* Now we are at (t2,x2), and the link is off. Try to move. */
     int dir = NDIRS*mersenne();
     int t3 = tdir(t2,dir), x3 = xdir(x2,dir);

     //printf("measure_susceptibility: Trying site (%d,%d), field %d\n",t3,x3,field[t3][x3]);

     if( mersenne() < 1 ){   //0.5 ){
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

     } else if( field[t3][x3] == 0 ) {
       /* Not using, makes no sense to accept moves every time */
       // First the current monomer becomes a Dirac link
       diraclink[t2][x2] = dir;
       field[t2][x2] = 0;


       //There should be a second link pointing at t3,x3
       int t4,x4; //New site
       for( int nu = 0; nu<NDIRS; nu++ ) {
         t4 = tdir(t3, nu), x4 = xdir(x3, nu);
         if( diraclink[t4][x4] == opp_dir(nu) ) if( !( t4==t2 && x4==x2) ){
           //Found a site pointing to t3,x3
           break;
         }
       }
       //Put the monomer there
       field[t4][x4] = SOURCE_MONOMER;
       diraclink[t4][x4] = NDIRS;

       //move to the new site
       t2 = t4; x2 = x4;

       //Now there will be a Dirac link pointing at the monomer.
       int t5,x5; //Problem site
       for( int nu = 0; nu<NDIRS; nu++ ) {
         t5 = tdir(t4, nu), x5 = xdir(x4, nu);
         if( diraclink[t5][x5] == opp_dir(nu) ){
           //Found a new site pointing to t4,x4
           break;
         }
       }
 
       //printf("Site %d %d pointing at the monomer\n",t5,x5);
       
       // Start a worm to find a new configuration 
       int done=1;
       do {
         //Pick a random direction to turn the link
         int mu = mersenne()*NDIRS;
         int t6 = tdir(t5, mu), x6 = xdir(x5, mu);
         if( field[t6][x6] == 0 ) { //Cannot point at a monomer
           //flip the link
           diraclink[t5][x5] = mu;
           //See if we now have two sites pointing to the same place
           done = 1;
           int t7,x7; //New site
           for( int nu = 0; nu<NDIRS; nu++ ) {
             t7 = tdir(t6, nu), x7 = xdir(x6, nu);
             if( diraclink[t7][x7] == opp_dir(nu) ){
               //Found a site pointing to t6,x6
               if( !( t7==t5 && x7==x5) ){
                 //And it's a new one, flip it next
                 done = 0;
                 break;
               }
             }
           }

           //move to the new site
           t5 = t7; x5 = x7;
         } else {
           done = 0;
         }
       } while( done == 0 );  //Quit if no new site is found
       
    } //neighbouring site (t3,x3)
   } //steps 
  // print_config();
  } //First site site (t2,x2)
 } //attempts


 printf("Susceptibility %g \n",(double)steps/(U*2*NDIRS*n_attempts));
  
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
      //if(empty==1) {
        if(diraclink[t][x]==TUP) printf(" v ");
        if(diraclink[t][x]==TDN) printf(" ^ ");
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

  /* "Warm up" the rng generator */
  for (i=0; i<543210; i++) mersenne();

  printf(" \n++++++++++++++++++++++++++++++++++++++++++\n");
  printf(" 4D free fermion, ( %d , %d ) lattice\n", NT, NX );
  printf(" %d updates per measurements\n", n_measure );
  printf(" m %f \n", m);
  printf(" U %f \n", U);
  printf(" mu %f \n", mu);
  printf(" Random seed %ld\n", seed );


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

  /* fill monomers and links */
  for (int t=0; t<NT; t++) for (int x=0; x<NX; x++) {
    field[t][x] = 0;
    diraclink[t][x] = 2*(x%2);
  }
  
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

      //print_config();
      
      /* Statistics */
      printf("MONOMERS %d \n", n_monomers);
      printf("LINKS %d \n", n_links);

      measure_charge();
      measure_susceptibility();

      gettimeofday(&start,NULL);
    }
  }

  printf(" ** simulation done\n");

  return(1);
}
















