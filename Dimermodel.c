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
double U;

/* Monomers and links
 * field stores both, 0 for empty, 1 for monomer and 2+dir for links
 */
int **field;

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
    }
    printf(" \n");
  }
  printf(" \n");
  //usleep(100000);
}


void measure_charge(int *c, int *q){
  int _c = 0;
  int _q = 0;
  for (int x=0; x<NX; x++){
    if( field[0][x] == LINK_TUP ){
      _q+= 2*((x)%2==1?1:-1);
    }
  }
  *c = _c;
  *q = _q;
}


n_links = VOLUME/2;

double measure_susceptibility_with_background( ){
  int nsteps = 0;
  int n_attempts = 1;
  double scale_factor = n_links/(2.*VOLUME*n_attempts);
  
  for(int i=0; i<n_attempts;i++){
    //Pick a site
    int t=0,x=0;
    if(n_links>0) do{
      t= mersenne()*NT, x=mersenne()*NX;
    } while(field[t][x] < LINK_TUP);
  
    //printf("STARTING susceptibility at (%d,%d)\n",t,x);
    if( field[t][x] >= LINK_TUP){
      //Replace the link with two source monomers
      int dir = field[t][x]-2;
      int t2 = tdir(t, dir), x2 = xdir(x, dir);

      field[t][x] = SOURCE_MONOMER;
      field[t2][x2] = SOURCE_MONOMER;
    
      for(;;){
        //print_config();
        nsteps++;
        int dir = mersenne()*NDIRS;
        int t3 = tdir(t2, dir), x3 = xdir(x2, dir);
      
        if( t3==t && x3==x ) {
           // Back at the original site, turn into a link 
           field[t][x]=0; field[t2][x2]=0;
           link_on(t2,x2,dir);
           break;

         } else if( field[t3][x3] >= LINK_TUP ) {
           // found a link, flip it  
           int linkdir = field[t3][x3] - LINK_TUP;
           int t4 = tdir(t3,linkdir), x4 = xdir(x3,linkdir);

           field[t2][x2] = 0;
           link_off(t3,x3,linkdir); link_on(t2,x2,dir);
           field[t4][x4] = SOURCE_MONOMER;

           t2 = t4; x2 = x4;

         } 
         //for(int n=0; n<10; n++) update_dirac_background();
      }  
    }
  }

  // print_config();
  return (double)nsteps*scale_factor;
}









/* Main function
 */
int main(int argc, char* argv[])
{
  #ifdef DEBUG
  feenableexcept(FE_INVALID | FE_OVERFLOW);
  #endif 

  int i,n_loops,n_average;
  long seed;

  /* Read in the input */
  printf(" Number of updates : ");
  scanf("%d",&n_loops);

  printf(" Average over : ");
  scanf("%d",&n_average);
  printf("%d\n",n_average);

  printf(" Random number : ");
  scanf("%ld",&seed);
  seed_mersenne( seed );


  /* "Warm up" the rng generator */
  for (i=0; i<543210; i++) mersenne();

  printf(" \n++++++++++++++++++++++++++++++++++++++++++\n");
  printf(" 2D Dimer model, ( %d , %d ) lattice\n", NT, NX );
  printf(" Random seed %ld\n", seed );

  /* Allocate location and field arrays */
  field = malloc( NT*sizeof(int *) );
  for (int t=0; t<NT; t++){
    field[t] = malloc( (NX+1)*sizeof(int) );
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

  /* fill monomers and links */
  for (int t=0; t<NT; t++) if(t%2==0) for (int x=0; x<NX; x++) {
    link_on(t,x,TUP);
  }
  
#ifdef OPENX
  for (int t=0; t<NT; t++) {
    field[t][NX] = EMPTY; //Site doesn't exist, no links or monomers, but not free either
  }
#endif
  
  /* and the update/measure loop */
  int sum_charge = 0;
  int sum_c2 = 0;
  int sum_q = 0;
  int sum_q2 = 0;
  double sum_susc_wb = 0;

  struct timeval start, end;
  double updatetime=0, measuretime = 0;
  for (i=1; i<n_loops+1; i++) {
    //print_config();

      /* Time */
      gettimeofday(&start,NULL);

      int c, q;
      measure_charge(&c, &q);
      sum_charge += c;
      sum_c2 += c*c;
      sum_q += q;
      sum_q2 += q*q; 

      sum_susc_wb += measure_susceptibility_with_background();
      gettimeofday(&end,NULL);
      measuretime += 1e6*(end.tv_sec-start.tv_sec) + end.tv_usec-start.tv_usec;

      if((i%(n_average))==0){
        printf("%d, %d measurements in %.3g seconds\n", i, n_average, 1e-6*measuretime);
        updatetime = 0; measuretime = 0;

        printf("CHARGE %g %g \n", (double)sum_charge/n_average, (double)sum_c2/n_average);
        printf("QCHI %g %g \n", (double)sum_q/n_average, (double)sum_q2/n_average);
        printf("SUSCEPTIBILITY %g \n", (double)sum_susc_wb/n_average);

        sum_charge = 0; sum_c2 = 0; sum_q = 0; sum_q2 = 0; sum_susc_wb = 0;
      }
      
      gettimeofday(&start,NULL);
    
  }

  printf(" ** simulation done\n");

  return(1);
}















