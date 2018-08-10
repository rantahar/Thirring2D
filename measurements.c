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
extern int    ***eta;   //Staggered eta matrix
double *Dinv;             //Storage for even to odd propagator
int    *evenlist, *oddlist;  //Lists of occupied sites
int    *unoccupied_evenlist, *unoccupied_oddlist;  //Lists of occupied sites
int *added_evensites, *added_oddsites, *removed_evenlist, *removed_oddlist;
int n_added_even,n_added_odd,n_removed_even,n_removed_odd;
extern double m;
extern double U;
extern double mu;

/* Maximum number of fluctuations from the background configuration */
int max_changes;

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
double *Ginv;


/* Neighbour index arrays, to be filled at the beginning
 */
extern int *tup,*xup,*tdn,*xdn;

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


/* Functions for fetching neighboring coordinates */
static inline int tdir(int t, int dir){
  if( dir == TUP ) return tup[t];
  if( dir == TDN ) return tdn[t];
  return(t);
}

static inline int xdir(int x, int dir){
  if( dir == XUP ) return xup[x];
  if( dir == XDN ) return xdn[x];
  return(x);
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








/* Functions for measurements
 */

/* Propagator
 */
void measure_propagator(){ 
  double prop[NT];
  double j[NT];
  double c[NT];
  //double q[NT];
  double boson[NT];
  for( int t1=0; t1<NT; t1++){
    prop[t1]=0; j[t1] = 0; c[t1]=0; //q[t1] = 0;
    boson[t1] = 0;
  }

  /*for( int t1=0;t1<NT;t1++) for( int x1=0;x1<NX;x1++) {
   field[t1][x1]=0;
  }
  field[3][0]=LINK_XUP;
  field[3][1]=LINK_XDN;
  field[0][3]=LINK_TDN;
  field[3][3]=LINK_TUP;
  */

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

      /*if( t1==0 && x1==0 ){
        printf("Dinv %g\n",propagator[0][1]);
        printf("Dinv %g\n",propagator[1][0]);
        printf("Dinv %g\n",propagator[3][0]);
      }*/

      for( int t2=0; t2<NT; t2++) prop[(t2-t1+NT)%NT] += propagator[t2][x1];
      for( int t2=0; t2<NT; t2++) boson[(t2-t1+NT)%NT] += propagator[t2][x1]*propagator[t2][x1];

      j[t1] += bc_up*exp(mu)*eta[t1][x1][0]*propagator[tup[t1]][x1];
      c[t1] += bc_up*((x1+t1)%2 ==0 ? 1:-1 )* exp(mu)*eta[t1][x1][0]*propagator[tup[t1]][x1];
      //q[t1] += bc_up*((x1+t1)%2 ==0 ? 1:-1 )* exp(mu)*eta[t1][x1][0]*propagator[tup[t1]][x1];

      j[tdn[t1]] += bc_dn*exp(-mu)* eta[tdn[t1]][x1][0]*propagator[tdn[t1]][x1];
      c[tdn[t1]] -= bc_dn*((x1+tdn[t1])%2 ==0 ? 1:-1 )* exp(-mu)* eta[tdn[t1]][x1][0]*propagator[tdn[t1]][x1];
      
    }

    if( field[t1][x1]== LINK_TUP ) { 
      double cl = ((t1+x1)%2==0 ? 1:-1 )*2;
      c[t1] += 2*cl;
      //q[t1] += cl;
    }
  }

  
  //for( int t2=0; t2<NT; t2++) printf("Propagator %d %g\n", t2, current_sign*prop[t2]/(VOLUME) );
  //for( int t2=0; t2<NT; t2++) printf("Boson %d %g\n", t2, current_sign*boson[t2]/(VOLUME) );
  for( int t1=0;t1<1;t1++) printf("Charge %d %g\n", t1, current_sign*j[t1]/2 );
  for( int t1=0;t1<1;t1++)  printf("Qchi %d %g\n", t1, current_sign*c[t1]/2 );
  //for( int t1=0;t1<1;t1++)  printf("qchi %d %g\n", t1, current_sign*q[t1] );
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
void measure_susceptibility( double * susc){
    
 static int init = 1;
 if( init == 1 ){
   #ifdef FLUCTUATION_MATRIX
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
   #endif
   #ifdef PROPAGATOR_MATRIX
   int ** field_copy = malloc( NT*sizeof(int *) );
   for (int t=0; t<NT; t++) field_copy[t] = malloc( (NX+1)*sizeof(int) );
   for (int t=0; t<NT; t++) for (int x=0; x<NX; x++) {
     field_copy[t][x] = field[t][x];
     field[t][x] = 0;
   }
   #ifdef OPENX
   for (int t=0; t<NT; t++) {
     field_copy[t][NX] = EMPTY; //Site doesn't exist, no links or monomers, but not free either
   }
   #endif
   calc_Dinv( );
   for (int t=0; t<NT; t++) for (int x=0; x<NX; x++) {
     field[t][x] = field_copy[t][x];
   }
   for (int t=0; t<NT; t++) free(field_copy[t]);
   free(field_copy);
   #endif
   
   init = 0;
 }
 update_linklists();
 update_background();
 
    
 int n = n_bc_monomers/2 + n_bc_links;
 int steps = 0;
 int n_attempts=10;
 double scale_factor = (double)n_links/(NDIRS*n_attempts*VOLUME);

 /* Do multiple attemps, these are cheap and the result is usually 0 */
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

 *susc = (double)steps*scale_factor;
  
}




























