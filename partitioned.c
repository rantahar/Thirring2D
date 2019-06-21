#define MAIN

#include "Thirring.h"


int **color;
int next_color=0;


void setup_lattice( long seed ){
    /* "Warm up" the rng generator */
    seed_mersenne( seed );
    for (int i=0; i<543210; i++) mersenne();

    field = malloc( NT*sizeof(int *) );
    color = malloc( NT*sizeof(int *) );
    eta = malloc( NT*sizeof(int *) );
    for (int t=0; t<NT; t++){
        field[t] = malloc( (NX+1)*sizeof(int) );
        color[t] = malloc( (NX+1)*sizeof(int) );
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

    /* fill up the index array */
    for(int i=0; i<NT; i++) {
      tup[i] = (i+1) % NT;
      tdn[i] = (i-1+NT) % NT;
    }
    for(int i=0; i<NX+1; i++) {
      xup[i] = (i+1) % NX;
      xdn[i] = (i-1+NX) % NX;
    }
    #ifdef OPENX  //Never match boundaries to actual sites
    xdn[0] = NX;
    xup[NX-1] = NX;
    xup[NX] = NX;
    #endif

    for(int t=0; t<NT; t++) for (int x=0; x<NX; x++) {
      field[t][x] = 0;
    }
    #ifdef OPENX
    for (int t=0; t<NT; t++) {
      field[t][NX] = EMPTY; //Site doesn't exist, no links or monomers, but not free either
    }
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
}


/* Determinant for given color */
double determinant( int target_color ){
  int *ipiv;
  int info;
  double *M;
  double det=1;

  struct timeval start, end;
  gettimeofday(&start,NULL);

  ipiv = malloc( VOLUME*sizeof(int) );
  M = malloc( VOLUME*VOLUME*sizeof(double) );


    int i=0,n=0;
    for (int t1=0; t1<NT; t1++) for (int x1=0; x1<NX; x1++)
    if(color[t1][x1]==target_color)
    {
        for (int t2=0; t2<NT; t2++) for (int x2=0; x2<NX; x2++)
        if(color[t2][x2]==target_color)
        {
            M[i] = fM_index( t1, x1, t2, x2, mu );
            //printf(" % .4f ",M[i]);
            i++;
        }
        //printf("\n");
        n++;
    }

    if(n>0){
        /* determinant from LU */
        LAPACK_dgetrf( &n, &n, M, &n, ipiv, &info );
        for( int a=0; a<n; a++) {
          det *= M[a*n+a];
          if(ipiv[a]!=a+1) det *= -1;
        }

        free(ipiv);
        free(M);
    } else {
        return 0;
    }

    gettimeofday(&end,NULL);
    int diff = 1e6*(end.tv_sec-start.tv_sec) + end.tv_usec-start.tv_usec;
    //printf("Calculated determinant in %.3g seconds\n", 1e-6*diff );

    return( det );
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
      if(field[t][x]==EMPTY) { empty = 0; printf("*"); }
      if(color[t][x] >= 0) printf(" %2d ", color[t][x]);
      if(color[t][x] == -1) printf("  x ");
    }
    printf(" \n");
  }
  printf(" \n");
}




int connected_with(int t, int x, int target_color){
    int connected = 0;
    for(int dir=0; dir<NDIRS; dir++){
        int t2 = tdir(t,dir), x2 = xdir(x,dir);
        if(color[t2][x2] == target_color){
            connected = 1;
            break;
        }
    }
    return connected;
}


double weight(int target_color);
int level = 0;
int verbose_level = 4;

double subconfigs(int target_color, int my_color,
    int sites[VOLUME][2], int nsites, int i0)
{
    double weight_sum = 0;

    for(int i=i0+1; i<nsites; i++) {
        int t = sites[i][0];
        int x = sites[i][1];
        if(connected_with(t, x, my_color)){
            color[t][x] = my_color;
            weight_sum += subconfigs(target_color, my_color, sites, nsites, i);

            double target_weight = weight(target_color);
            if(target_weight>0){
                double split_weight = weight(my_color);
                if(level < verbose_level){
                  print_config();
                  printf("%d weight %g %g\n", level, target_weight, split_weight);
                }
                target_weight *= target_weight;
            }

            weight_sum += target_weight;
            if(level < verbose_level)
                printf("%d subweights %g %g\n", level, weight_sum, target_weight);
            color[t][x] = target_color;
        }
    }

    return weight_sum;
}

double weight(int target_color){
    int sites[VOLUME][2];
    int nsites=0;
    double color_weight = 0;
    level +=1;

    for(int t=0; t<NT; t++) for(int x=0; x<NX; x++)
    if(color[t][x] == target_color)
    {
        sites[nsites][0] = t;
        sites[nsites][1] = x;
        nsites++;
    }

    if(level < verbose_level)
        printf(" %d sites\n", nsites);

    if(nsites > 0 && nsites%2==0){
      color_weight = determinant(target_color);
      if(level < verbose_level)
        printf("%d full weight %g\n", level, color_weight);

      /* See if we can split off two areas */
      int my_color = next_color;
      next_color += 1;
      int t = sites[0][0];
      int x = sites[0][1];
      color[t][x] = my_color;
      double split_weight = subconfigs(target_color, my_color, sites,nsites, 0);
      color_weight -= split_weight;
      if(level < verbose_level){
        printf("%d split weight %g\n", level, split_weight);
        printf("%d subtracted weight %g\n", level, color_weight);
      }

      for(int t=0; t<NT; t++) for(int x=0; x<NX; x++)
      if(color[t][x] == my_color)
      {
          color[t][x] = target_color;
      }
      next_color -= 1;
    }
    level -=1;
    return color_weight;
}


int main(int argc, char* argv[])
{
    U = 0.3;
    m = 0.1;
    mu = 0.0;
    long seed = 15423;
    setup_lattice(seed);

    for( int t=0; t<NT; t++) for( int x=0; x<NX; x++){
        int t2 = t/2;
        int x2 = x/3;
        color[t][x] = t2*NT/2+x2;
        next_color = color[t][x] +1;
    }

    print_config();

    double det = weight( 0 );
    printf("Determinant %f\n", det);

    return 0;
}
