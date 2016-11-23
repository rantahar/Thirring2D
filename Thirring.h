#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mersenne.h"
#include <lapacke/lapacke.h>
#include <time.h>
#include <sys/time.h>

/* Lattice size, adjust */
#define NT 32
#define NX 32
#define ND 2
#define NDIRS (2*ND)

#define VOLUME (NT*NX)

#define MAX_CHANGES 40

