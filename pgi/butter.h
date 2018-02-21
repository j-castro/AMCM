#ifndef butter_2017_manatee
#define butter_2017_manatee

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

double *binomial_mult( int n, double *p );
double *trinomial_mult( int n, double *b, double *c );

double *dcof_bwlp( int n, double fcf );
double *dcof_bwhp( int n, double fcf );

double *ccof_bwlp( int n );
double *ccof_bwhp( int n );

double sf_bwhp( int n, double fcf );

void butterHP(int, double, double*, double*);
void filterS(int, double*,double*,int,double*,double*);

#endif
