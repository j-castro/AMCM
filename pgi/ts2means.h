/*
Parallel C implementation of the ts2means algorithm, originally
implemented in MATLAB by Professor Arturo Camacho.

JORGE CASTRO C. <jcastro@cenat.ac.cr>
CNCA - CENAT

*/
#ifndef ts2means_2017_CNCA
#define ts2means_2017_CNCA

#include <omp.h>        //add -fopenmp to compiler options
#include <stdio.h>
#include <math.h>
//#include <stdlib.h>
#include <time.h>
#include <sndfile.h>    //http://www.mega-nerd.com/libsndfile/ audio processing library
#include "vec.h"
#include <string.h>

#ifndef max
	#define max( a, b ) ( ((a) > (b)) ? (a) : (b) )
#endif

//minimum samples per window
#define _LOWERLIMIT_ 4

typedef struct {double inicio; double fin;} intervalo;


vec ts2meansP(vec,double,intervalo,double,double,char*);
vec ts2means(vec,double,intervalo,double,double,char*);
vec ts2means2(vec);
vec ts2means3(vec,double);
void ts2means4(double*,int,double,double*); //wraper for codes based on primitive double*
vec w2means(vec,vec);
vec getWSvec(intervalo,double,double,int);
vec classifySig(vec*, matrix*,matrix*,matrix*);
void postProcesing(vec*,vec*,char*);
int verificacion(double, double, double, char*,intervalo);


#endif
