/* Copyright (c) 2009-2011 Kyle Gorman
*
*  Permission is hereby granted, free of charge, to any person obtaining a copy
*  of this software and associated documentation files (the "Software"), to deal
*  in the Software without restriction, including without limitation the rights
*  to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
*  copies of the Software, and to permit persons to whom the Software is
*  furnished to do so, subject to the following conditions:
*
*  The above copyright notice and this permission notice shall be included in
*  all copies or substantial portions of the Software.
*
*  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
*  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
*  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
*  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
*  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
*  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
*  THE SOFTWARE.
*
*  vec: some data structures for swipe
*  Kyle Gorman <kgorman@ling.upenn.edu>
*
*  This is version 1.0., i.e. I think I got all obvious stuff working ideally.
*/
#ifndef changoloco2017 
#define changoloco2017
//getSubvACC, autocorrvACC and powvACC work with double pointers (for GPU compatibility)

// vec stuff 
typedef struct                   { int x; double * v; } vec;

vec                            	makev(int);
vec                            	zerov(int);
vec                            	onesv(int);
vec                            	nansv(int);
vec                            	makeFillv(double,double,int);
vec				multv(vec,double);
//#pragma acc routine seq 
vec				roundv(vec);
vec                            	makeHannv(int);
vec                            	copyv(vec);
vec                            	copyv2(double *, int);
vec                            	addElementv(vec,double);
vec				dotpv(vec,vec);
//#pragma acc routine seq
vec                            	getSubv(vec,int,int);
vec                            	distancev(vec,vec);
//#pragma acc routine seq
vec 				powv(vec,int);
//#pragma acc routine seq
vec 				autocorrv(vec,int,int);
/*#pragma acc routine seq
void                            getSubvACC(double*,int,int,int,double*,int);
#pragma acc routine seq
void 				powvACC(double*,int,int,double*,int);
#pragma acc routine seq
void 				autocorrvACC(double*,int,int,int,double*,int);*/
//#pragma acc routine seq 
vec				thresholdv(vec,double);
//#pragma acc routine seq
vec				removeSSv(vec, int);
//#pragma acc routine seq 
vec				removeLSv(vec, int);
//#pragma acc routine seq
vec				softThresv(vec,double);

void 				multv2(vec&, double);
void   				resizev(vec *, int);
//#pragma acc routine seq
void                            absv(vec*);

double                          maxvv(vec);
/*#pragma acc routine seq
double                          maxvvACC(double*,int);*/
double                          minvv(vec);
double 				meanv(vec);
/*#pragma acc routine seq
double 				meanvACC(double*,int);
double                          meanvACCint(int*,int); */
double				madv(vec);
int                             maxv(vec);
int                             minv(vec);
int                             bisectv(vec, double);
int                             bilookv(vec, double, int);
int                             equalv(vec,vec);

void                              freev(vec);
void                              printv(vec);

// intvec stuff
typedef struct                    { int x; int* v; } intvec;

intvec                         makeiv(int);
intvec                         zeroiv(int);
intvec                         onesiv(int);
intvec                         copyiv(intvec);

vec                            iv2v(intvec);

int                               maxiv(intvec);
int                               miniv(intvec);
int                               bisectiv(intvec, int);
int                               bilookiv(intvec, int, int);

void                              freeiv(intvec);
void                              printiv(intvec);

// matrix stuff
typedef struct                   { int x; int y; double** m; } matrix;

matrix                            makem(int, int);
matrix                            zerom(int, int);
matrix                            onesm(int, int);
matrix                            nansm(int, int);
//vec* 						  mat2vect(matrix);
matrix                            copym(matrix);
void                              freem(matrix);
void                              printm(matrix);
double                            maxm(matrix);
void                              copyRow(matrix*, vec,int);
void                              copyRow2(matrix&, vec,int);
vec				  getRowm(matrix, int);
matrix 				  add(matrix, matrix);
void				  movAvrM(matrix&, int); //TODO check if it is better to use a pointer

// intmatrix stuff
typedef struct                    { int x; int y; int** m; } intmatrix;

intmatrix                         makeim(int, int);
intmatrix                         zeroim(int, int);
intmatrix                         onesim(int, int);
intmatrix                         copyim(intmatrix);

matrix                            im2m(intmatrix); // cast

void                              freeim(intmatrix);
void                              printim(intmatrix);

// prime sieve
#define P                          1
#define NP                         0

#define PRIME(x)                   (x == 1)

int                                sieve(intvec);
intvec                          primes(int);

// cubic spline
#define YP1                        2.
#define YPN                        2.

vec                             spline(vec, vec);
double                             splinv(vec, vec, vec, double, int);
vec                              interp1(vec, vec, vec);

// polynomial fitting
vec                             polyfit(vec, vec, int);
double                             polyval(vec, double);
#endif

