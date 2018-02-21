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
*  vec: some data structures and mathematical function
*  Kyle Gorman <kgorman@ling.upenn.edu>
*
*  If for some reason you didn't get vec.h, which you'll need to #include
*  vec, it is available at:
*
*               http://ling.upenn.edu/~kgorman/c/vec.h
*
*  This is version 1.0., i.e. I think I got all obvious stuff working ideally.
*/

#include <math.h>
#include <stdio.h>
#include <limits.h>
#include <string.h>
#include <stdlib.h>

#include "vec.h"

#ifndef NAN
    #define NAN     sqrt(-1.)
#endif

#ifndef pi
    #define pi 3.1415926535897932384626433832795
#endif

// create a vec of size xSz
vec makev(int xSz) {
    vec nw_vec;
    nw_vec.x = xSz;
    nw_vec.v = new double[xSz];
    return(nw_vec);
}



// make a vec of zeros of size xSz
vec zerov(int xSz) {

    vec nw_vec;
    nw_vec.x = xSz;
    nw_vec.v = new double[xSz]();

    return(nw_vec);
}

// make a vec of ones of size xSz
vec onesv(int xSz) {
    int i;
    vec nw_vec = makev(xSz);
    for (i = 0; i < nw_vec.x; i++) {
        nw_vec.v[i] = 1.;
    }
    return(nw_vec);
}

//make a vec starting in 'start', in steps of 'increment' with 'nIter' iterations
vec makeFillv(double start,double increment, int nIter){

    vec nw_vec;
    int i;

    if (nIter > 0 && increment > 0){   //parameters check

        nw_vec.x = nIter;
        nw_vec.v = new double[nIter];

        for (i = 0; i < nIter; i++) {
            nw_vec.v[i] = start + (double)i*increment;
        }
    }

    else {
        fprintf(stderr, "Wrong parameters values: 'nIter' and 'increment' must be greater than 0\n");
    }

    return(nw_vec);
}

//Multiplication of a vector by a constant
vec multv(vec in_vec, double alpha){
    vec nw_vec = makev(in_vec.x);
    
    for(int i=0; i < in_vec.x; i++) {
        nw_vec.v[i] = in_vec.v[i] * alpha;
    }

    return nw_vec;
}

//Multiplication of a vector by a constant
void multv2(vec & in_vec, double alpha){
    
    for(int i=0; i < in_vec.x; i++) {
        in_vec.v[i] = in_vec.v[i] * alpha;
    }
}

//Round each element of a vector to the nearest integer
vec roundv(vec in_vec){
    vec nw_vec = makev(in_vec.x);
    
    for(int i=0; i < in_vec.x; i++) {
        nw_vec.v[i] = llround(in_vec.v[i]);
    }

    return nw_vec;
}


//make a Hann window
vec makeHannv(int xSz){

    vec nw_vec;
    int i;

    nw_vec.x = xSz;
    nw_vec.v = new double[xSz];

	for (i = 1; i <= xSz; i++) {
		nw_vec.v[i-1] = .5 * (1 - cos((2 * pi * i) / (xSz + 1)));
	}

    return(nw_vec);
}

// make a vec of NaNs of size xSz
vec nansv(int xSz) {
    int i;
    vec nw_vec = makev(xSz);
    for (i = 0; i < nw_vec.x; i++) {
        nw_vec.v[i] = NAN;
    }
    return(nw_vec);
}

//Copy a vector
vec copyv(vec yr_vec) {
    vec nw_vec = makev(yr_vec.x);
    //memcpy(nw_vec.v, yr_vec.v, sizeof(double) * yr_vec.x);
    for (int i =0; i < yr_vec.x; i++) {
	nw_vec.v[i] = yr_vec.v[i];
    }
    return(nw_vec);
}

//Copy a double * into a vec type
vec copyv2(double * in, int sz) {
   vec nw_vec = makev(sz);

    for (int i =0; i < sz; i++) {
	    nw_vec.v[i] = in[i];
    }
    return(nw_vec);
}


// free the memory associated with the vec
void freev(vec yr_vec) {
    delete[] yr_vec.v;
}

// print the vec
void printv(vec yr_vec) {
    int i;
    for (i = 0; i < yr_vec.x; i++) {
        printf("%f\n", yr_vec.v[i]);
    }
}

//return a new vec with this new value added
vec addElementv(vec v, double value){

    vec nw_vec = makev(v.x + 1);
    int i;
    for (i = 0;i < v.x ;i++){

        nw_vec.v[i] = v.v[i];

    }

    nw_vec.v[nw_vec.x -1] = value;
    freev(v);

    return nw_vec;
}

//Apply soft thresholding rule to vector "in" using threshold "t"
vec softThresv(vec in,double t) {
    
    vec out = makev(in.x);

    for (int i=0; i < in.x; i++) {
	if (fabs(in.v[i]) <= t)
	    out.v[i] = 0;
	else if (in.v[i] > t)
	    out.v[i] = in.v[i] - t;
	else 
	    out.v[i] = in.v[i] + t;
    }

    return out;
}

// return the index of the maximum value of the vec
int maxv(vec yr_vec) {
    int i;
    int index;
    double val = SHRT_MIN;
    for (i = 0; i < yr_vec.x; i++) {
        if (yr_vec.v[i] > val) {
            val = yr_vec.v[i];
            index = i;
        }
    }
    return(index);
}

// return the value of the maximum value of the vec
double maxvv(vec yr_vec) {
    int i;
    int index;
    double val = SHRT_MIN;
    for (i = 0; i < yr_vec.x; i++) {
        if (yr_vec.v[i] > val) {
            val = yr_vec.v[i];
            index = i;
        }
    }
    return yr_vec.v[index];
}


// return the value of the maximum value of the double * in
/*double maxvvACC(double* in, int sz) {
    int i;
    int index;
    double val = SHRT_MIN;
    for (i = 0; i < sz; i++) {
        if (in[i] > val) {
            val = in[i];
            index = i;
        }
    }
    return in[index];
}*/
// return the index of the minimum value of the vec
int minv(vec yr_vec) {
    int i;
    int index;
    double val = SHRT_MAX;
    for (i = 0; i < yr_vec.x; i++) {
        if (yr_vec.v[i] < val) {
            val = yr_vec.v[i];
            index = i;
        }
    }
    return(index);
}

// return the value of the minimum value of the vec
double minvv(vec yr_vec) {
    int i;
    int index;
    double val = SHRT_MAX;
    for (i = 0; i < yr_vec.x; i++) {
        if (yr_vec.v[i] < val) {
            val = yr_vec.v[i];
            index = i;
        }
    }
    return yr_vec.v[index];;
}

// return the median absolute deviation for vector "vec"
double madv(vec in) {

    double temp, med, mad;
    vec temp_vec = copyv(in);

    //Compute median of "in"
    for (int i = 0; i < in.x; i++) {
	
	for(int j = 1; j <   in.x-i; j++) {
	    if (temp_vec.v[j-1] > temp_vec.v[j]){ //exchange positions
		temp = temp_vec.v[j];
		temp_vec.v[j] = temp_vec.v[j-1];
		temp_vec.v[j-1] = temp; 
	    } 
	}
    }

    med = temp_vec.v[temp_vec.x/2];
    /*printf("Ordered vector:\n");
    printv(temp_vec);
    printf("\nMedian: %f\n",med);*/

    //Compute absolute deviation of the samples
    for (int i=0; i < in.x; i++) {
	temp_vec.v[i] = fabs(in.v[i] - med);
    }
    
    /*printf("Absolute deviation vector:\n");
    printv(temp_vec);*/

    //Compute median of the absolute deviation
    for (int i = 0; i < in.x; i++) {
	
	for(int j = 1; j <  in.x-i; j++) {
	    if (temp_vec.v[j-1] > temp_vec.v[j]){ //exchange positions
		temp = temp_vec.v[j];
		temp_vec.v[j] = temp_vec.v[j-1];
		temp_vec.v[j-1] = temp; 
	    } 
	}
    }

    mad = temp_vec.v[in.x/2];
    /*printf("\nMAD IN: %f\n",mad);*/
    freev(temp_vec);
    return mad;
}

void resizev(vec *x, int size){

	vec y = copyv(*x);
	freev(*x);
	*x = zerov(size);
	int i;
	for(i = 0; i < size; ++i){
		x->v[i] = y.v[i];
	}
	freev(y);
	//x->v = (double*)realloc(x->v,sizeof(double)*size);
}

//get a subvec from index: i1 to i2
vec getSubv(vec x, int i1, int i2){

    vec nw_vec;
    int i;

    if ((i2 >= i1) && (i2 < x.x) && (i1 >= 0)) { //Parameters check

        nw_vec = zerov((i2-i1) + 1);

        for(i = i1; i <= i2; ++i){
            nw_vec.v[i-i1] = x.v[i];
        }
    }
    else {
        //printf("values x.x: %d ; i1: %d ; i2: %d\n",x.x,i1,i2);
        //fprintf(stderr, "Wrong parameters values, follow the restriccions: \n1) 'i2' >= 'i1' \n2) 'i2' < size of x \n3) 'i1' >= 0\n");
    }

	return(nw_vec);
}


/*void powvACC(double * in, int szIn, int n, double * out, int szOut) {

    if (szIn <= szOut) {
        for(int i = 0;i < szIn; i++) {
	    out[i] = pow(in[i], n);
        }
    }
}*/

//Compute the sample autocorrelation of vector "in" from 
//first lag "fl" to last lag "ll"
/*void autocorrvACC(double* in, int szIn, int fl, int ll, double * out, int szOut){
  
    if ((ll-fl+1) <= szOut) {

        double mean = meanvACC(in,szIn);  	//Sample mean 
        double var = 0;		//Sample variance

        //Compute sample variance
        for (int i = 0; i < szIn; i++) {
	    var += (in[i]-mean)*(in[i]-mean);
        }
        var = var / szIn;

        //Compute sample-ACF
        for(int i = fl;i <= ll; i++) {// Lag loop	
	    out[i-fl] = 0;	
	    for(int j = i;j < szIn; j++) { //Dot product
	        out[i-fl] += (in[j] - mean) * (in[j-i] - mean); 
            }
	    out[i-fl] = out[i-fl] / (var * szIn);
        }
    }
}*/

//Get a sub-vector of array "x" from index "i1" to "i2"
//and store it in array "y"
/*void getSubvACC(double* x, int szX, int i1, int i2, double *y, int szY){

    if ((i2 >= i1) && (i2 < szX) && (i1 >= 0) && ((i2-i1+1) <= szY)) { //Parameters check

        for(int i = i1; i <= i2; i++){
            y[i-i1] = x[i];
        }
    }
}*/

vec				dotpv(vec,vec);

//returns the dot product between vec1 and vec2
vec dotpv(vec v1, vec v2){

    vec out = makev(v1.x);
    if (v1.x != v2.x) {
    	fprintf(stderr, "ERROR: vector sizes must be the same to apply a dot product");
    }
    else {
	for (int i=0; i<v1.x; i++) {
	    out.v[i] = v1.v[i] * v2.v[i];
	}
    }

    return out;
}

//returns distance vector between v1 and v2 
vec distancev(vec v1,vec v2){

    vec distance = zerov(v1.x);
    int i;

    if (v1.x == v2.x) {     //size check

        for(i = 0;i < v1.x; i++){
            distance.v[i] = fabs(v1.v[i] - v2.v[i]);
        }
    }
     else {

       fprintf(stderr, "Error: both vecs (v1,v2) must have the same length \n");

     }

    return distance;
}

//Applies a threshold "th" to all elements in vector "in" 
vec thresholdv(vec in, double th){

    vec out = zerov(in.x);
    for (int i = 0; i < in.x; i++) {
	if (in.v[i] > th)
	    out.v[i] = 1; 
    }

    return out;
}

//Remove sequences of 0's and 1's shorter than "minL" in
//a binary vector "bv"
vec removeSSv(vec bv, int minL){
    vec out = copyv(bv);
    int seqS = 0; //Sequence start index

    //Remove small regions of zeros
    for(int i=1;i < out.x; i++) {

	if ((out.v[i] - out.v[i-1]) != 0) {//is a region border?
	    
	    if ((out.v[i] - out.v[i-1]) == -1) {//is the start of a region of zeros?
		seqS = i;  //update the sequence start
	    }
	    else {//is the end of a region of zeros?
		if ((i - seqS) < minL) { //check minimum length restricction
		    for (int j = seqS; j < i;j++) //Replace zeros by ones
			out.v[j] = 1;	    
		}
	    }
	}	
    }

    seqS = 0;
    //remove small regions of ones
    for(int i=1;i < out.x; i++) {

	if ((out.v[i] - out.v[i-1]) != 0) {//is a region border?
	    
	    if ((out.v[i] - out.v[i-1]) == 1) {//is the start of a region of ones?
		seqS = i;  //update the sequence start
	    }
	    else {//is the end of a region of ones?
		if ((i - seqS) < minL) { //check minimum length restricction
		    for (int j = seqS; j < i;j++) //Replace ones by zeros
			out.v[j] = 0;   
		}
	    }
	}	
    }

    return out;
}

//Remove sequences of 1's longer than "maxL" in
//a binary vector "bv"
vec removeLSv(vec bv, int maxL){
    vec out = copyv(bv);
    int seqS = 0; //Sequence start index
    
    //remove long regions of ones
    for(int i=1;i < out.x; i++) {

	if ((out.v[i] - out.v[i-1]) != 0) {//is a region border?
	    
	    if ((out.v[i] - out.v[i-1]) == 1) {//is the start of a region of ones?
		seqS = i;  //update the sequence start
	    }
	    else {//is the end of a region of ones?
		if ((i - seqS) > maxL) { //check minimum length restricction
		    for (int j = seqS; j < i; j++) //Replace ones by zeros
			out.v[j] = 0;   
		}
	    }
	}	
    }

    return out;
}


//return absolute valor of the vec
void absv(vec * x){

    int i;
    for(i = 0; i < x->x; i++){
        x->v[i] = fabs(x->v[i]);
    }
}

//Return 1 if the vecs are equal
int equalv(vec v1, vec v2){

    int eq = 1;
    int i = 0;

    if (v1.x == v2.x) {//check sizes first
        while (eq == 1 && i < v1.x) {
            if (v1.v[i] != v2.v[i]) {
                eq = 0;
            }
            i++;
        }
    }
    else {
        eq = 0;
    }

    return eq;

}


// Return the mean of the vector "in"
double meanv(vec in) {
    
    double sum = 0;
    for(int i = 0;i < in.x; i++) {
	sum = sum + in.v[i];
    }

    return sum / in.x;
}

// Return the mean of double array "in"
/*double meanvACC(double* in, int sz) {
    
    double sum = 0;
    for(int i = 0;i < sz; i++) {
	sum = sum + in[i];
    }

    return sum / sz;
}*/

// Return the mean of int array "in"
/*double meanvACCint(int* in, int sz) {
    
    double sum = 0;
    for(int i = 0;i < sz; i++) {
	sum = sum + in[i];
    }

    return sum / sz;
}*/

//Return the vector "in" powered to "n"
vec powv(vec in, int n) {

    vec out = zerov(in.x);
    for(int i = 0;i < in.x; i++) {
	out.v[i] = pow(in.v[i], n);
    }

    return out;

}

//Compute the sample autocorrelation of vector "in" from 
//first lag "fl" to last lag "ll"
vec autocorrv(vec in, int fl, int ll){
  
    vec out = zerov(ll-fl+1); 	//Check error. in fl and ll
    double mean = meanv(in);  	//Sample mean 
    double var = 0;		//Sample variance

    //printf("Mean:  %f\n", mean);

    //Compute sample variance
    for (int i = 0; i < in.x; i++) {
	var += (in.v[i]-mean)*(in.v[i]-mean);
    }
    var = var / in.x;
    
    //printf("Variance:  %f\n", var);

    //Compute sample-ACF
    for(int i = fl;i <= ll; i++) {// Lag loop	
	for(int j = i;j < in.x; j++) { //Dot product
	    out.v[i-fl] += (in.v[j] - mean) * (in.v[j-i] - mean); 
        }
	out.v[i-fl] = out.v[i-fl] / (var * in.x);
    }

    return out;
}

// find the bisection index of the vec for key
int bisectv(vec yr_vec, double key) {
    int md;
    int lo = 1;
    int hi = yr_vec.x;
    while (hi - lo > 1) {
        md = (hi + lo) >> 1;
        if (yr_vec.v[md] > key) {
            hi = md;
        }
        else {
            lo = md;
        }
    }
    return(hi);
}

// like bisectv(), but the minimum starting value is passed as an argument. This
// is good for multiple bisection calls for forming a new vec when the
// queries are a non-constant interval; but make sure to use bisectv() the
// first time.
int bilookv(vec yr_vec, double key, int lo) {
    int md;
    int hi = yr_vec.x;
    lo--;
    while (hi - lo > 1) {
        md = (hi + lo) >> 1;
        if (yr_vec.v[md] > key) {
            hi = md;
        }
        else {
            lo = md;
        }
    }
    return(hi);
}

// intvec versions of the above

intvec makeiv(int xSz) {
    intvec nw_vec;
    nw_vec.x = xSz;
    nw_vec.v = new int[xSz];
    return(nw_vec);
}

intvec zeroiv(int xSz) {
    int i;
    intvec nw_vec = makeiv(xSz);
    for (i = 0; i < nw_vec.x; i++) {
        nw_vec.v[i] = 0;
    }
    return(nw_vec);
}

intvec onesiv(int xSz) {
    int i;
    intvec nw_vec = makeiv(xSz);
    for (i = 0; i < nw_vec.x; i++) {
        nw_vec.v[i] = 1;
    }
    return(nw_vec);
}

intvec copyiv(intvec yr_vec) {
    intvec nw_vec = makeiv(yr_vec.x);
    memcpy(nw_vec.v, yr_vec.v, sizeof(int) * nw_vec.x);
    return(nw_vec);
}

// convert an intvec into a vec using implicit casts to double
vec iv2v(intvec yr_vec) {
    int i;
    vec nw_vec = makev(yr_vec.x);
    for (i = 0; i < yr_vec.x; i++) {
        nw_vec.v[i] = yr_vec.v[i];
    }
    return(nw_vec);
}

void freeiv(intvec yr_vec) {
    delete[] yr_vec.v;
}

void printiv(intvec yr_vec) {
    int i;
    for (i = 0; i < yr_vec.x; i++) {
        printf("%d\n", yr_vec.v[i]);
    }
}

int maxiv(intvec yr_vec) {
    int i;
    int index;
    int val = SHRT_MIN;
    for (i = 0; i < yr_vec.x; i++) {
        if (yr_vec.v[i] > val) {
            val = yr_vec.v[i];
            index = i;
        }
    }
    return(index);
}

int miniv(intvec yr_vec) {
    int i;
    int index;
    int val = SHRT_MAX;
    for (i = 0; i < yr_vec.x; i++) {
        if (yr_vec.v[i] < val) {
            val = yr_vec.v[i];
            index = i;
        }
    }
    return(index);
}

int bisectiv(intvec yr_vec, int key) {
    int md;
    int lo = 1;
    int hi = yr_vec.x;
    while (hi - lo > 1) {
        md = (hi + lo) >> 1;
        if (yr_vec.v[md] > key) {
            hi = md;
        }
        else {
            lo = md;
        }
    }
    return(hi);
}

int bilookiv(intvec yr_vec, int key, int lo) {
    int md;
    int hi = yr_vec.x;
    lo--;
    while (hi - lo > 1) {
        md = (hi + lo) >> 1;
        if (yr_vec.v[md] > key) {
            hi = md;
        }
        else {
            lo = md;
        }
    }
    return(hi);
}

// matrix versions of the above

//TODO check this memory
matrix makem(int xSz, int ySz) {
    int i;
    matrix nw_matrix;
    nw_matrix.x = xSz;
    nw_matrix.y = ySz;
    nw_matrix.m = new double*[xSz];
    for (i = 0; i < nw_matrix.x; i++) {
        nw_matrix.m[i] =  new double[ySz];
    }
    return(nw_matrix);
}

matrix zerom(int xSz, int ySz) {
    int i;
    int j;
    matrix nw_matrix = makem(xSz, ySz);
    for (i = 0; i < nw_matrix.x; i++) {
        for (j = 0; j < nw_matrix.y; j++) {
            nw_matrix.m[i][j] = 0.;
        }
    }
    return(nw_matrix);
}

matrix onesm(int xSz, int ySz) {
    int i;
    int j;
    matrix nw_matrix = makem(xSz, ySz);
    for (i = 0; i < nw_matrix.x; i++) {
        for (j = 0; j < nw_matrix.y; j++) {
            nw_matrix.m[i][j] = 1.;
        }
    }
    return(nw_matrix);
}

matrix nansm(int xSz, int ySz) {
    int i;
    int j;
    matrix nw_matrix = makem(xSz, ySz);
    for (i = 0; i < nw_matrix.x; i++) {
        for (j = 0; j < nw_matrix.y; j++) {
            nw_matrix.m[i][j] = NAN;
        }
    }
    return(nw_matrix);
}

//TODO check this method for posibble error
matrix copym(matrix yr_matrix) {
    int i;
    matrix nw_matrix = makem(yr_matrix.x, yr_matrix.y);
    for (i = 0; i < yr_matrix.x; i++) { // does not assume contiguous memory
        memcpy(nw_matrix.m[i], yr_matrix.m[i], sizeof(double) * yr_matrix.y);
    }
    return(nw_matrix);
}

//TODO check this method for posibble error
void freem(matrix yr_matrix) {
    int i;
    for (i = 0; i < yr_matrix.x; i++) {
        delete[] yr_matrix.m[i];
    }
    delete yr_matrix.m;
}

void printm(matrix yr_matrix) {
    int i;
    int j;
    for (i = 0; i < yr_matrix.x; i++) {
        for (j = 0; j < yr_matrix.y; j++) {
            printf("%f\t", yr_matrix.m[i][j]);
        }
        printf("\n");
    }
}

// intmatrix versions of the above


// return the value of the maximum value of the vec
double maxm(matrix I) {
    int i,j;
    //int ind1, ind2;
    double val = SHRT_MIN;
    for (i = 0; i < I.x; i++) {
        for (j = 0; j < I.y; j++) {
            if (I.m[i][j] > val) {
                val = I.m[i][j];
            }
        }
    }
    return val;
}

intmatrix makeim(int xSz, int ySz) {
    intmatrix nw_matrix;
    nw_matrix.x = xSz;
    nw_matrix.y = ySz;
    nw_matrix.m = new int*[xSz];
    int i;
    for (i = 0; i < nw_matrix.x; i++) {
        nw_matrix.m[i] = new int[ySz];
    }
    return(nw_matrix);
}

intmatrix zeroim(int xSz, int ySz) {
    int i;
    int j;
    intmatrix nw_matrix = makeim(xSz, ySz);
    for (i = 0; i < nw_matrix.x; i++) {
        for (j = 0; j < nw_matrix.y; j++) {
            nw_matrix.m[i][j] = 0;
        }
    }
    return(nw_matrix);
}

intmatrix onesim(int xSz, int ySz) {
    int i;
    int j;
    intmatrix nw_matrix = makeim(xSz, ySz);
    for (i = 0; i < nw_matrix.x; i++) {
        for (j = 0; j < nw_matrix.y; j++) {
            nw_matrix.m[i][j] = 1;
        }
    }
    return(nw_matrix);
}

intmatrix copyim(intmatrix yr_matrix) {
    int i;
    intmatrix nw_matrix = makeim(yr_matrix.x, yr_matrix.y);
    for (i = 0; i < yr_matrix.x; i++) { // NB: does not assume contiguous memory
        memcpy(nw_matrix.m[i], yr_matrix.m[i], sizeof(int) * yr_matrix.y);
    }
    return(nw_matrix);
}

matrix im2m(intmatrix yr_matrix) {
    int i;
    int j;
    matrix nw_matrix = makem(yr_matrix.x, yr_matrix.y);
    for (i = 0; i < yr_matrix.x; i++) {
        for (j = 0; j < yr_matrix.y; j++) {
            nw_matrix.m[i][j] = yr_matrix.m[i][j];
        }
    }
    return nw_matrix;
}

//copy the vec v in the row 'i' of matrix A
void copyRow(matrix* A, vec v,int i){

    int j;

    if (i >= 0 && i < A->x && v.x == A->y){

        for(j = 0;j < v.x; j++){

            A->m[i][j] = v.v[j];
        }

    }
    else{

        fprintf(stderr, "Error: the index of the row 'i' is out of bounds and/or the size of vec 'v' and the number of colums in 'A' are diferent\n");

    }
}

//copy the vec v in the row 'i' of matrix A
void copyRow2(matrix & A, vec v,int i){

    int j;

    if (i >= 0 && i < A.x && v.x == A.y){

        for(j = 0;j < v.x; j++){

            A.m[i][j] = v.v[j];
        }

    }
    else{

        fprintf(stderr, "Error: the index of the row 'i' is out of bounds and/or the size of vec 'v' and the number of colums in 'A' are diferent\n");

    }
}


//Get a row from a matrix
vec getRowm(matrix I, int row){
	vec out = zerov(I.y);
	for (int i=0; i < I.y; i++) {
	    out.v[i] = I.m[row][i];
	}

	return out;
}


matrix add(matrix a, matrix b){
	int fi,co;
	matrix c = zerom(a.x, a.y);
	for(fi = 0; fi < a.x; fi++){
		for(co = 0; co < a.y; co++){
			c.m[fi][co] = a.m[fi][co] + b.m[fi][co];
		}
	}

	return c;

}

/*vec* mat2vect(matrix a){
	vec* x = (vec*)malloc(sizeof(vec)*a.x);
	int i;
	int j;
	for(i = 0; i < a.x; i++){
		x[i] = makev(a.y);
		for(j = 0; j < a.y; j++){
			x[i].v[j] = a.m[i][j];
		}
	}
	return x;
}*/

//Apply a row-wise moving average on matrix "I"
void movAvrM(matrix & I, int sz){
    //printf("Size inside mvAvrM: %d",sz);
    vec temp = zerov(I.y);
    int cntr;
    double sum;
    for (int i = 0; i < I.x; i++) { //Row index
	temp = getRowm(I,i);
        for (int j = 0; j < I.y; j++) { //Column index
	    cntr = 0;
	    sum = 0;
	    for (int k = 0; k < sz; k++) { //m.a. index 
		if ((j + k -sz/2) >= 0 && (j + k -sz/2) < I.y) {//Check boundaries
		    sum = sum + temp.v[j + k -sz/2];
		    cntr = cntr + 1;		
		}
	    }
	    I.m[i][j] = sum  / cntr;
	}
    }

}


void freeim(intmatrix yr_matrix) {
    int i;
    for (i = 0; i < yr_matrix.x; i++) {
        delete [] yr_matrix.m[i];
    }
    delete yr_matrix.m;
}

void printim(intmatrix yr_matrix) {
    int i;
    int j;
    for (i = 0; i < yr_matrix.x; i++) {
        for (j = 0; j < yr_matrix.y; j++) {
            printf("%d\t", yr_matrix.m[i][j]);
        }
        printf("\n");
    }
}

// a naive Sieve of Erasthones for prime numbers
int sieve(intvec ones) {
    int i;
    int j;
    int k = 0;
    int sp = floor(sqrt(ones.x));
    ones.v[0] = NP; // Because 1 is not prime (though sometimes we wish it was)
    for (i = 1; i < sp; i++) {
        if PRIME(ones.v[i]) {
            for (j = i + i + 1; j < ones.x; j += i + 1) {
                ones.v[j] = NP; // Mark it not prime
            }
            k++;
        }
    }
    for (i = sp; i < ones.x; i++) { // Now we're only counting
        if PRIME(ones.v[i]) {
            k++;
        }
    }
    return(k);
}

intvec primes(int n) {
    int i;
    int j = 0;
    intvec myOnes = onesiv(n);
    intvec myPrimes = makeiv(sieve(myOnes)); // size of the # of primes
    for (i = 0; i < myOnes.x; i++) { // could start at 1, unless we're hacking
        if PRIME(myOnes.v[i]) {
            myPrimes.v[j++] = i + 1;
        }
    }
    freeiv(myOnes);
    return(myPrimes);
}


// cubic spline function, based on Numerical Recipes in C, 2nd ed.
vec spline(vec x, vec y) {
    int i;
    int j;
    double p;
    double qn;
    double sig;
    vec y2 = makev(x.x);
    double* u = new double[x.x -1];//malloc((unsigned) (x.x - 1) * sizeof(double));
    y2.v[0] = -.5; // Left boundary
    u[0] = (3. / (x.v[1] - x.v[0])) * ((y.v[1] - y.v[0]) /
                                       (x.v[1] - x.v[0]) - YP1);
    for (i = 1; i < x.x - 1; i++) { // Decomp loop
        sig = (x.v[i] - x.v[i - 1]) / (x.v[i + 1] - x.v[i - 1]);
        p = sig * y2.v[i - 1] + 2.;
        y2.v[i] = (sig - 1.) / p;
        u[i] = (y.v[i + 1] - y.v[i]) / (x.v[i + 1] - x.v[i]) -
                                  (y.v[i] - y.v[i - 1]) / (x.v[i] - x.v[i - 1]);
        u[i] = (6. * u[i] / (x.v[i + 1] - x.v[i - 1]) - sig * u[i - 1]) / p;
    }
    qn = .5; // Right boundary
    y2.v[y2.x - 1] = ((3. / (x.v[x.x - 1] - x.v[x.x - 2])) * (YPN -
                               (y.v[y.x - 1] - y.v[y.x -  2]) / (x.v[x.x - 1] -
                                    x.v[x.x - 2])) - qn * u[x.x - 2]) /
                                             (qn * y2.v[y2.x - 2] + 1.);
    for (j = x.x - 2; j >= 0; j--) { // Backsubstitution loop
        y2.v[j] = y2.v[j] * y2.v[j + 1] + u[j];
    }
    delete [] u;
    return(y2);
}

// query the cubic spline
double splinv(vec x, vec y, vec y2, double val, int hi) {
    double h;
    double b;
    double a;
    int lo = hi - 1; // find hi linearly, or using bisectv()
    h = x.v[hi] - x.v[lo];
    a = (x.v[hi] - val) / h;
    b = (val - x.v[lo]) / h;
    return(a * y.v[lo] + b * y.v[hi] + ((a * a * a - a) * y2.v[lo] *
                                  (b * b * b - b) * y2.v[hi]) * (h * h) / 6.);
}

// polynomial fitting with CLAPACK: solves poly(A, m) * X = B
vec polyfit(vec A, vec B, int order) {
    int i;
    int j;
    //int info;
    order++; // I find it intuitive this way...
    double* Ap = new double [order * A.x]; // Build up the A matrix
    for (i = 0; i < order; i++) {                      // as a vec in column-
        for (j = 0; j < A.x; j++) {                    // major-order.
            Ap[i * A.x + j] = pow(A.v[j], order - i - 1); // Mimics MATLAB
        }
    }
    vec Bp = makev(order >= B.x ? order : B.x);
    for (i = 0; i < B.x; i++) {
        Bp.v[i] = B.v[i];
    }
    i = 1; // nrhs, j is info
    j = A.x + order; // lwork
    double* work = new double[j];
    //dgels_("N", &A.x, &order, &i, Ap, &B.x, Bp.v, &order, work, &j, &info);
    delete [] Ap;
    //if (info < 0) {
    //    fprintf(stderr, "LAPACK routine dgels() returns error: %d\n", info);
    //    exit(EXIT_FAILURE);
    //}
    //else {
        return(Bp);
    //}
}

/*
*Returns the interpolated function
*@param f, the function to interpolate
*@param g
*@param x, x values corresponding to the function f
*@return y, interpolated function
*/
vec interp1(vec f, vec g, vec x){
	int i = 0;  //indice
	int a = 0;
	double x0 ;
	double x1;
	double y0;
	double y1;

	vec y = makev(x.x);

	x0 = 0;
	x1 = f.v[i];
	y0 = 0;
	y1 = g.v[i];//This section calculates when x's indexes are inferior than the ones in f
	while(x.v[a]<f.v[i] && a < x.x){
		y.v[a] = y0 + (((x.v[a]- x0)*y1 -(x.v[a] -x0)*y0)/(x1-x0));//Lineal interpolation formula
		a++;
	}
	//This section calculates when the values from x are in f's the range.
	while(i < f.x && a < x.x){
		if(f.v[i] < x.v[a]){
			i++;
		}else{
			x0 = f.v[i-1];
			x1 = f.v[i];
			y0 = g.v[i-1];
			y1 = g.v[i];
			y.v[a] = y0 + (((x.v[a]- x0)*y1 -(x.v[a] -x0)*y0)/(x1-x0));//Formula de interpolacion lineal
			a++;
		}
	}
	x0 = f.v[i-1];
	x1 = 0;
	y0 = g.v[i-1];
	y1 = 0;
	if(a < x.x){
		y.v[a] = y0 + (((x.v[a]- x0)*y1 -(x.v[a] -x0)*y0)/(x1-x0));//Lineal interpolation formula
		a++;
	}
	//This section calculates when the values of x are above the ones from f
	while(a < x.x){
		y.v[a] = 0;
		a++;
	}

	return y;
}


// given a vec of coefficients and a value for x, evaluate the polynomial
double polyval(vec coefs, double val) {
    int i;
    double sum = 0.;
    for (i = 0; i < coefs.x; i++) {
        sum += coefs.v[i] * pow(val, coefs.x  - i - 1);
    }
    return(sum);
}

// some test code
#ifdef DEBUG
int main(void) {

    int i
    int j;

    printf("vec example\n");
    vec a = makev(10);
    for (i = 0; i < a.x; i++) {
        a.v[i] = i * i;
    }
    vec b = copyv(a);
    printv(b);
    freev(a);
    freev(b);
    printf("\n");

    printf("INTvec example\n");
    intvec c = makeiv(10);
    for (i = 0; i < c.x; i++) {
        c.v[i] = i * i;
    }
    intvec d = copyiv(c);
    printiv(d);
    freeiv(c);
    freeiv(d);
    printf("\n");

    printf("more INTvec\n");
    intvec c1 = zeroiv(10);
    printiv(c1);
    freeiv(c1);
    intvec d1 = onesiv(10);
    printiv(d1);
    freeiv(d1);
    printf("\n");

    printf("MATRIX example\n");
    matrix e = makem(20, 3);
    for (i = 0; i < e.x; i++) {
        for (j = 0; j < e.y; j++) {
            e.m[i][j] = i * i + j;
        }
    }
    matrix f = copym(e);
    printm(f);
    freem(e);
    freem(f);
    printf("\n");

    printf("INTMATRIX example\n");
    intmatrix g = makeim(20, 3);
    for (i = 0; i < g.x; i++) {
        for (j = 0; j < g.y; j++) {
            g.m[i][j] = i * i + j;
        }
    }
    intmatrix h = copyim(g);
    printim(h);
    freeim(g);
    freeim(h);
    printf("\n");

    printf("SIEVE example (input: 23)\n");
    printiv(primes(23));
    printf("\n");

    printf("BILOOK example\n");
    vec fives = makev(300);
    for (i = 0; i < fives.x; i++) {
        fives.v[i] = (i + 10) * 5.;
    }
    vec twenties = makev(100);
    for (i = 0; i < twenties.x; i++) {
        twenties.v[i] = i * 20.;
    }
    printf("searching for values of vec fives in twenties...\n");
    printf("fives (sz:%d): %f < x < %f\n", fives.x, fives.v[minv(fives)],
                                                    fives.v[maxv(fives)]);
    printf("twenties (sz:%d): %f < x < %f\n", twenties.x,
                                              twenties.v[minv(twenties)],
                                              twenties.v[maxv(twenties)]);
    int hi = bisectv(twenties, fives.v[14]);
    for (i = 15; i < 30; i++) {
        hi = bilookv(twenties, fives.v[i], hi - 1);
        printf("twenties[%d] %f <= fives[%d] %f < twenties[%d] %f\n", hi - 1,
                         twenties.v[hi - 1], i, fives.v[i], hi, twenties.v[hi]);
    }
    freev(fives);
    freev(twenties);
    printf("\n");

    printf("POLY example\n");
    vec x = makev(4);
    vec y = makev(4);
    x.v[0] = 3.0;
    x.v[1] = 1.5;
    x.v[2] = 4.0;
    x.v[3] = 2.;
    y.v[0] = 2.5;
    y.v[1] = 3.1;
    y.v[2] = 2.1;
    y.v[3] = 1.0;
    printv(polyfit(x, y, 4));
    printf("\nOctave sez: -0.683446 5.276186 -10.846127 -0.092885 13.295935\n");
    printf("%f\n", polyval(polyfit(x, y, 4), 3));
    printf("Octave sez: 2.5\n\n");

}


#endif
