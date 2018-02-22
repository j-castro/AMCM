/*
 * Test file for modwt
*/

#include <omp.h>
#include <sndfile.h>
#include <iostream>
#include "wavelib/wavelib.h" //wavelib RAFAT - C wavelet library
#include <fstream>
#include "butter.h"	//Butterworth filter
#include "ts2means.h"	//dynamic clustering

using namespace std;

// return the value of the maximum value of the double * in
double maxvvACCpgi(double* in, int sz) {
    int i;
    int index = 0;
    double val = in[0];
    for (i = 1; i < sz; i++) { 
        if (in[i] > val) {
            val = in[i];
            index = i;
        }
    }
    return in[index];
}

// return the value of the manimum value of the double * in
double minvvACCpgi(double* in, int sz) {
    int i;
    int index = 0;
    double val = in[0];
    for (i = 1; i < sz; i++) {
        if (in[i] < val) {
            val = in[i];
            index = i;
        }
    }
    return in[index];
}

// Return the mean of double array "in"
#pragma acc routine vector
double meanvACCpgi(double* in, int sz) {
    
    double sum = 0;
    #pragma acc loop reduction(+:sum)
    for(int i = 0;i < sz; i++) {
	    sum = sum + in[i];
    }

    return sum / sz;
}

// Return the mean of int array "in"
//TODO parallelize
double meanvACCintpgi(int* in, int sz) {
    
    double sum = 0;
    for(int i = 0;i < sz; i++) {
	sum = sum + in[i];
    }

    return sum / sz;
}

//Get a sub-vector of array "x" from index "i1" to "i2"
//and store it in array "y"
#pragma acc routine vector
void getSubvACCpgi(double* x, int szX, int i1, int i2, double *y, int szY){

    if ((i2 >= i1) && (i2 < szX) && (i1 >= 0) && ((i2-i1+1) <= szY)) { //Parameters check
        #pragma acc loop vector
        for(int i = i1; i <= i2; i++){
            y[i-i1] = x[i];
        }
    }
}

#pragma acc routine vector
void modwtACC2pgi(int len_X, int level_N, int len_filt, double * lpd, int lpd_len, double * hpd, double *inp, double * cA, double * cD, double * filt,double * paramsOUT) {
	
  int i;
  int J;
  int temp_len;
  int iter;
  int M;
  int lenacc;

	temp_len = len_X;
	J = level_N;
	M = 1;

#pragma acc loop vector 
	for (i = 0; i < temp_len; i++) {
		paramsOUT[i] = inp[i];
	}

  lenacc = (J + 1) * temp_len; //levLens[J + 1];
	for (iter = 0; iter < J; iter++) {
		lenacc = lenacc - temp_len;
		if (iter > 0) {
			M = 2 * M;
		}

//PASTED CODE - PASTED CODE - PASTED CODE - PASTED CODE--------------
//PASTED CODE - PASTED CODE - PASTED CODE - PASTED CODE--------------
    int l;
    int i2;
    int t;
    double s;

    s = sqrt(2.0);
    #pragma acc loop vector
	  for (i2 = 0; i2 < lpd_len; i2++) {
		  filt[i2] = double(lpd[i2] / s);
		  filt[lpd_len + i2] = double(hpd[i2]  / s);
    }

#pragma acc loop vector
	  for (i2 = 0; i2 < temp_len; i2++) {
		  t = i2;
		  cA[i2] = filt[0] * paramsOUT[t]; 
		  cD[i2] = filt[lpd_len] * paramsOUT[t];
      for (l = 1; l < lpd_len; l++) {
			  t = t - M;
			  while (t >= temp_len) {
				  t = t - temp_len;
			  }
			  while (t < 0) {
				  t = t + temp_len;
			  }

			  cA[i2] =  cA[i2] + filt[l] * paramsOUT[t];
			  cD[i2] =  cD[i2] + filt[lpd_len + l] * paramsOUT[t];

      }
    }
//PASTED CODE - PASTED CODE - PASTED CODE - PASTED CODE--------------
//PASTED CODE - PASTED CODE - PASTED CODE - PASTED CODE--------------
#pragma acc loop vector
    for (i = 0; i < temp_len; i++) {
		    paramsOUT[i] = cA[i];
		    paramsOUT[lenacc + i] = cD[i];
		}
	}

}

//Compute the sample autocorrelation of vector "in" from 
//first lag "fl" to last lag "ll"
#pragma acc routine vector
void autocorrvACCpgi(double* in, int szIn, int fl, int ll, double * out, int szOut){
  
    if ((ll-fl+1) <= szOut) {

        //double mean = meanvACC(in,szIn);  	//Sample mean
        
        //To prevent OpenACC bug---------------------
        double mean = 0;
        #pragma acc loop reduction(+:mean)
        for(int i = 0;i < szIn; i++) {
	        mean = mean + in[i];
        }

        mean = mean / szIn;
        //-------------------------------------------
         
        double var = 0;		//Sample variance

        //Compute sample variance
        #pragma acc loop reduction(+:var)
        for (int i = 0; i < szIn; i++) {
	        var += (in[i]-mean)*(in[i]-mean);
        }
        var = var / szIn;

        //Compute sample-ACF
        #pragma acc loop vector
        for(int i = fl;i <= ll; i++) {// Lag loop
	        out[i-fl] = 0;

	        for(int j = i;j < szIn; j++) { //Dot product
	          out[i-fl] += (in[j] - mean) * (in[j-i] - mean); 
          }
	        out[i-fl] = out[i-fl] / (var * szIn);
        }
    }
}

#pragma acc routine vector
void powvACCpgi(double * in, int szIn, int n, double * out, int szOut) {

    if (szIn <= szOut) {
        #pragma acc loop vector
        for(int i = 0;i < szIn; i++) {
	        out[i] = pow(in[i], n);
        }
    }
}


//Get a row from a matrix
void getRowm(double ** I,int x, int y, int row, double * out){
  if (row < x) {
  	for (int i=0; i < y; i++) {
  	    out[i] = I[row][i];
  	}
  }
}

//Apply a row-wise moving average on matrix "I"
void movAvrM(double ** I,int x, int y, int sz){

  double * temp = new double[y];
  int cntr;
  double sum;
  for (int i = 0; i < x; i++) { //Row index
    getRowm(I,x,y,i, temp);
    for (int j = 0; j < y; j++) { //Column index
      cntr = 0;
      sum = 0;
	    for (int k = 0; k < sz; k++) { //m.a. index 
        if ((j + k -sz/2) >= 0 && (j + k -sz/2) < y) {//Check boundaries
          sum = sum + temp[j + k -sz/2];
          cntr = cntr + 1;		
		    }
      }
      I[i][j] = sum  / cntr;
    }
  }

}

//Copy the content of one vector into another
void copyv(double * in, int sz, double * out) {
  
  for (int i =0; i < sz; i++) {
	  out[i] = in[i];
  }
  
}

//Applies a threshold "th" to all elements in double * "in"
//and returns the result in a variable "out" 
void thresholdv(double * in, int sz, double th, double * out){

    for (int i = 0; i < sz; i++) {
      if (in[i] > th){
	      out[i] = 1;
      }
      else {
        out[i] = 0;
      }
    }
}

//Multiply v1 and v2 and store the result in v1
void dotpv(double * v1, int sz1, double * v2, int sz2){

  if (sz1 != sz2) {
    	fprintf(stderr, "ERROR: vector sizes must be the same to apply a dot product");
    }
  else {
	  for (int i=0; i<sz1; i++) {
	     v1[i] = v1[i] * v2[i];
	  }
  }
  
}

//Remove sequences of 0's and 1's shorter than "minL" in
//a binary double * "bv"
void removeSSv(double* bv,int sz, int minL){

  int seqS = 0; //Sequence start index

  //Remove small regions of zeros
  for(int i=1;i < sz; i++) {

	  if ((bv[i] - bv[i-1]) != 0) {//is a region border?
	    
	    if ((bv[i] - bv[i-1]) == -1) {//is the start of a region of zeros?
		    seqS = i;  //update the sequence start
	    }
	    else {//is the end of a region of zeros?
		    if ((i - seqS) < minL) { //check minimum length restricction
		      for (int j = seqS; j < i;j++) //Replace zeros by ones
			      bv[j] = 1;	    
		      }
	      }
	    }	
    }

    seqS = 0;
    //remove small regions of ones
    for(int i=1;i < sz; i++) {

	    if ((bv[i] - bv[i-1]) != 0) {//is a region border?
	    
	      if ((bv[i] - bv[i-1]) == 1) {//is the start of a region of ones?
		      seqS = i;  //update the sequence start
	      }
	      else {//is the end of a region of ones?
		      if ((i - seqS) < minL) { //check minimum length restricction
		        for (int j = seqS; j < i;j++) //Replace ones by zeros
			        bv[j] = 0;   
		        }
	        }
	      }	
    }

}

//Remove sequences of 1's longer than "maxL" in
//a binary double * "bv"
void removeLSv(double * bv,int sz, int maxL){

  int seqS = 0; //Sequence start index
    
    //remove long regions of ones
  for(int i=1;i < sz; i++) {

	if ((bv[i] - bv[i-1]) != 0) {//is a region border?
	    
    if ((bv[i] - bv[i-1]) == 1) {//is the start of a region of ones?
	    seqS = i;  //update the sequence start
    }
    else {//is the end of a region of ones?
	    if ((i - seqS) > maxL) { //check maximum length restricction
	      for (int j = seqS; j < i; j++) //Replace ones by zeros
		      bv[j] = 0;   
	      }
      }
    }	
  }

}

//copy the double* v in the row 'i' of matrix A
void copyRow(double** A, int xSz, int ySz,double* v, int vSz,int i){

    int j;

    if (i >= 0 && i < xSz && vSz == ySz){

        for(j = 0;j < vSz; j++){

            A[i][j] = v[j];
        }

    }
    else{

        fprintf(stderr, "Error: the index of the row 'i' is out of bounds and/or the size of vec 'v' and the number of colums in 'A' are diferent\n");

    }
}

// return the median absolute deviation for double* in, temp_vecACC is a temporary vector
#pragma acc routine vector
double madvACC(double* in, int szIn, double * temp_vecACC) {

    double temp, med, mad;
    
    #pragma acc loop vector
    for (int i=0; i < szIn; i++) {
      temp_vecACC[i] = in[i];
    }

    //Compute median of "in"
    for (int i = 0; i < szIn; i++) {	
    	for(int j = 1; j <   szIn-i; j++) {
        if ( temp_vecACC[j-1] >  temp_vecACC[j]){ //exchange positions
    		  temp =  temp_vecACC[j];
   		    temp_vecACC[j] =  temp_vecACC[j-1];
   		    temp_vecACC[j-1] = temp; 
   	    } 
    	}
    }

    med =  temp_vecACC[szIn/2];
    /*printf("Ordered vector:\n");
    printv(temp_vec);
    printf("\nMedian: %f\n",med);*/

    //Compute absolute deviation of the samples
    #pragma acc loop vector
    for (int i=0; i < szIn; i++) {
	    temp_vecACC[i] = fabs(in[i] - med);
    }
    
    /*printf("Absolute deviation vector:\n");
    printv(temp_vec);*/

    //Compute median of the absolute deviation
    for (int i = 0; i < szIn; i++) {	
	    for(int j = 1; j <  szIn-i; j++) {
	      if (temp_vecACC[j-1] > temp_vecACC[j]){ //exchange positions
		      temp = temp_vecACC[j];
		      temp_vecACC[j] = temp_vecACC[j-1];
		      temp_vecACC[j-1] = temp; 
	      } 
	    }
    }

    mad = temp_vecACC[szIn/2];
    /*printf("\nMAD IN: %f\n",mad);*/

    return mad;
}

//Multiplication of a vector by a constant
#pragma acc routine vector
void multv2ACC(double * in_vec, int sz, double alpha){   
    #pragma acc loop vector
    for(int i=0; i < sz; i++) {
        in_vec[i] = in_vec[i] * alpha;
    }
}

//Apply soft thresholding rule to vector "in" using threshold "t"
#pragma acc routine vector
void softThresvACC(double * in, int sz, double t) {    
  #pragma acc loop vector
  for (int i=0; i < sz; i++) {
	  if (fabs(in[i]) <= t)
	    in[i] = 0;
	  else if (in[i] > t)
	    in[i] = in[i] - t;
	  else 
	    in[i] = in[i] + t;
    }
}


/*INPUTS:
 *  len_X: Length of input signal X (wt->siglength)
 *  level_N: Number of levels        (wt->J;)
 *  lpd: Low Pass Details     (wt->wave->lpd)
 *  lpd_len: lpd length       (wt->wave->lpd_len)
 *  hpd: High Pass Details    (wt->wave->hpd)
 *  params: wt-object parameters  (wt->params)
 *  N: is not a parameter actually but refered here as the size of the signal in time domain
 *DYNAMIC MEMORY INPUTS:
 *  filt: a temporal filter vector 
 *  tempX:a temporary vector of size N (X)
 *OUTPUTS:
 *  output: is the vector of size N with aprox coef of first level and the output of the method (dwtop and cA)
*/
//Inverse wavelet transform TODO: match parameters with call
#pragma acc routine vector
void imodwtACC(int len_X, int level_N, double * lpd, int lpd_len, double * hpd, double * params, double * filt, double* tempX, double *output) {
	int iter, i, j, U;
	int lenacc,M;
  double * cD;
	U = 2;
	lenacc = len_X;
	M = (int)pow(2.0, (double)level_N - 1.0);

	for (iter = 0; iter < level_N; ++iter) {
		if (iter > 0) {
			M = M / 2;
		}
   
    cD = params + lenacc;  

//PASTED CODE ---------------------------------------------------------------------------------------------------------------------
  	int len_avg, i, l, t;
  	double s;
  	len_avg = lpd_len;
  
  	s = sqrt(2.0);
   #pragma acc loop vector
  	for (i = 0; i < len_avg; ++i) {
  		filt[i] = lpd[i] / s;
  		filt[len_avg + i] = hpd[i] / s;
  	}
  
    #pragma acc loop vector
  	for (i = 0; i < len_X; ++i) {
  		t = i;
  		tempX[i] = (filt[0] * output[t]) + (filt[len_avg] * cD[t]);
  		for (l = 1; l < len_avg; l++) {
  			t += M;
  			while (t >= len_X) {
  				t -= len_X;
  			}
  			while (t < 0) {
  				t += len_X;
  			}
  
  			tempX[i] += (filt[l] * output[t]) + (filt[len_avg + l] * cD[t]);
  
  		}
  	}
//PASTED CODE ---------------------------------------------------------------------------------------------------------------------
    #pragma acc loop vector
		for (j = 0; j < len_X; ++j) {
			output[j] = tempX[j];
		}

		lenacc += len_X;
	}
}


int main (int argc, char *argv[]) {
    
  double gtBeg = omp_get_wtime(); //Global time
    
  //Variables***************************************************
  double * x;                   	//noisy signal
  int xSz;                  //size of variable x
  double fs;                 	//sample frequency (Hz)
  double fc = 2000;          	//cut-off frequency (Hz)
      double tt = 80;		//transient time (milliseconds)
  double wsz = 3;             //window size (milliseconds)
  int ln = 4;               	//number of levels
 	double energy;			        //signal energy computed at each level/window
 	double univ;			//universal threshold
  double startP, endP, durationP; //time 
  //***********************************************************
     
  //Open a WAV*************************************************
  FILE* wavf;
  SF_INFO info;
  char* fileName;
  
  if (argc == 2) {
    fileName = argv[1];
	  printf("\n<<<Opening file: %s>>>\n\n",argv[1]);
  }
  else {
    fileName = "aud/9manatees.wav";
    printf("\n<<<Using default file: aud/9manatees.wav>>>\n\n");
  }

  wavf = fopen(fileName, "r");
    
  if(wavf != NULL){   //Succesful wav read

    SNDFILE* source = sf_open_fd(fileno(wavf), SFM_READ, &info, 1);	

    fs = info.samplerate;
    xSz = info.frames;
    x = new double[xSz];
    sf_read_double(source, x, xSz);
    sf_close(source);
    fclose(wavf);

    //Denoising method begins*****************************************************        
    //printf("\nExecuting parallel wavelet... \n\n");		
	      
    //Process constants-------------------
    wave_object family = wave_init("db8");  //wavelet family
    double ZERO_CONST =  pow(10,-14);       //represents minimum signal energy
   	double AC_MIN_CONST = 0.1;             	//call minimun autocorrelation-rms value
   	double MAN_MAX_DUR = 1 * fs;           	//call maximun duration in # of samples
    double wszNS = llround(wsz*0.001*fs);  	//window size in # samples
    double ttNS = llround(tt*0.001 * fs);  	//transients duration in # samples
    double maSize = llround(ttNS / wszNS);   //moving average order (# of windows)(22)
    int wszINT = wszNS;
    int minVD = llround(ttNS / wszNS);       	//minimun vocalization duration (# of windows) (25)
	  int maxVD = llround(MAN_MAX_DUR / wszNS);   	//maximun vocalization duration (# of windows)	
    //--------------------------------------
    
   	//Initialize Variables-------------------------------------------------
    double * y = new double[xSz];	//output signal
    double alagBegNS = round(0.21 *0.001*fs);     //autocorrelation lag beggining
    double alagEndNS = round(1.25 *0.001*fs);     //autoorrelation lag end
 		double * inter_vec = new double[xSz];	//intermediate output vector
    double dt = 1. / fs;		//time step in seconds
  	double lastIdx =0;		//last index of processed portion of the signal
  	wt_object wt = wt_init(family,"modwt", wszNS, ln); //undecimated wavelet transform object
  
  	//DEBUG - Print basic info -----------------------------------
  	/*printf("\nProcess constants ----------------- \n");
    printf("Signal length: %d\n", xSz);
  	printf("wszNS: %g samples\n", wszNS);
    printf("Alag in number of smaples: [%f,%f]\n",alagBegNS,alagEndNS);
  	printf("----------------------------------- \n\n");*/
  	//--------------------------------------------------------------------	

  	//Allocate all memory before OpenACC region
    
    //Size of vectors in number of samples 	
  	int svSz = wszINT;     //size of temporal subvector
    int wsvSz = wszINT*(ln+1); //size of wavelet coefficients of the temporal subvector
  	int wvSz = xSz*(ln+1); //size of wavelet coefficients of the entire signal
    int avSz = alagEndNS - alagBegNS + 1; //size of autocorrelation vector
    int tfbSz = wszINT;
    
    //Size of blocks and windows
    int bSz = 128; //block size in number of windows
    int nWs = floor(xSz/wszINT);  //number of windows
    int nBlocks = ceil( floor(xSz/wszINT) / bSz );//number of blocks 
    int resWsz = nWs%bSz; //residual of the # of windows divided by the block size
        
    
  	//OpenACC variables-------------------------------------------------------------------------------------
  	double * restrict sub_vecACC = new double[svSz]; //a temporal subvector of the input signal      
  	double * restrict inter_vecACC = new double[xSz]; //= x; //a pointer to the original input //TODO RESTORE!!
    double * restrict w_vecACC = new double[wvSz]; //UDWT of the entire signal
    double * restrict tfbACC = new double[tfbSz];   //a temporal vector to store one frequency channel of sub-vector sub_vecACC
    double * restrict ACF_vecACC = new double[avSz] ;        //vector that store the autocorrelation function result         
    //printf("\nDebug checkpoint #0\n");
    double ** restrict rms_vecACC = new double*[ln+1];       //rms value matrix of ACF at each level of each time window
    //------------------------------------------------------------------------------------------------------
     
    //High Pass filter-------------------------------------------------------------------------------
    int fOrd = 4;		//filter order
  	double * bb = (double *)calloc( fOrd+1, sizeof(double) );	//TF numerator
  	double * aa = (double *)calloc( fOrd+1, sizeof(double) );	//TF denominator
  	
  	butterHP(fOrd, fc/(fs/2),bb ,aa); 		//Get butterworth high pass filter coefficients
  	filterS(fOrd,bb,aa,xSz,x,inter_vecACC); 	//apply high-pass filter INITIALIZATION OF INTER_VECACC???
  
  	free(bb);
  	free(aa);    
    //----------------------------------------------------------------------------------------------- 
     
    	//printf("\nDebug checkpoint #1\n");
    //Initialization of matrix  
    for(int i = 0; i < ln+1; i++){
      rms_vecACC[i] = new double[nWs];
      /*for (int j = 0; j < nWs; j++) { Not necessary
        rms_vecACC[i][j] = 0;
      }*/
    }
     
    //Initialization of wavelet coefficient's vector 
    for (int yy = 0; yy < xSz*(ln+1); yy++){
      w_vecACC[yy] = 0;
    }

    //Input parameters of modwtACC2-------------------------------------
    int len_X = wt->siglength;          //Length of input signal X
    int level_N = wt->J;                //Number of levels
    int len_filt = wt->wave->filtlength;//Wave-filter length     
    int * levLens = wt->length;         //Length of each level
    double * lpd = wt->wave->lpd;       //Low Pass Details   
    int lpd_len = wt->wave->lpd_len;    //lpd length     
    double * hpd = wt->wave->hpd;       //High Pass Details   
    //Dynamic memory inputs: allocate and free memory outside ACC-loop 
    double * cA = new double[len_X]; //free(cA);
    double * cD = new double[len_X]; //free(cD);
    double * filt = new double[2 * lpd_len]; //free(filt)
    //Output parameters:
    double * paramsOUT = new double[wsvSz];       
    //------------------------------------------------------------------

    
    //----------------------------------------------------------------
    /*printf("\nBEFORE ACC REGION ----------\n");
    printf("Number of windows: %d\n",nWs);
    printf("Number of blocks: %d\n",nBlocks);
    printf("Residual of windows: %d\n",resWsz);
    printf("w_vecACC size: %d\n", xSz*(ln+1));
    printf("w_vecACC update operation's size %d\n", (xSz - wszINT)*(ln+1));
    printf("lpd length: %d\n", lpd_len);
    printf("LPD and HPD mean values: (%f,%f)\n", meanvACC(lpd,lpd_len),meanvACC(hpd,lpd_len));
    printf("LPD and HPD max: (%f,%f)\n",maxvvACC(lpd,lpd_len),maxvvACC(hpd,lpd_len));
    printf("inter_vecACC max value is: %f\n", maxvvACC(inter_vecACC,xSz));
    printf("inter_vecACC mean value is: %f\n", meanvACC(inter_vecACC,xSz));*/
    
    //Set levLens before cycle-------------------------
    levLens[level_N +1]= (level_N +1)*len_X;
    for (int i = 0; i < level_N+1; i++) {
        levLens[i] = len_X;
    }   
    //--------------------------------------------------       
    //printf("levLens mean: %f\n", meanvACCint(levLens,level_N + 2));        
    //----------------------------------------------------------------------------


    //DEBUG - READ ORIGINAL WAVELET COEFFICIENTS----------------------------------------
    ifstream myFile; //DEBUG CODE 
	      myFile.open ("manateeOUT.txt"); //DEBUG CODE 
    int lineCnt = 0;
    string line;
    while (getline(myFile,line))
        ++lineCnt;
    myFile.clear();
    myFile.seekg(0,ios::beg);

    double * fileCoef = new double[lineCnt];
    for (int i=0; i <  lineCnt; i++) {
         myFile >>  fileCoef[i]; //DEBUG CODE
    }
	      myFile.close(); //DEBUG CODE
    
    //----------------------------------------------------------------------------------
    // DEBUG - WRITE WAVELET COEFFICIENTS-----------------------------------------------
   	ofstream myfile2;
    myfile2.open ("manateeOUT2.txt"); //DEBUG CODE*/
    startP = omp_get_wtime();
    //----------------------------------------------------------------------------------

//1 - Wavelet-ACF stage------------------------------------------------------------------------------------------------------------------------------
#pragma acc enter data create(w_vecACC[:wvSz]) 
#pragma acc data copyin(hpd[:lpd_len],lpd[:lpd_len]), copyout(rms_vecACC[:ln+1][:nWs]), create(inter_vecACC[:xSz]), present(w_vecACC) //**ASYNC COPYIN**          
{
    for (int q = 0; q < nBlocks; q++) { //For each block of windows
      int begI = q * bSz * wszINT; //block beggining
      int endI = min(begI + (bSz * wszINT),((nBlocks-1) * bSz * wszINT) + (resWsz * wszINT)); //block end

#pragma acc update device(inter_vecACC[begI:endI-begI]) async (q%4) //**ASYNC COPYIN**             
#pragma acc parallel loop gang private(sub_vecACC[:svSz],paramsOUT[:wsvSz],cA[:len_X],cD[:len_X],filt[:2*lpd_len]) async (q%4) //gang, vector_length(128) //**ASYNC COPYIN**          
      for (int i = begI ;i < endI; i = i + wszINT) { //For each window compute the ACF-rms value in the wavelet domain
      		    
        int val = i + wszINT -1; //added to avoid a pgi compiler's bug
        //Get subvector	(window)	
  		  getSubvACCpgi(inter_vecACC, xSz,i,val,sub_vecACC, svSz);
  
  		  //Apply UDWT to window
        modwtACC2pgi(len_X,level_N,len_filt,lpd,lpd_len,hpd,sub_vecACC,cA,cD,filt,paramsOUT);
                      
  		  //Store UDWT coefficients
        #pragma acc loop
        for (int j=0; j < wsvSz; j++) {
  	      w_vecACC[i*(ln+1) + j] = paramsOUT[j];
        }
                 
        //ACF and rms value computation for each level
        #pragma acc loop private(tfbACC[:tfbSz],ACF_vecACC[:avSz],energy,i)           
        for (int j=0;j < ln+1; j++) {
  		    
  	      energy = 0; 
    			
          #pragma acc loop reduction(+:energy)
    		  for (int k=0; k < wszINT; k++) { //Obtain level j's coefs and estimate signal energy		
    			  tfbACC[k] = paramsOUT[k + (j*wszINT)];
    			  energy += pow(tfbACC[k],2); 
    		  }         
    		  energy = sqrt(energy / wszNS);  //signal energy in window "i" and level "j"
  
    		  if (energy < ZERO_CONST) { //if minimun energy is not reached, asign zero score
    			  rms_vecACC[j][i/wszINT] = 0; //zero score to the actual window and level 
    		  }
    		  else { 					                                                
            autocorrvACCpgi(tfbACC,tfbSz,alagBegNS,alagEndNS, ACF_vecACC,avSz); //Autocorrelation function
            powvACCpgi(ACF_vecACC,avSz,2,ACF_vecACC,avSz);		
            rms_vecACC[j][i/wszINT] = sqrt(meanvACCpgi(ACF_vecACC,avSz)); //store rms value of ACF on [alag[0] alag[1]] range                                   
    		  }
        }
      }     				                                  
          //Update block of data
      #pragma acc update self(w_vecACC[begI*(ln+1):(endI-begI)*(ln+1)]) async (q%4) //,rms_vecACC[:ln+1][begI/wszINT:(endI-begI)/wszINT]
    }
    #pragma acc wait //ASYNC wait 
}
//----------------------------------------------------------------------------------------------------------------------------------------------------
  
   	endP = omp_get_wtime();
    durationP = endP - startP;
    
    printf("Parallel execution time of OpenACC region #1: %f seconds\n", durationP);                            
    /*printf("LPD and HPD mean values: (%f,%f)\n", meanvACC(lpd,lpd_len),meanvACC(hpd,lpd_len));
    printf("LPD and HPD max: (%f,%f)\n",maxvvACC(lpd,lpd_len),maxvvACC(hpd,lpd_len));
    printf("inter_vecACC max value is: %f\n", maxvvACC(inter_vecACC,xSz));
    printf("inter_vecACC mean value is: %f\n", meanvACC(inter_vecACC,xSz));
    printf("levLens mean: %f\n", meanvACCint(levLens,level_N + 2));            

    printf("w_vecACC max value is: %f\n", maxvvACC(w_vecACC, wvSz));
    printf("w_vecACC min value is: %f\n", minvvACC(w_vecACC, wvSz));	    
    printf("w_vecACC mean value is: %f\n", meanvACC(w_vecACC, wvSz));*/
    //printf("line_Cnt is: %d and rms_vecACC size is: %d\n",lineCnt, nWs*(ln+1));         
    
    startP = omp_get_wtime();
    //2 - Dynamic clustering and removal of transients stage-------------------------------------------------------------------------  	    
    movAvrM(rms_vecACC,ln+1,nWs, maSize); //transient attenuation using a moving average filter
    double * temp_vecACC1 = new double[nWs];
    double * temp_vecACC2 = new double[nWs];
        
    for (int j=0;j < ln+1; j++) {

  		getRowm(rms_vecACC,ln+1,nWs,j,temp_vecACC1); //extract values for level j
  		copyv(temp_vecACC1,nWs,temp_vecACC2);        //backup the extracted row
  		ts2means4(temp_vecACC1, nWs,dt,temp_vecACC1); //cluster the rms values of the whole signal into two clases 0 and 1		
  		thresholdv(temp_vecACC2,nWs,AC_MIN_CONST,temp_vecACC2);//Apply a minimun threshold to sub_vec2
  		dotpv(temp_vecACC1,nWs,temp_vecACC2,nWs); //multiply sub_vec and subvec2 and store in subvec		
  		removeSSv(temp_vecACC1,nWs, minVD); //Remove short-duration sequences of 1's or 0's less than minVD
  		removeLSv(temp_vecACC1,nWs, maxVD); //Remove long-duration sequences of 1's      		
  		copyRow(rms_vecACC,ln+1,nWs,temp_vecACC1,nWs,j); //store "" in row "i" of matrix "C"		

    }

    //Delete memory for temp vectors
    delete [] temp_vecACC1;
    delete [] temp_vecACC2;
    //--------------------------------------------------------------------------------------------------------------------------------
    endP = omp_get_wtime();
    durationP = endP - startP;
    printf("Parallel execution time of OpenACC region #2: %f seconds\n", durationP);
    
    //TODO: parallelize and Debug, run asynchronously, quit initialization of variables in loop
    int val1, val2;
    double * temp_vecACC3 = new double[svSz];
    double* tempX = new double[len_X];
     
    startP = omp_get_wtime();
    //3- Wavelet thresholding ----------------------------------------------------------------------------------------------------------
#pragma acc data copyin(hpd[:lpd_len],lpd[:lpd_len],rms_vecACC[:ln+1][:nWs]), copyout(y[:xSz]), present(w_vecACC) //**ASYNC COPYIN**          
{    	
    for (int q = 0; q < nBlocks; q++) { //For each block of windows
      int begI = q * bSz * wszINT; //block beggining
      int endI = min(begI + (bSz * wszINT),((nBlocks-1) * bSz * wszINT) + (resWsz * wszINT)); //block end
#pragma acc parallel loop gang private(sub_vecACC[:svSz],paramsOUT[:wsvSz],temp_vecACC3[:svSz],filt[:2*lpd_len],tempX[len_X])        
      for (int i = begI ;i < endI; i = i + wszINT) { //For each window compute the ACF-rms value in the wavelet domain    
	    
        val1 = i*(ln+1);           //added to avoid a pgi compiler's bug
        val2 = (i+wszINT)*(ln+1)-1; //added to avoid a pgi compiler's bug
        
  	    //3.1 Get UDWT coefficients
   		  getSubvACCpgi(w_vecACC, wvSz,val1,val2,paramsOUT, wsvSz);       
  		
        val1 = ln*wszINT;
        val2 = ln*wszINT + wszINT-1;
   	    
         //3.2 Compute Universal threshold
        getSubvACCpgi(paramsOUT, wsvSz,val1,val2,sub_vecACC, svSz); //Obtain highest level coefficients                  	       
        univ = madvACC(sub_vecACC, svSz, temp_vecACC3) * sqrt(2*log(wszNS)); // 0.6745 decide if using this factor
  		
  	    //3.3 Apply soft thresholding on wavelet coefficients
  	    for (int j=0; j < ln+1; j++) { //for each level
  		
  		    //get all wavelet coef. from level j
          #pragma acc loop 
  		    for (int k=0; k < wszINT; k++)
  		      sub_vecACC[k] = paramsOUT[k + (j*wszINT)];
  				
  		    //Determine threshold to apply
  		    if (rms_vecACC[j][i/wszINT] == 0) {//if belongs to class 0
  		      multv2ACC(sub_vecACC,svSz,0); //set all coefficients to zero (infinite threshold)
  		    }
  		    else {
            softThresvACC(sub_vecACC,svSz,univ); //apply soft thresholding rule using universal threshold
  		    }
  
  		    //Copyback all thresholded wavelet coef. from level j to sub_vector
          #pragma acc loop    
  		    for (int k=0; k < wszINT; k++)
  		      paramsOUT[k + (j*wszINT)] = sub_vecACC[k];
  
  	    }		
     
  	    //3.4 Apply Inverse UDWT
        #pragma acc loop   
        for (int j=0; j < svSz; j++) { //Copy first level coefficients to sub_vecACC (Initialization for imodwtACC)
  		    sub_vecACC[j] = paramsOUT[j];
  	    }
  	    
        imodwtACC(len_X,level_N, lpd,lpd_len, hpd, paramsOUT, filt, tempX, sub_vecACC); //Inverse UDWT		
  	    
  	    //concatenate result to output vector
        #pragma acc loop   
  	    for (int j=0; j< svSz;j++) {
  		    y[i+j] = sub_vecACC[j];
  	    }
  	    
  	    //lastIdx = i + svSz;
      }
	  }
}
	//-----------------------------------------------------------------------------------------------------------------------------------------------
 #pragma acc exit data delete(w_vecACC)
     
    lastIdx = ((nBlocks-1) * bSz * wszINT) + (resWsz * wszINT);
	  //Get last signal samples
	  for (int i=lastIdx; i < xSz; i++) {
	    y[i] = inter_vecACC[i];
	  }  
    
    endP = omp_get_wtime();
    durationP = endP - startP;
    printf("Parallel execution time of OpenACC region #3: %f seconds\n", durationP);
        
    //Compare wavelet coefficients----------------
    double mae = 0;
    double sigSum = 0;
    double sigSum2 = 0;
    for (int i = 0; i < xSz; i++) { //Comparison for 1 dimensional vectors
       myfile2 << y[i] << "\n"; 
       mae += fabs(fileCoef[i] - y[i]);
       sigSum += fabs(fileCoef[i]);
       sigSum2 += fabs(y[i]);
    }
    /*for (int j = 0;j<ln+1; j++) {
      for (int i = 0;i<nWs; i++) { //Comparison for 2 dimensional vectors
        myfile2 << rms_vecACC[j][i] << "\n"; 
        mae += fabs(fileCoef[j*nWs + i] -  rms_vecACC[j][i]);
        sigSum += fabs(fileCoef[j*nWs + i]);
        sigSum2 += fabs(rms_vecACC[j][i]);
      }
    }*/
    
    sigSum = sigSum / lineCnt;
    sigSum2 = sigSum2 / lineCnt;
    mae = mae / lineCnt;
    printf("MAE: %f\n", mae);
    printf("Original SigSum: %f\n", sigSum);
    printf("Obtained SigSum: %f\n", sigSum2);
    //--------------------------------------------
    myfile2.close();
    
    
    delete [] fileCoef;// */
    
    //write wavfile with results
    SNDFILE *outF = sf_open("pgimexaParResult.wav",SFM_WRITE, &info);       	
	  sf_write_double(outF, y, xSz);
	  sf_close (outF);
    
    
    //Free used memory in ACC cycle
    //Delete rms structure
    for (int i = 0;i < ln+1; i++){
      delete [] rms_vecACC[i];
    }
    delete [] rms_vecACC;
    delete [] temp_vecACC3;
    delete [] tempX;
    delete [] tfbACC;
    delete [] inter_vec;
    delete [] sub_vecACC;
    delete [] paramsOUT; 
    delete [] cA;
    delete [] cD;
    delete [] filt;
    delete [] w_vecACC; 
    delete [] y;
    wt_free(wt);
           
    //---------------------------------------------------------------------------
    //---------------------------------------------------------------------------

    delete[] x;

  }
  double gtEnd = omp_get_wtime(); //Global time
  double gtDur = gtEnd - gtBeg;
  printf("Global execution time: %f seconds\n", gtDur);                            
     
  return 0;
}
