/*
 *                            COPYRIGHT
 *
 *  liir - Recursive digital filter functions
 *  Copyright (C) 2007 Exstrom Laboratories LLC
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  A copy of the GNU General Public License is available on the internet at:
 *
 *  http://www.gnu.org/copyleft/gpl.html
 *
 *  or you can write to:
 *
 *  The Free Software Foundation, Inc.
 *  675 Mass Ave
 *  Cambridge, MA 02139, USA
 *
 *  You can contact Exstrom Laboratories LLC via Email at:
 *
 *  stefan(AT)exstrom.com
 *
 *  or you can write to:
 *
 *  Exstrom Laboratories LLC
 *  P.O. Box 7651
 *  Longmont, CO 80501, USA
 *
 */

#include "butter.h"

/**********************************************************************
  binomial_mult - multiplies a series of binomials together and returns
  the coefficients of the resulting polynomial.
  
  The multiplication has the following form:
  
  (x+p[0])*(x+p[1])*...*(x+p[n-1])

  The p[i] coefficients are assumed to be complex and are passed to the 
  function as a pointer to an array of doubles of length 2n.

  The resulting polynomial has the following form:
  
  x^n + a[0]*x^n-1 + a[1]*x^n-2 + ... +a[n-2]*x + a[n-1]
  
  The a[i] coefficients can in general be complex but should in most
  cases turn out to be real. The a[i] coefficients are returned by the
  function as a pointer to an array of doubles of length 2n. Storage
  for the array is allocated by the function and should be freed by the
  calling program when no longer needed.
  
  Function arguments:
  
  n  -  The number of binomials to multiply
  p  -  Pointer to an array of doubles where p[2i] (i=0...n-1) is
        assumed to be the real part of the coefficient of the ith binomial
        and p[2i+1] is assumed to be the imaginary part. The overall size
        of the array is then 2n.
*/

double *binomial_mult( int n, double *p )
{
    int i, j;
    double *a;

    a = (double *)calloc( 2 * n, sizeof(double) );
    if( a == NULL ) return( NULL );

    for( i = 0; i < n; ++i )
    {
	for( j = i; j > 0; --j )
	{
	    a[2*j] += p[2*i] * a[2*(j-1)] - p[2*i+1] * a[2*(j-1)+1];
	    a[2*j+1] += p[2*i] * a[2*(j-1)+1] + p[2*i+1] * a[2*(j-1)];
	}
	a[0] += p[2*i];
	a[1] += p[2*i+1];
    }
    return( a );
}


/**********************************************************************
  trinomial_mult - multiplies a series of trinomials together and returns
  the coefficients of the resulting polynomial.
  
  The multiplication has the following form:

  (x^2 + b[0]x + c[0])*(x^2 + b[1]x + c[1])*...*(x^2 + b[n-1]x + c[n-1])

  The b[i] and c[i] coefficients are assumed to be complex and are passed
  to the function as a pointers to arrays of doubles of length 2n. The real
  part of the coefficients are stored in the even numbered elements of the
  array and the imaginary parts are stored in the odd numbered elements.

  The resulting polynomial has the following form:
  
  x^2n + a[0]*x^2n-1 + a[1]*x^2n-2 + ... +a[2n-2]*x + a[2n-1]
  
  The a[i] coefficients can in general be complex but should in most cases
  turn out to be real. The a[i] coefficients are returned by the function as
  a pointer to an array of doubles of length 4n. The real and imaginary
  parts are stored, respectively, in the even and odd elements of the array.
  Storage for the array is allocated by the function and should be freed by
  the calling program when no longer needed.
  
  Function arguments:
  
  n  -  The number of trinomials to multiply
  b  -  Pointer to an array of doubles of length 2n.
  c  -  Pointer to an array of doubles of length 2n.
*/

double *trinomial_mult( int n, double *b, double *c )
{
    int i, j;
    double *a;

    a = (double *)calloc( 4 * n, sizeof(double) );
    if( a == NULL ) return( NULL );

    a[2] = c[0];
    a[3] = c[1];
    a[0] = b[0];
    a[1] = b[1];
  
    for( i = 1; i < n; ++i )
    {
	a[2*(2*i+1)]   += c[2*i]*a[2*(2*i-1)]   - c[2*i+1]*a[2*(2*i-1)+1];
	a[2*(2*i+1)+1] += c[2*i]*a[2*(2*i-1)+1] + c[2*i+1]*a[2*(2*i-1)];

	for( j = 2*i; j > 1; --j )
	{
	    a[2*j]   += b[2*i] * a[2*(j-1)]   - b[2*i+1] * a[2*(j-1)+1] + 
		c[2*i] * a[2*(j-2)]   - c[2*i+1] * a[2*(j-2)+1];
	    a[2*j+1] += b[2*i] * a[2*(j-1)+1] + b[2*i+1] * a[2*(j-1)] +
		c[2*i] * a[2*(j-2)+1] + c[2*i+1] * a[2*(j-2)];
	}

	a[2] += b[2*i] * a[0] - b[2*i+1] * a[1] + c[2*i];
	a[3] += b[2*i] * a[1] + b[2*i+1] * a[0] + c[2*i+1];
	a[0] += b[2*i];
	a[1] += b[2*i+1];
    }

    return( a );
}


/**********************************************************************
  dcof_bwlp - calculates the d coefficients for a butterworth lowpass 
  filter. The coefficients are returned as an array of doubles.

*/

double *dcof_bwlp( int n, double fcf )
{
    int k;            // loop variables
    double theta;     // M_PI * fcf / 2.0
    double st;        // sine of theta
    double ct;        // cosine of theta
    double parg;      // pole angle
    double sparg;     // sine of the pole angle
    double cparg;     // cosine of the pole angle
    double a;         // workspace variable
    double *rcof;     // binomial coefficients
    double *dcof;     // dk coefficients

    rcof = (double *)calloc( 2 * n, sizeof(double) );
    if( rcof == NULL ) return( NULL );

    theta = M_PI * fcf;
    st = sin(theta);
    ct = cos(theta);

    for( k = 0; k < n; ++k )
    {
	parg = M_PI * (double)(2*k+1)/(double)(2*n);
	sparg = sin(parg);
	cparg = cos(parg);
	a = 1.0 + st*sparg;
	rcof[2*k] = -ct/a;
	rcof[2*k+1] = -st*cparg/a;
    }

    dcof = binomial_mult( n, rcof );
    free( rcof );

    dcof[1] = dcof[0];
    dcof[0] = 1.0;
    for( k = 3; k <= n; ++k )
        dcof[k] = dcof[2*k-2];
    return( dcof );
}

/**********************************************************************
  dcof_bwhp - calculates the d coefficients for a butterworth highpass 
  filter. The coefficients are returned as an array of doubles.

*/

double *dcof_bwhp( int n, double fcf )
{
    return( dcof_bwlp( n, fcf ) );
}


/**********************************************************************
  ccof_bwlp - calculates the c coefficients for a butterworth lowpass 
  filter. The coefficients are returned as an array of integers.

*/

double *ncof_bwlp( int n )
{
    double *ccof;
    int m;
    int i;

    ccof = (double *)calloc( n+1, sizeof(double) );
    if( ccof == NULL ) return( NULL );

    ccof[0] = 1;
    ccof[1] = n;
    m =  n/2;
    for( i=2; i <= m; ++i)
    {
        ccof[i] = (double)((n-i+1)*((int)ccof[i-1]/i));
        ccof[n-i]= ccof[i];
    }
    ccof[n-1] = n;
    ccof[n] = 1;

    return( ccof );
}

/**********************************************************************
  ccof_bwhp - calculates the c coefficients for a butterworth highpass 
  filter. The coefficients are returned as an array of integers.

*/

double *ncof_bwhp( int n )
{
    double *ccof;
    int i;

    ccof = ncof_bwlp( n );
    if( ccof == NULL ) return( NULL );

    for( i = 0; i <= n; ++i)
        if( i % 2 ) ccof[i] = -ccof[i];

    return( ccof );
}


/**********************************************************************
  sf_bwhp - calculates the scaling factor for a butterworth highpass filter.
  The scaling factor is what the c coefficients must be multiplied by so
  that the filter response has a maximum value of 1.

*/

double sf_bwhp( int n, double fcf )
{
    int m, k;         // loop variables
    double omega;     // M_PI * fcf
    double fomega;    // function of omega
    double parg0;     // zeroth pole angle
    double sf;        // scaling factor

    omega = M_PI * fcf;
    fomega = sin(omega);
    parg0 = M_PI / (double)(2*n);

    m = n / 2;
    sf = 1.0;
    for( k = 0; k < n/2; ++k )
        sf *= 1.0 + fomega * sin((double)(2*k+1)*parg0);

    fomega = cos(omega / 2.0);

    if( n % 2 ) sf *= fomega + sin(omega / 2.0);
    sf = pow( fomega, n ) / sf;

    return(sf);
}

void filterS(int ord, double *b,  double *a, int np, double *x, double *y)
{
    int i,j;
    y[0]=b[0] * x[0];
    for (i=1;i<ord+1;i++)
    {
        y[i]=0.0;
        for (j=0;j<i+1;j++)
            y[i]=y[i]+b[j]*x[i-j];
        for (j=0;j<i;j++)
            y[i]=y[i]-a[j+1]*y[i-j-1];
    }
    for (i=ord+1;i<np;i++)
    {
        y[i]=0.0;
        for (j=0;j<ord+1;j++)
            y[i]=y[i]+b[j]*x[i-j];
        for (j=0;j<ord;j++)
            y[i]=y[i]-a[j+1]*y[i-j-1];
    }
   
}

void butterHP(int filterOrder, double cutoff, double * NumC, double * DenC) {

	double * NC = ncof_bwhp( filterOrder );
	double * DC = dcof_bwhp( filterOrder, cutoff );
	double normC = sf_bwhp(filterOrder , cutoff); //normalization coefficient
	
	/*printf("Filter Coefficients\nNum \t\t Den\n");*/
	for(int k = 0; k <= filterOrder; k++)
    	{
		NumC[k] = NC[k] * normC;
		DenC[k] = DC[k];		
		/*printf("%f \t %f \n",NumC[k],DenC[k]);*/
    	}
	
}
