/*********************************************************************

    File:           spline.c
    
    Created:        Wed Nov 11 14:24:14 EST 1992
    
    Author:         Stephen Walker
                    CSIRO Division of Oceanography
   
    Purpose:	    Spline fitting routine similar to those
		    in "Numerical Recipes in C" by Press et al.
		    Here more error checking is done and arrays
		    are zero based.
    
    Revisions:      none

    $Id: spline.c 2761 2011-07-01 04:35:25Z gor171 $

*********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include "sjwlib.h"


/** Calculate the second derivative.
  * Taken from "Numerical Recipes in C" by Press et al.
  * Here more error checking is done and arrays
  * are zero based.
  * @param x input x values
  * @param y input y values
  * @param n number of input points
  * @param derivspec flag true if end derivatives specified
  * @param start_deriv first derivative at start
  * @param end_deriv first derivative at end
  * @param ydd poiinter to returned second derivative values
  */
void
dspline(
double *x,		/* input x values */
double *y,		/* input y values */
long int n,		/* number of input points */
int derivspec,		/* flag true if end derivatives specified */
double start_deriv,	/* first derivative at start */
double end_deriv,	/* first derivative at end */
double *ydd		/* output second derivative values */
)
{
	long i,k;
	double p,qn,sig,un;
	double *u;
    
	if( n < 3 )
	    quit("dspline: not enough points\n");
	for(i=0; i<n-1; i++)
	    if( x[i+1] <= x[i] )
		quit("dspline: x values not monotonic increasing\n");

	if( (u=(double *)malloc((size_t)(n-1)*sizeof(double))) == NULL )
	    quit("dspline: Can't allocate memory\n");
	if( !derivspec )
		ydd[0]=u[0]=0.0;
	else {
		ydd[0] = -0.5;
		u[0]=(3.0/(x[1]-x[0]))*((y[1]-y[0])/(x[1]-x[0])-start_deriv);
	}
	for (i=1;i<n-1;i++) {
		sig=(x[i]-x[i-1])/(x[i+1]-x[i-1]);
		p=sig*ydd[i-1]+2.0;
		ydd[i]=(sig-1.0)/p;
		u[i]=(y[i+1]-y[i])/(x[i+1]-x[i]) - (y[i]-y[i-1])/(x[i]-x[i-1]);
		u[i]=(6.0*u[i]/(x[i+1]-x[i-1])-sig*u[i-1])/p;
	}
	if( !derivspec )
		qn=un=0.0;
	else {
		qn=0.5;
		un=(3.0/(x[n-1]-x[n-2]))*(end_deriv-(y[n-1]-y[n-2])/(x[n-1]-x[n-2]));
	}
	ydd[n-1]=(un-qn*u[n-2])/(qn*ydd[n-2]+1.0);
	for (k=n-2;k>=0;k--)
		ydd[k]=ydd[k]*ydd[k+1]+u[k];
	free(u);
}


/** Calculate spline'd value at a known x position.
  * Taken from "Numerical Recipes in C" by Press et al.
  * Here more error checking is done and arrays
  * are zero based.
  *
  * @param xa original x values
  * @param ya original y values
  * @param ydd second derivative values calculated from dspline
  * @param n number of input points
  * @param x x value at evaluation point
  * @param y pointer to returned y value
  * @see dspline
  */
void
dsplint(
double *xa,		/* original x values */
double *ya,		/* original y values */
double *ydd,		/* second derivative values calculated from above */
long int n,		/* number of input points */
double x,		/* x value at evaluation point */
double *y		/* pointer to returned y value */
)
{
	int lo,hi,k;
	double h,b,a;
    
	if( n < 3 )
	    quit("dsplint: not enough points\n");
	    
	lo=0;
	hi=(int)n-1;
	if( x < xa[lo] || x > xa[hi] )
	    quit("dsplint: evaluation point outside original data\n");

	while (hi-lo > 1) {
		k=(hi+lo) >> 1;
		if (xa[k] > x) hi=k;
		else lo=k;
	}
	h=xa[hi]-xa[lo];
	if (h == 0.0)
	    quit("splint: bad input x data\n");
	a=(xa[hi]-x)/h;
	b=(x-xa[lo])/h;
	*y=a*ya[lo]+b*ya[hi]+((a*a*a-a)*ydd[lo]+(b*b*b-b)*ydd[hi])*(h*h)/6.0;
}
