/*

    File:           weight_fn.c
    
    Created:        Thu Jun  1 09:15:04 EST 1995
    
    Author:         Stephen Walker
                    CSIRO Division of Oceanography
    
    Purpose:        This file contains routines to
		    calculate various functional forms.
		    All routines are of the form

		    double function(x,scale);

                    All functions have unit value at x=0,
		    and in general decrease to zero
		    away from x=0 over a distance of scale.

    Arguments:      x - where to evaluate the function
		    scale - length scale over which function
			    decreases to zero

    Returns:        double     (1.0 at x=0)
    
    Revisions:      None

    $Id: weight_fn.c 2761 2011-07-01 04:35:25Z gor171 $

*/

#include <stdio.h>
#include <math.h>
#include <sjwlib.h>

/** Compute's a tophat weighting function.
  *
  * @param x where to evaluate function.
  * @param scale length scale for function.
  * @return evaluate value.
  */
double
w_tophat(
double x,	/* Where to evaluate function */
double scale	/* Length scale for function */
)
{
    if( fabs(x) <= scale )
	return(1.0);
    else
	return(0.0);
}

/** Compute's a linear weighting function.
  *
  * @param x where to evaluate function.
  * @param scale length scale for function.
  * @return evaluate value.
  */
double
w_linear(
double x,	/* Where to evaluate function */
double scale	/* Length scale for function */
)
{
    if( fabs(x) <= scale )
	return(1.0 - x/scale);
    else
	return(0.0);
}

/** Compute's a parabolic weighting function.
  *
  * @param x where to evaluate function.
  * @param scale length scale for function.
  * @return evaluate value.
  */
double
w_parabolic(
double x,	/* Where to evaluate function */
double scale	/* Length scale for function */
)
{
    if( fabs(x) <= scale )
	return(1.0 - (x*x)/(scale*scale));
    else
	return(0.0);
}

/** Compute's a gaussian weighting function.
  *
  * @param x where to evaluate function.
  * @param scale length scale for function.
  * @return evaluate value.
  */
double
w_gaussian(
double x,	/* Where to evaluate function */
double scale	/* Length scale for function */
)
{
    return( exp( - x*x/(2*scale*scale)));
}
