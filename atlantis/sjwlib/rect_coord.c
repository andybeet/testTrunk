/*********************************************************************

    File:           rect_coord.c
    
    Created:        Thu Aug 13 19:12:21 EST 1992
    
    Author:         Stephen Walker
                    CSIRO Division of Oceanography
    
    Purpose:        Routines to calculate grid coordinates for
		    m3d model grids with rectangular geometries
    
    Arguments:      lots - see below

    Returns:        void
    
    Revisions:      This routines is a standalone version
		    of that previously used in gengrid, which
		    used global variables to get the grid parameters

    $Id: rect_coord.c 2761 2011-07-01 04:35:25Z gor171 $

*********************************************************************/

#include <math.h>
#include <stdio.h>
#include "sjwlib.h"

#define DEG2RAD		(M_PI/180.0)

/** Calculate coordinates for rectangular grid.
  *
  * @param x where to store grid x values
  * @param y where to store grid y values
  * @param h1where to store h1 metric values
  * @param h2 where to store h2 metric values
  * @param a1 where to store a1 angle values
  * @param a2 where to store a2 angle values
  * @param nce1 number of cells in e1 direction
  * @param nce2 number of cells in e2 direction
  * @param x00 x origin offset
  * @param y00 y origin offset
  * @param rotn angle (degrees) between East and e1 axis
  * @param xinc cell size in x direction
  * @param yinc cell size in y direction
  */
void
rect_coord(
double **x,	/* where to store grid x values */
double **y,	/* where to store grid y values */
double **h1,	/* where to store h1 metric values */
double **h2,	/* where to store h2 metric values */
double **a1,	/* where to store a1 angle values */
double **a2,	/* where to store a2 angle values */
long int nce1,	/* number of cells in e1 direction */
long int nce2,	/* number of cells in e2 direction */
double x00,	/* x origin offset */
double y00,	/* y origin offset */
double rotn,	/* angle (degrees) between East and e1 axis */
double xinc,	/* cell size in x direction */
double yinc	/* cell size in y direction */
)
{
    long i,j;
    double xval,yval;
    double sinth;
    double costh;

    sinth = sin(rotn*DEG2RAD);
    costh = cos(rotn*DEG2RAD);
    for(j=0; j<nce2+1; j++) {
	yval = (double)j*yinc;
        for(i=0; i<nce1+1; i++) {
	    xval = (double)i*xinc;
            x[j][i] = x00 + xval*costh - yval*sinth;
            y[j][i] = y00 + xval*sinth + yval*costh;
	    /***** Analytic calculation deleted
	    h1[j][i] = fabs(xinc);
	    h2[j][i] = fabs(yinc);
	    *******/
        }
    }

    /* Calculate h1 and h2 numerically */
    gridmetric(x,y,nce1,nce2,h1,h2);
    /* calculate a1, a2 numerically */
    gridangle(x,y,nce1,nce2,a1,a2);
}
