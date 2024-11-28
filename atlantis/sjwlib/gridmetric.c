/*********************************************************************

    File:           gridmetric.c
    
    Created:        Thu Sep 10 12:45:09 EST 1992
    
    Author:         Stephen Walker
                    CSIRO Division of Oceanography
    
    Purpose:        Calculate grid metrics numerically
    
    Arguments:      see below

    Returns:        void
    
    Revisions:      none

    $Id: gridmetric.c 2761 2011-07-01 04:35:25Z gor171 $

*********************************************************************/

#include <math.h>
#include <stdio.h>
#include "sjwlib.h"

#define DEG2RAD(d) ((d)*M_PI/180.0)
#define RAD2DEG(r) ((r)*180.0/M_PI)

#define RADIUS 6370997.0
#define ECC 0.0
#define GEODESIC(x1, y1, x2, y2) (inv_sodanos(DEG2RAD(x1), DEG2RAD(y1),\
                                 DEG2RAD(x2), DEG2RAD(y2),\
                                 RADIUS, ECC))

/** Calculate grid metrics numerically.
  *
  * @param x grid x coordinates
  * @param y grid y coordinates
  * @param nce1 number of cells in e1 direction
  * @param nce2 number of cells in e2 direction
  * @param h1 storage for h1 values
  * @param h2 storage for h2 values
  */
void gridmetric(
double **x,		/* grid x coordinates */
double **y,		/* grid y coordinates */
long int nce1,		/* number of cells in e1 direction */
long int nce2,		/* number of cells in e2 direction */
double **h1,		/* storage for h1 values */
double **h2		/* storage for h2 values */
)
{
    long i,j;
    double dx;
    double dy;

    if( h1 != NULL ) {
        /* interior h1 values */
        for(j=0; j<nce2+1; j++)
	    for(i=1; i<nce1; i++) {
	        dx = x[j][i+1] - x[j][i-1];
	        dy = y[j][i+1] - y[j][i-1];
	        h1[j][i] = _hypot(dx,dy)/2.0;
	    }
        /* boundary h1 values */
        for(j=0; j<nce2+1; j++) {
	    /* h1 value when i == 0 */
	    dx = x[j][1] - x[j][0];
	    dy = y[j][1] - y[j][0];
            h1[j][0] = _hypot(dx,dy);
	    /* h1 value when i == nce1 */
	    dx = x[j][nce1] - x[j][nce1-1];
	    dy = y[j][nce1] - y[j][nce1-1];
            h1[j][nce1] = _hypot(dx,dy);
        }
    }

    if( h2 != NULL ) {
        /* interior h2 values */
        for(j=1; j<nce2; j++)
	    for(i=0; i<nce1+1; i++) {
	        dx = x[j+1][i] - x[j-1][i];
	        dy = y[j+1][i] - y[j-1][i];
	        h2[j][i] = _hypot(dx,dy)/2.0;
	    }
        /* boundary h2 values */
        for(i=0; i<nce1+1; i++) {
	    /* h2 value when j == 0 */
	    dx = x[1][i] - x[0][i];
	    dy = y[1][i] - y[0][i];
	    h2[0][i] = _hypot(dx,dy);
	    /* h2 value when j == nce2 */
	    dx = x[nce2][i] - x[nce2-1][i];
	    dy = y[nce2][i] - y[nce2-1][i];
	    h2[nce2][i] = _hypot(dx,dy);
        }
    }
}


/** Calculate grid metrics numerically for a geographic grid.
  *
  * @param x grid longitude coordinates (degrees)
  * @param y grid latitude coordinates (degrees)
  * @param nce1 number of cells in e1 direction
  * @param nce2 number of cells in e2 direction
  * @param h1 storage for h1 values
  * @param h2 storage for h2 values
  */
void geog_gridmetric(double **x, double **y, int nce1, int nce2,
                      double **h1, double **h2)
{
    int i,j;

    if( h1 != NULL ) {
        /* interior h1 values */
        for(j=0; j<nce2+1; j++)
            for(i=1; i<nce1; i++)
                h1[j][i] = GEODESIC(x[j][i+1], y[j][i+1], x[j][i-1], y[j][i-1])/2.0;

        /* boundary h1 values */
        for(j=0; j<nce2+1; j++) {
            /* h1 value when i == 0 */
            h1[j][0] = GEODESIC(x[j][1], y[j][1], x[j][0], y[j][0]);
            /* h1 value when i == nce1 */
            h1[j][nce1] = GEODESIC(x[j][nce1], y[j][nce1],
                                 x[j][nce1-1], y[j][nce1-1]);
        }
    }

    if( h2 != NULL ) {
        /* interior h2 values */
        for(j=1; j<nce2; j++)
            for(i=0; i<nce1+1; i++)
                h2[j][i] = GEODESIC(x[j+1][i], y[j+1][i], x[j-1][i], y[j-1][i])/2.0;

        /* boundary h2 values */
        for(i=0; i<nce1+1; i++) {
            /* h2 value when j == 0 */
            h2[0][i] = GEODESIC(x[1][i], y[1][i], x[0][i], y[0][i]);
            /* h2 value when j == nce2 */
            h2[nce2][i] = GEODESIC(x[nce2][i], y[nce2][i],
                              x[nce2-1][i], y[nce2-1][i]);
        }
    }
}
