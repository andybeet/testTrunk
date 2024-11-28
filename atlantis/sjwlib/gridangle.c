/*

    File:           gridangle.c
    
    Created:        Thu Sep 10 13:50:09 EST 1992
    
    Author:         Stephen Walker
                    CSIRO Division of Oceanography
    
    Purpose:        Calculate grid angle numerically
    
    Arguments:      see below

    Returns:        void
    
    Revisions:      none

    $Id: gridangle.c 2761 2011-07-01 04:35:25Z gor171 $

*/

#include <math.h>
#include <stdio.h>
#include "sjwlib.h"

#define DEG2RAD(d) ((d)*M_PI/180.0)
#define RAD2DEG(r) ((r)*180.0/M_PI)

#define RADIUS 6370997.0
#define ECC 0.0
#define GANGLE(x1, y1, x2, y2) (M_PI/2 - inv_sodanos_angles(\
				DEG2RAD(x1), DEG2RAD(y1),\
                                DEG2RAD(x2), DEG2RAD(y2),\
                                RADIUS, ECC, NULL))

/** Calculate grid angle numerically.
  *
  * @param x grid x coordinates
  * @param y grid y coordinates
  * @param nce1 number of cells in e1 direction
  * @param nce2 number of cells in e2 direction
  * @param a1 storage for angles between x and e1
  * @param a2 storage for angles between x and e2
  */
void gridangle(
double **x,		/* grid x coordinates */
double **y,		/* grid y coordinates */
long int nce1,		/* number of cells in e1 direction */
long int nce2,		/* number of cells in e2 direction */
double **a1,		/* storage for angles between x and e1 */
double **a2		/* storage for angles between x and e2 */
)
{
    long i,j;
    double dx;
    double dy;

    if( a1 != NULL ) {
        /* interior a1 values */
        for(j=0; j<nce2+1; j++)
	    for(i=1; i<nce1; i++) {
	        dx = x[j][i+1] - x[j][i-1];
	        dy = y[j][i+1] - y[j][i-1];
                a1[j][i] = atan2(dy,dx);
	    }
        /* boundary a1 values */
        for(j=0; j<nce2+1; j++) {
	    /* a1 value when i == 0 */
	    dx = x[j][1] - x[j][0];
	    dy = y[j][1] - y[j][0];
            a1[j][0] = atan2(dy,dx);
	    /* a1 value when i == nce1 */
	    dx = x[j][nce1] - x[j][nce1-1];
	    dy = y[j][nce1] - y[j][nce1-1];
            a1[j][nce1] = atan2(dy,dx);
        }
    }

    if( a2 != NULL ) {
        /* interior a2 values */
        for(j=1; a2 != NULL && j<nce2; j++)
	    for(i=0; i<nce1+1; i++) {
	        dx = x[j+1][i] - x[j-1][i];
	        dy = y[j+1][i] - y[j-1][i];
                a2[j][i] = atan2(dy,dx);
	    }
        /* boundary a2 values */
        for(i=0; i<nce1+1; i++) {
	    /* a2 value when j == 0 */
	    dx = x[1][i] - x[0][i];
	    dy = y[1][i] - y[0][i];
            a2[0][i] = atan2(dy,dx);
	    /* a2 value when j == nce2 */
	    dx = x[nce2][i] - x[nce2-1][i];
	    dy = y[nce2][i] - y[nce2-1][i];
            a2[nce2][i] = atan2(dy,dx);
        }
    }
}

/** Calculate grid angles numerically for a geographic grid.
  *
  * @param x grid longitude coordinates (degrees)
  * @param y grid latitude coordinates (degrees)
  * @param nce1 number of cells in e1 direction
  * @param nce2 number of cells in e2 direction
  * @param a1 storage for angles between x and e1
  * @param a2 storage for angles between x and e2
  */
void geog_gridangle(double **x, double **y, int nce1, int nce2,
                      double **a1, double **a2)
{
    int i,j;
 
    if( a1 != NULL ) {
        /* interior a1 values */
        for(j=0; j<nce2+1; j++)
            for(i=0; i<nce1; i++)
                a1[j][i] = GANGLE(x[j][i], y[j][i], x[j][i+1], y[j][i+1]);

        /* boundary a1 values */
        for(j=0; j<nce2+1; j++) {
            /* a1 value when i == nce1 */
            a1[j][nce1] = GANGLE(x[j][nce1], y[j][nce1],
                                 x[j][nce1-1], y[j][nce1-1]) - M_PI;
        }
    }

    if( a2 != NULL ) {
        /* interior a2 values */
        for(j=0; j<nce2; j++)
            for(i=0; i<nce1+1; i++)
                a2[j][i] = GANGLE(x[j][i], y[j][i], x[j+1][i], y[j+1][i]);

        /* boundary a2 values */
        for(i=0; i<nce1+1; i++) {
            /* a2 value when j == nce2 */
            a2[nce2][i] = GANGLE(x[nce2][i], y[nce2][i],
                              x[nce2-1][i], y[nce2-1][i]) - M_PI;
        }
    }
}
