/*

    File:           ellipt_coord.c
    
    Created:        Thu Aug 13 19:12:21 EST 1992
    
    Author:         Stephen Walker
                    CSIRO Division of Oceanography
    
    Purpose:        Routines to calculate grid coordinates for
		    m3d model grids with elliptic geometries
    
    Arguments:      lots - see each routine

    Returns:        void
    
    Revisions:      This routine is a standalone version
		    of that previously used in gengrid, which
		    used global variables to get the grid parameters

    $Id: ellipt_coord.c 2761 2011-07-01 04:35:25Z gor171 $

*/

#include <math.h>
#include <stdio.h>
#include "sjwlib.h"

#define	DEG2RAD		(M_PI/180.0)

/** Calculates grid coordinates for model grids with elliptic
  * geometries.
  *
  * @param x where to store grid x values
  * @param y where to store grid y values
  * @param h1 where to store h1 metric values
  * @param h2 where to store h2 metric values
  * @param a1 where to store a1 angle values
  * @param a2 where to store a2 angle values
  * @param nce1 number of cells in e1 direction
  * @param nce2 number of cells in e2 direction
  * @param x00 x origin offset
  * @param y00 y origin offset
  * @param rotn angle (degrees) between East and sigma=1 axis
  * @param ella distance from origin to either focus point
  * @param taumax max tau value ( 1 >= taumax > taumin >= -1 )
  * @param taumin min tau value ( 1 >= taumax > taumin >= -1 )
  * @param nsimm number of symmetric cells in e2 reflected about sigma = 1
  */
void ellipt_coord(
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
double rotn,	/* angle (degrees) between East and sigma=1 axis */
double ella,	/* distance from origin to either focus point */
double taumax,	/* max tau value ( 1 >= taumax > taumin >= -1 ) */
double taumin,	/* max tau value ( 1 >= taumax > taumin >= -1 ) */
long int nsimm	/* number of symmetric cells in e2 reflected about sigma = 1 */
)
{
    int i,j;
    double dtau;
    double sinth;
    double costh;
    double th;
    double sigma;
    double tau;
    double xval;
    double yval;
    double v;

    th = rotn*DEG2RAD;
    sinth = sin(th);
    costh = cos(th);

    /* tau spacing */
    dtau = (taumax - taumin)/(double)nce1;
    sigma = 1.0;
    for(j=0; j<nce2+1; j++) {
	for(i=0; i<nce1+1; i++) {
	    tau = taumin + i*dtau;
	    /* (x,y) before rotation and offset */
	    x[j][i] = ella*sigma*tau;
	    y[j][i] = ella*sqrt((sigma*sigma-1)*(1-tau*tau));
	}
	/* calculate new sigma to make cells roughly square on y axis */
	v = sigma*dtau + sqrt(sigma*sigma-1);
	sigma = sqrt(v*v+1); 
    }

    /* shift and reflect as necessary for symmetric cells */
    for(j=(int)nce2; j>=(int)nsimm; j--)
	for(i=0; i<nce1+1; i++) {
	    x[j][i] = x[j-nsimm][i];
	    y[j][i] = y[j-nsimm][i];
	}
    for(j=0; j<nsimm; j++)
	for(i=0; i<nce1+1; i++) {
	    x[j][i] = x[2*nsimm-j][i];
	    y[j][i] = -y[2*nsimm-j][i];
	}

    /* do rotation and offset */
    for(j=0; j<nce2+1; j++)
	for(i=0; i<nce1+1; i++) {
	    xval = x[j][i];
	    yval = y[j][i];
            x[j][i] = x00 + xval*costh - yval*sinth;
            y[j][i] = y00 + xval*sinth + yval*costh;
        }

    /* calculate h1 and h2 values numerically */
    gridmetric(x,y,nce1,nce2,h1,h2);
    /* calculate a1 and a2 angle values numerically */
    gridangle(x,y,nce1,nce2,a1,a2);
}
