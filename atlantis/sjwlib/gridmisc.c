/*

    File:           gridmisc.c
    
    Created:        Mon Nov 13 16:16:58 EST 1995
    
    Author:         Stephen Walker
                    CSIRO Division of Oceanography
    
    Purpose:        Miscellaneous grid routines
    
    Arguments:      

    Returns:        void
    
    Revisions:      none

    $Id: gridmisc.c 2761 2011-07-01 04:35:25Z gor171 $

*/

#include <stdio.h>
#include <math.h>
#include "sjwlib.h"

void
CentreToCorner(int nce1, int nce2, double **cx, double **cy, double **gx, double **gy)
{
    int i, j;

    /* Generate interior grid coordinates by interpolation */
    /* Interior points */
    for(j=1; j<nce2; j++) {
        for(i=1; i<nce1; i++) {
            gx[j][i] = (cx[j][i]+cx[j][i-1]+cx[j-1][i]+cx[j-1][i-1])/4;
            gy[j][i] = (cy[j][i]+cy[j][i-1]+cy[j-1][i]+cy[j-1][i-1])/4;
        }
    }
    /* Left and right sides of grid */
    for(j=1; j<nce2; j++) {
        gx[j][0] = cx[j][0]+cx[j-1][0] - gx[j][1];
        gy[j][0] = cy[j][0]+cy[j-1][0] - gy[j][1];
        gx[j][nce1] = cx[j][nce1-1]+cx[j-1][nce1-1] - gx[j][nce1-1];
        gy[j][nce1] = cy[j][nce1-1]+cy[j-1][nce1-1] - gy[j][nce1-1];
    }
    
    /* Bottom and top sides of grid */
    for(i=1; i<nce1; i++) {
        gx[0][i] = cx[0][i]+cx[0][i-1] - gx[1][i];
        gy[0][i] = cy[0][i]+cy[0][i-1] - gy[1][i];
        gx[nce2][i] = cx[nce2-1][i]+cx[nce2-1][i-1] - gx[nce2-1][i];
        gy[nce2][i] = cy[nce2-1][i]+cy[nce2-1][i-1] - gy[nce2-1][i];
    }
    
    /* Corners */
    gx[0][0] = gx[1][0] + gx[0][1] - gx[1][1];
    gy[0][0] = gy[1][0] + gy[0][1] - gy[1][1];
    gx[nce2][0] = gx[nce2][1] + gx[nce2-1][0] - gx[nce2-1][1];
    gy[nce2][0] = gy[nce2][1] + gy[nce2-1][0] - gy[nce2-1][1];
    gx[0][nce1] = gx[0][nce1-1] + gx[1][nce1] - gx[1][nce1-1];
    gy[0][nce1] = gy[0][nce1-1] + gy[1][nce1] - gy[1][nce1-1];
    gx[nce2][nce1] = gx[nce2][nce1-1] + gx[nce2-1][nce1] - gx[nce2-1][nce1-1];
    gy[nce2][nce1] = gy[nce2][nce1-1] + gy[nce2-1][nce1] - gy[nce2-1][nce1-1];
}
