/*********************************************************************

    File:           xytoij.h

    Created:        Thr Jan 23 14:00:00 EST 1997

    Author:         Daniel Delbourgo/Stephen Walker/Jason Waring
                    CSIRO Division of Oceanography

    Purpose:        calculates the indices (i,j) of a
		    topologically rectangular grid
		    cell containing the point
		    (x,y)

    Arguments:      x - world x coordinate (double)
                    y - world y coordinate (double)
		    i - pointer to storage for i
		    j - pointer to storage for j

    Returns:        1 if point is inside any model cell
		    0 if point is outside all model cells

    Revisions:      110898 JRW
		    Added IJtoXY conversion and
		    fractional XYtoIJ and IJtoXY conversion.

    $Id: xytoij.h 2761 2011-07-01 04:35:25Z gor171 $

*********************************************************************/

#ifndef _SJW_XYTOIJ_H
#define _SJW_XYTOIJ_H

#ifndef _SJWLIB_H
#include "sjwlib.h"
#endif

typedef struct leaf
{
	polyline *boundary;
	long mini;
	long maxi;
	long minj;
	long maxj;
	struct leaf	*half1;
	struct leaf	*half2;
} leaf;

typedef struct
{
	long leaves;
	polyline *outline;
	leaf	*trunk;

/**** PRIVATE - used for housekeeping and efficiency reasons.
 ****/
	int lasti;
	int lastj;
	polyline *lastcell;

	double **gridx;	/* Array of X coords (of size (nce1+1)*(nce2+1)). */
	double **gridy;	/* Array of Y coords (of size (nce1+1)*(nce2+1)). */
	int    nce1;	/* Number of cells in e1 direction. */
	int    nce2;	/* Number of cells in e2 direction. */

} tree;

/* Proto-types */
tree *xytoij_init(double **gx, double **gy, int nce1, int nce2);
int xytoij(tree *partition, double x, double y, int *ival, int *jval);
int ijtoxy(tree *partition, int ival, int jval, double *x, double *y);
int xytofij(tree *partition, double x, double y, double *ival, double *jval);
int fijtoxy(tree *partition, double ival, double jval, double *x, double *y);

#endif
