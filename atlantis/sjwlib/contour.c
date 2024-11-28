/*

    File:           contour.c
    
    Created:        Tue Apr 20 10:43:46 EST 1993
    
    Author:         Stephen Walker
                    CSIRO Division of Oceanography
    
    Purpose:        Contour plotting routines compatible with my
		    postscript library
    
    Arguments:      See below

    Revisions:      15/12/1993
		    Changed scheme to collect line segments and then
		    join them to form long contour polylines. Added
		    labelling of contour values.

		    21/3/1994 SJW
		    Removed any calls to GMT stuff. Made compatible with
		    my postscript library

		    28/11/1997 SJW
		    Altered labelling algorithm to place labels at
		    1/4, 1/2 or 3/4 along contour lines

    $Id: contour.c 2761 2011-07-01 04:35:25Z gor171 $

*/

#include <ctype.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <sjwlib.h>

/* structure used by contouring routines */
struct pv {
    double x;
    double y;
    double val;
};

/* prototypes */
void conbox(struct pv *bl, struct pv *br, struct pv *tl, struct pv *tr,
	    double cval, double badval);
void contri(struct pv *p1, struct pv *p2, struct pv *p3, double cval);
int con_between(struct pv *p1, struct pv *p2, double val, struct pv *at);
void labelcontour(polyline pl, char *cval);
int lencmp( const void *pl1, const void *pl2);

/* Storage for contour line segments */
#ifndef THINK_C
#define    MAXSEGS    5000
#else
#define    MAXSEGS    500
#endif
polyline *seglist = NULL;
long    nseg;

/* polyline length comparison routine */
int
lencmp(const void *p1, const void *p2)
{
    polyline  *pl1 = (polyline *)p1;
    polyline  *pl2 = (polyline *)p2;

    if( pl1->np > pl2->np )
	return(-1);
    if( pl1->np < pl2->np )
	return(1);
    return(0);
}


/** Contour plotting routines compatible with Stephen Walker's
  * postscript library.
  *
  * @param val 2d array of values
  * @param x 2d array of x-coordinates
  * @param y 2d array of y-coordinates
  * @param nx number of points in x direction
  * @param ny number of points in y direction
  * @param cval pointer to the contour value
  * @param badval designated bad value.
  */
void
contour(
double **val,		/* 2d array of values */
double **x,		/* 2d array of x-coordinates */
double **y,		/* 2d array of y-coordinates */
long int nx, 		/* number of points in x direction */
long int ny, 		/* number of points in y direction */
char *cval, 		/* contour value */
double badval
)
{
    long i;
    long j;
    long n;
    struct pv bl;
    struct pv br;
    struct pv tl;
    struct pv tr;
    
    /* Allocate memory for segment list and set it all to zero */
    if( !seglist && (seglist=(polyline *)calloc(MAXSEGS,sizeof(polyline))) == NULL )
        quit("contour: No memory for segment list\n");

    /* Indicate contour value in file */
    printf("%% %s contour\n", cval);

    /* loop through data array, finding contour line segments */
    for(j=0; j<ny-1; j++)
	for(i=0; i<nx-1; i++) {
	    bl.x = x[j][i]; bl.y = y[j][i]; bl.val = val[j][i];
	    br.x = x[j][i+1]; br.y = y[j][i+1]; br.val = val[j][i+1];
	    tl.x = x[j+1][i]; tl.y = y[j+1][i]; tl.val = val[j+1][i];
	    tr.x = x[j+1][i+1]; tr.y = y[j+1][i+1]; tr.val = val[j+1][i+1];
            conbox(&bl,&br,&tl,&tr,atof(cval),badval);
	}
    
    /* join line segments to form entire contour lines */
    while( joinsegments(seglist,nseg,1e-10) > 0 ) /* loop */;
    while( joinsegments(seglist,nseg,1e-8) > 0 ) /* loop */;
    while( joinsegments(seglist,nseg,1e-6) > 0 ) /* loop */;

    /* despike contours */
    for(n=0; n<nseg; n++)
	if( seglist[n].np > 3 )
	    despike(&seglist[n],1e-6);
     
    /* sort contour lines into length order, longest first */
    qsort(seglist,(size_t)nseg,sizeof(polyline),lencmp);

    /* Loop over the contour lines, plotting, and labelling */
    for(n=0; n<nseg; n++) {
	if( seglist[n].np > 1 ) {
	    plotpolyline(stdout,seglist[n]);
	    /* label those lines which are at least one sixth as long
	     * as the longest (first) one
	     */
	    if( seglist[n].np > seglist[0].np/6 )
		labelcontour(seglist[n], cval);
	}
    }

    /* Clear the polyline list */
    for(n=0; n<nseg; n++)
	if( seglist[n].np > 0 )
	    clearpolyline(&seglist[n]);
    nseg = 0;
}

#define NODATA(x) (isnan(x) || fabs(x-badval) < 1e-6)

/* finds contour lines in a box with value cval */
void
conbox(struct pv *bl, struct pv *br, struct pv *tl,
       struct pv *tr, double cval, double badval)
{
    struct pv mid;
    double sum;
    int n;
    
    /* calculate midpoint of box */
    mid.x = (bl->x + br->x + tl->x + tr->x)/4;
    mid.y = (bl->y + br->y + tl->y + tr->y)/4;
    sum = 0.0;
    n = 0;
    if( !NODATA(bl->val) ) {
	sum += bl->val;
	n++;
    }
    if( !NODATA(br->val) ) {
	sum += br->val;
	n++;
    }
    if( !NODATA(tl->val) ) {
	sum += tl->val;
	n++;
    }
    if( !NODATA(tr->val) ) {
	sum += tr->val;
	n++;
    }

    if( n > 2 ) {
	mid.val = sum/n;
	/* contour triangles making up box */
	if( !NODATA(tl->val) && !NODATA(bl->val) )
	    contri(&mid, tl, bl, cval);
	if( !NODATA(bl->val) && !NODATA(br->val) )
	    contri(&mid, bl, br, cval);
	if( !NODATA(br->val) && !NODATA(tr->val) )
	    contri(&mid, br, tr, cval);
	if( !NODATA(tr->val) && !NODATA(tl->val) )
	    contri(&mid, tr, tl, cval);
    }
}

/* finds contour lines in a triangle with value cval */
void
contri(struct pv *p1, struct pv *p2, struct pv *p3, double cval)
{
    struct pv where;
    int n;
    dpoint p[3];
    
    /* Check that there is room for more segments */
    if( nseg >= MAXSEGS ) {
	fprintf(stderr,"Contour segments may be missing - list full\n");
	return;
    }

    /* Special case - all points equal - no contour */
    if( p1->val == p2->val && p2->val == p3->val )
	return;
    
    /* Special case - contour exactly along a side -
     * the best approach here is to jitter the contour
     * value slightly - other techniques lead to trouble
     */
    if( (cval == p1->val && cval == p2->val) ||
        (cval == p2->val && cval == p3->val) ||
        (cval == p3->val && cval == p1->val)  ) {
	if( cval == 0.0 ) cval = 1e-12;
	else cval *= 1.000000000001;
    }

    /* Normal case - check for contour intersections on each side */
    n = 0;
    if( con_between(p1,p2,cval,&where) ) {
        p[n].x = where.x;
        p[n].y = where.y;
        n++;
    }
    if( con_between(p2,p3,cval,&where) ) {
        p[n].x = where.x;
        p[n].y = where.y;
        n++;
    }
    if( con_between(p3,p1,cval,&where) ) {
        p[n].x = where.x;
        p[n].y = where.y;
        n++;
    }

    /* store segment if found and not zero length */
    if( n == 2 && !(p[0].x == p[1].x && p[0].y == p[1].y) ) {
        addtoend(&seglist[nseg],p[0]);
        addtoend(&seglist[nseg],p[1]);
        nseg++;
    }
    
    /* deal with special case where contour passes through
     * a corner and also the opposite side
     */
    else if( n == 3 ) {
        if( p[0].x == p[1].x && p[0].y == p[1].y ) {
	    addtoend(&seglist[nseg],p[1]);
	    addtoend(&seglist[nseg],p[2]);
	    nseg++;
	}
        else if( p[0].x == p[2].x && p[0].y == p[2].y ) {
	    addtoend(&seglist[nseg],p[0]);
	    addtoend(&seglist[nseg],p[1]);
	    nseg++;
	}
        else if( p[1].x == p[2].x && p[1].y == p[2].y ) {
	    addtoend(&seglist[nseg],p[0]);
	    addtoend(&seglist[nseg],p[1]);
	    nseg++;
	}
    }
}	

int
con_between(struct pv *p1, struct pv *p2, double val, struct pv *at)
{
    double frac;
    
    if( (p1->val == p2->val)             ||
        (val < p1->val && val < p2->val) ||
        (val > p1->val && val > p2->val) )
	return(0);
    
    frac = (val - p1->val)/(p2->val - p1->val);
    at->x = p1->x + (p2->x - p1->x)*frac;
    at->y = p1->y + (p2->y - p1->y)*frac;
    at->val = val;
    return(1);
}

void
labelcontour(polyline pl, char *cval)
{
    int i;
    int mid;
    linepoint *ptr;
    double x1, y1, x2, y2;
    static int pos = 1;

    /* reject really small lines */
    if( pl.np <= 4 )
	return;
    
    /* find the position of the segment to be labelled */
    mid = pos*(int)pl.np/4 - 1;
    for(ptr=pl.start, i=0; ptr!=NULL && i++<mid; ptr = ptr->next) /* loop */;

    x1 = ptr->p.x;
    y1 = ptr->p.y;

    /* FIX - fudge to avoid some sort of rare rounding error bug
     * SJW 12/3/1996
     */
    if( !ptr->next )
	return;
    x2 = ptr->next->p.x;
    y2 = ptr->next->p.y; 
    printf("(%s)\n",cval); 
    printf("%e %e dat2page %e %e dat2page\n",x2,y2,x1,y1);
    printf("exch 3 1 roll sub 3 1 roll sub atan\n"); 
    printf("dup 180 gt { 360 sub} if\n");
    printf("dup 90 gt { 180 sub} if\n");
    printf("dup -90 lt { 180 add} if\n");
    printf("%e %e dat2page txt_CC\n",(x1+x2)/2, (y1+y2)/2);

    pos++;
    if( pos > 3 ) pos = 1;
}
