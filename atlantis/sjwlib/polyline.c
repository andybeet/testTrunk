/*

    File:           polyline.c

    Created:        Mon Aug 24 15:33:54 EST 1992

    Author:         Stephen Walker
                    CSIRO Division of Oceanography

    Purpose:        Library routines for points, rectangles
		    and polylines

    Arguments:      See each routine

    Returns:        ditto

    Revisions:      none

    $Id: polyline.c 3369 2012-08-29 06:16:46Z gor171 $

*/

#include <ctype.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "sjwlib.h"
#include <atlantisboxmodel.h>
#include <atUtilLib.h>


/** Create a new instance of a dpoint.
  *
  * @param x X coordinate.
  * @param y Y coordinate.
  * @return pointer to a dpoint.
  */
dpoint *createdpoint(double x, double y)
{
	dpoint *p;

	if( (p=(dpoint *)malloc(sizeof(dpoint))) == NULL )
	    quit("createdpoint: memory allocation failure\n");
	p->x = x;
	p->y = y;
	return(p);
}

/** Create a line point from a dpoint.
  *
  * @param p dpoint.
  * @return pointer to a linepoint.
  */
linepoint
*createlinepoint(dpoint p)
{
	linepoint *ptr;

	if( (ptr=(linepoint *)malloc(sizeof(linepoint))) == NULL )
	    quit("createlinepoint: memory allocation failure\n");
	ptr->p = p;
	ptr->next = NULL;
	return(ptr);
}


/** Create an empy polyline.
  *
  * @return pointer to a polyline structure.
  */
polyline *
createpolyline(void)
{
	polyline *pl;

	if( (pl=(polyline *)malloc(sizeof(polyline))) == NULL )
	    quit("createpolyline: memory allocation failure\n");
	pl->np = 0;
	pl->start = NULL;
	pl->end = NULL;
	pl->bound = def_rect(0.0,0.0,0.0,0.0);
	return(pl);
}

/**
 * Free the memory associated with a dynamically-allocated polyline
 * (such as one created via createpolyline).
 */
void
freepolyline(polyline *pl)
{
    linepoint *lp = pl->start, *hold;
    while ((hold = lp)) {
        lp = lp->next;
        free(hold);
    }
    free(pl);
}

/** Create a rectangle.
  *
  * @param x1 X coordinate of bottom left corner.
  * @param y1 Y coordinate of bottom left corner.
  * @param x2 X coordinate of top right corner.
  * @param y2 Y coordinate of top right corner.
  * @return pointer to a rect structure.
  */
rect
*createrect(double x1, double y1, double x2, double y2)
{
	rect *r;
	if( (r=(rect *)malloc(sizeof(rect))) == NULL )
	    quit("createrect: memory allocation failure\n");
	r->botleft.x = x1;
	r->botleft.y = y1;
	r->topright.x = x2;
	r->topright.y = y2;
	return(r);
}


/** Creates a dpoint instance.
  *
  * @param x X coordinate.
  * @param y Y coordinate.
  * @return dpoint structure.
  */
dpoint
def_dpoint(double x, double y)
{
	dpoint p;

	p.x = x;
	p.y = y;
	return(p);
}


/** Creates a rect instance.
  *
  * @param x1 X coordinate of bottom left corner.
  * @param y1 Y coordinate of bottom left corner.
  * @param x2 X coordinate of top right corner.
  * @param y2 Y coordinate of top right corner.
  * @return rect structure.
  */
rect
def_rect(double x1, double y1, double x2, double y2)
{
	rect r;
	r.botleft.x = x1;
	r.botleft.y = y1;
	r.topright.x = x2;
	r.topright.y = y2;
	return(r);
}


/** Tests whether two pointer are the same. A threshold distance
  * is used to define the tolerance within which two point
  * maybe considered the same.
  *
  * @param p1 first point
  * @param p2 second point
  * @param dist threshold distance.
  * @return non-zero if successful.
  */
int
samepoint(dpoint p1, dpoint p2, double dist)
{
	return( (fabs(p1.x-p2.x)<=dist) && (fabs(p1.y-p2.y)<=dist) );
}


/** Tests whether a point is within a rectangle.
  *
  * @param p point to test
  * @param r rectangle.
  * @return non-zero if successful.
  */
int
pointinrect(dpoint p, rect r)
{
	if( p.x >= r.botleft.x &&
	    p.x <= r.topright.x &&
	    p.y >= r.botleft.y &&
	    p.y <= r.topright.y )
	    return(1);
	return(0);
}

/** Tests whether a point is inside a polygon.
  * Here the polyline is assumed to be closed: an extra
  * line segment is assumed from the end point back to the
  * start point
  *
  * @param p point to test
  * @param pl closed polyline to test against.
  * @return non-zero if successful.
  */
int
pointinpoly(dpoint p, polyline pl)
{
	linepoint *lp1;
	linepoint *lp2;
	dpoint p1;
	dpoint p2;
	long cnum;
	double x;

	if( pl.np <= 1 )
	    return(0);
	if( !pointinrect(p,pl.bound) )
	    return(0);

	/* loop over each line segment, counting crossings
	 * by polygon of the horizontal line from point p to
	 * +infinity
	 * Full crossings increment cnum by 2.
	 * Half crossings (where a segment starts or ends on the
	 * horizontal line) increment or decrement cnum by 1,
	 * depending on the direction of the approach to the line.
	 * This allows crossings to be counted correctly even when
	 * part of the polygon just touches the horizontal line
	 * but doesn't really cross it.
	 */
	cnum = 0;
	for(lp1=pl.start; lp1!=NULL; lp1 = lp1->next) {
	    lp2 = lp1->next;
	    if( lp2 == NULL )
		lp2 = pl.start;
	    /* shift origin to p */
	    p1.x = lp1->p.x - p.x;
	    p1.y = lp1->p.y - p.y;
	    p2.x = lp2->p.x - p.x;
	    p2.y = lp2->p.y - p.y;
	    if( p1.y == 0.0 && p2.y == 0.0 ) {
		/* polygon segment is part of x axis - ignore it unless the
		 * point is within the segment
		 */
		if( (p1.x<=0.0 && p2.x>=0.0) ||
		    (p1.x>=0.0 && p2.x<=0.0) )
		    return(1);      /* point is on polygon boundary */
	    }
	    else if( p1.y == 0.0 ) {
		/* polygon segment starts on x axis */
		if( p1.x == 0.0 )
		    return(1);      /* point is on polygon boundary */
		if( p1.x > 0.0 )
		    /* record a 'half' crossing with appropriate sign */
		    cnum += (p2.y>0.0) ? 1 : -1;
	    }
	    else if( p2.y == 0.0 ) {
		/* polygon segment ends on x axis */
		if( p2.x == 0.0 )
		    return(1);      /* point is on polygon boundary */
		if( p2.x > 0.0 )
		    /* record a 'half' crossing with appropriate sign */
		    cnum += (p1.y<0.0) ? 1 : -1;
	    }
	    else if( (p1.y<0.0 && p2.y>0.0)||(p1.y>0.0 && p2.y<0.0) ) {
		/* segment crosses x axis, determine if intercept >= 0 */
		if( p1.x >0 && p2.x > 0)
		    cnum += 2;
		else if( (p1.x<=0 && p2.x>=0) ||
			 (p1.x>=0.0 && p2.x<=0.0) ) {
		    /* have to do some maths here */
                    x = p1.x - (p2.x-p1.x)*p1.y/(p2.y-p1.y);
		    if( x == 0.0 )
			return(1);    /* point is on polygon boundary */
		    if( x > 0.0 )
			cnum += 2;
		}
	    }
	}

	/* cnum should now be two times the actual number of
	 * crossings. If it is odd, something has gone wrong
	 */
	if( cnum%2 )
	    warn("pointinpoly: Error in counting crossings\n");

	/* divide by 2 to get the actual number of crossings -
	 * if odd, the point is inside
	 */
        return(int)( (cnum/2)%2 );
}


/** Tests whether a point is one on a polyline. This routine
  * test the control points only.
  *
  * @param p point to test
  * @param pl polyline to test against.
  * @return non-zero if successful.
  */
int
pointonpolyline(dpoint p, polyline pl)
{
	linepoint *ptr;

	if( pl.np <= 0 )
	    return(0);
	for( ptr=pl.start; ptr!=NULL; ptr = ptr->next ) {
	    if( ptr->p.x == p.x &&
		ptr->p.y == p.y )
		return(1);
	}
	return(0);
}


/** Tests whether two rectangles intersect.
  *
  * @param r1 first rectangle.
  * @param r2 second rectangle.
  * @return non-zero if successful.
  */
int
sectrect(rect r1, rect r2)
{
	double x1;
	double y1;
	double x2;
	double y2;

	x1 = max(r1.botleft.x, r2.botleft.x);
	y1 = max(r1.botleft.y, r2.botleft.y);
	x2 = min(r1.topright.x, r2.topright.x);
	y2 = min(r1.topright.y, r2.topright.y);
	if( x1 < x2 && y1 < y2 )
	    return(1);
	return(0);
}


/** Clear a polyline structure of all memory allocated.
  *
  * @param pl pointer to the polyline.
  */
void
clearpolyline(polyline *pl)
{
    linepoint *ptr;
    linepoint *next;

    if( pl == NULL || pl->np <= 0 )
    	quit("clearpolyline: NULL or empty line\n");

    pl->np = 0;
    pl->bound = def_rect(0.0,0.0,0.0,0.0);
    pl->end = NULL;
    for(ptr=pl->start; ptr!=NULL; ) {
	next = ptr->next;
	free(ptr);
	ptr = next;
    }
    pl->start = NULL;
}


/** Reads a polyline from a file stream and store
  * within a polyline.
  *
  * @param fp file stream. open and readable.
  * @param pl pointer to returned polyline.
  * @return the number of points in the polyline or 0.
  */
long
readpolyline(FILE *fp, polyline *pl)
{
    dpoint p;
    char line[MAXLINELEN];
    int more;

    /* Skip comments and blank lines */
    if( nextline(line,MAXLINELEN,fp) == 0 )
	return(0);

    pl->np = 0;
    pl->bound = def_rect(0.0,0.0,0.0,0.0);
    pl->start = NULL;
    pl->end = NULL;
    more = 1;
    while( more && sscanf(line,"%lf %lf",&p.x,&p.y) == 2 ) {
	addtoend(pl,p);
	if( fgets(line,MAXLINELEN,fp) == NULL )
	    more = 0;
    }
    return(pl->np);
}


/** Write a polyline definition to a file stream.
  *
  * @param fp file stream. open and writable.
  * @param pl pointer to polyline.
  */
void
printpolyline(FILE *fp, polyline *pl)
{
	linepoint *ptr;

	if( pl == NULL || pl->np <= 0 )
	    return;

	fprintf(fp,"# %ld\n",pl->np);
	for(ptr=pl->start; ptr != NULL; ptr = ptr->next )
	    fprintf(fp,"%.4f %.4f\n",ptr->p.x, ptr->p.y);
}

/** Resample a polyline. Only including points more than
  * threshold distance apart.
  * The original start and end points are preserved.
  *
  * @param pl pointer to polyline structure.
  * @param res threshold distance.
  * @return a new resampled polyline.
  */
polyline *
resamplepolyline(polyline *pl, double res)
{
	polyline *new;
	linepoint *ptr;

	if( pl == NULL || pl->np <= 0 )
	    quit("resamplepolyline: NULL polyline\n");

	/* make new subsampled polyline */
	new = createpolyline();
	/* add original start point */
	addtoend(new,pl->start->p);
	for(ptr=pl->start->next; ptr!=NULL; ptr = ptr->next )
	    /* add each point to end of new line if it is more
	     * than res away from last point added
	     */
	    if( !samepoint(ptr->p, new->end->p, res) )
	        addtoend(new,ptr->p);
	/* add original end point if we haven't already */
	if( !samepoint(pl->end->p, new->end->p, 0.0) )
	    addtoend(new, pl->end->p);
	return(new);
}

/** Smooth a polyline. Fits a spline to each of x and y as
  * functions of distance along the line.
  *
  * @param pl pointer to polyline.
  * @param ns number of points in smoothed polyline.
  * @return new smoothed polyline.
  */
polyline *
smoothpolyline(polyline *pl, long int ns)
{
    polyline *new;
    linepoint *ptr;
    double *x1;	/* original x values */
    double *y1;	/* original y values */
    double *d;	/* distance from start values */
    double *ydd;	/* spline computed second derivatives */
    double *x2;	/* smoothed x values */
    double *y2;	/* smoothed y values */
    double x;
    double y;
    double prevx;
    double prevy;
    double dist;
    long i;
    long np;

    /* check for sensible polyline */
    if( pl == NULL )
        quit("smoothpolyline: NULL polyline\n");
    if( pl->np < 3 )
        quit("smoothpolyline: polyline must have 3 or more points\n");

    /* allocate memory for x1, y1, dist, x2, y2 */
    x1 = alloc1d(pl->np);
    y1 = alloc1d(pl->np);
    d = alloc1d(pl->np);
    ydd = alloc1d(pl->np);
    x2 = alloc1d(ns);
    y2 = alloc1d(ns);

    /* load x1, y1, d */
    i = 0;
    prevx = pl->start->p.x;
    prevy = pl->start->p.y;
    for(ptr=pl->start; ptr != NULL; ptr = ptr->next ) {
	x = ptr->p.x;
	y = ptr->p.y;
	x1[i] = x;
	y1[i] = y;
	if( i > 0 )
            d[i] = d[i-1] + sqrt((x-prevx)*(x-prevx)+(y-prevy)*(y-prevy));
	else
            d[i] = 0.0;
	prevx = x;
	prevy = y;
	i++;
    }
    np = i;

    /* check if closed */
    if( samepoint(pl->start->p, pl->end->p, 0.001) ) {
	/* calculate periodic boundary condition */
	/* NOT IMPLEMENTED YET */
	;
    }

    /* fit spline to x */
    dspline(d,x1,np,0,1e31,1e31,ydd);
    /* loop to evaluate interpolated x values */
    for(i=0; i<ns; i++) {
	/* calculate distance, evaluate spline */
	if( i< ns-1 ) dist = (double)i*d[pl->np -1] / (double)(ns-1);
	else dist = d[pl->np -1];
	dsplint(d,x1,ydd,np,dist,&x2[i]);
    }

    /* fit spline to y */
    dspline(d,y1,np,0,1e31,1e31,ydd);
    /* loop to evaluate interpolated y values */
    for(i=0; i<ns; i++) {
	/* calculate distance, evaluate spline */
	if( i< ns-1) dist = (double)i*d[pl->np -1] / (double)(ns-1);
	else dist = d[pl->np -1];
	dsplint(d,y1,ydd,np,dist,&y2[i]);
    }

    /* create new smoothed polyline */
    new = createpolyline();
    for(i=0; i<ns; i++)
	addtoend(new,def_dpoint(x2[i],y2[i]));

    /* free memory for x1, y1, dist, x2, y2 */
    free1d(x1);
    free1d(y1);
    free1d(d);
    free1d(x2);
    free1d(y2);

    return(new);
}



/** Add a point to the start of a polyline.
  *
  * @param pl pointer to polyline.
  * @param point to add.
  */
void
addtostart(polyline *pl, dpoint p)
{
        linepoint *ptr;

	/* create new linepoint */
	ptr = createlinepoint(p);

	if( pl->np <= 0 ) {
	    /* no points already in polyline, so add this one */
	    pl->start = ptr;
	    pl->end = ptr;
	    pl->np = 1;
	    pl->bound = def_rect(p.x,p.y,p.x,p.y);
	}
	else {
	    /* point new linepoint to old start */
	    ptr->next = pl->start;
	    /* point start to new linepoint */
	    pl->start = ptr;
	    /* update point count */
	    pl->np += 1;
	    /* reset bounds as necessary */
	    if( p.x < pl->bound.botleft.x )
		pl->bound.botleft.x = p.x;
	    if( p.x > pl->bound.topright.x )
		pl->bound.topright.x = p.x;
	    if( p.y < pl->bound.botleft.y )
		pl->bound.botleft.y = p.y;
	    if( p.y > pl->bound.topright.y )
		pl->bound.topright.y = p.y;
	}
}


/** Add a point to end of polyline.
  *
  * @param pl pointer to polyline.
  * @param p point to add to end.
  */
void
addtoend(polyline *pl, dpoint p)
{
        linepoint *ptr;

	/* create new linepoint */
	ptr = createlinepoint(p);

	if( pl->np <= 0 ) {
	    /* no points already in polyline, so add this one */
	    pl->start = ptr;
	    pl->end = ptr;
	    pl->np = 1;
	    pl->bound = def_rect(p.x,p.y,p.x,p.y);
	}
	else {
	    /* point old end linepoint to new end */
	    pl->end->next = ptr;
	    /* point end pointer to new linepoint */
	    pl->end = ptr;
	    /* update point count */
	    pl->np += 1;
	    /* reset bounds as necessary */
	    if( p.x < pl->bound.botleft.x )
		pl->bound.botleft.x = p.x;
	    if( p.x > pl->bound.topright.x )
		pl->bound.topright.x = p.x;
	    if( p.y < pl->bound.botleft.y )
		pl->bound.botleft.y = p.y;
	    if( p.y > pl->bound.topright.y )
		pl->bound.topright.y = p.y;
	}
}


/** Reverse the order of a polyline.
  *
  * @param pl polyline to reverse.
  * @return new reversed polyline.
  */
polyline *
reverseline(polyline *pl)
{
	polyline *rev;
	linepoint *ptr;

	if( pl == NULL || pl->np <= 0 )
	    quit("reverseline: NULL polyline\n");

	/* make new polyline in reverse order */
	rev = createpolyline();
	for(ptr=pl->start; ptr!=NULL; ptr = ptr->next )
	    addtostart(rev,ptr->p);
	return(rev);
}


/* delete point from polyline */
void
deletefromline(polyline *pl, dpoint p)
{
	quit("deletefromline not implemented yet\n");
}

/** Join polylines. The list of polylines must have common start
  * or end points.
  * This routine only makes 1 pass through the list, so needs to
  * be called repeatedly to form maximal length lines.
  *
  * @param list list of polylines to join.
  * @param n number of polylines in list.
  * @param gap threshold distance to test start/end points.
  * @return number of joins.
  */
long
joinsegments(polyline *list, long int n, double gap)
{
    long joins;
    long i, j;
    polyline *pl;
    polyline *new;
    linepoint *ptr;

    joins = 0;
    for(i=0; i<n; i++) {
	pl = &list[i];
	if( pl->np <= 0 )
	    continue;
	/* check if segment closed - if so, don't allow joins */
	if( pl->np > 2 && samepoint(pl->start->p, pl->end->p, gap) )
	    continue;
	for(j=i+1; j<n; j++) {
	    if( list[j].np <= 0 )
		continue;
	    /* check if segment closed - if so, don't allow joins */
	    if( list[j].np > 2 && samepoint(list[j].start->p, list[j].end->p, gap) )
		continue;
	    if( samepoint(list[j].start->p, pl->end->p, gap) ) {
		/* add points to end of pl */
		for(ptr=list[j].start->next; ptr!=NULL; ptr = ptr->next )
		    addtoend(pl, ptr->p);
		clearpolyline(&list[j]);
		joins++;
	    }
	    else if( samepoint(list[j].start->p, pl->start->p, gap) ) {
		/* add points to start of pl */
		for(ptr=list[j].start->next; ptr!=NULL; ptr = ptr->next )
		    addtostart(pl,ptr->p);
		clearpolyline(&list[j]);
		joins++;
	    }
	    else if( samepoint(list[j].end->p, pl->start->p, gap) ) {
		new = reverseline(&list[j]);
		/* add points to start of pl */
		for(ptr=new->start->next; ptr!=NULL; ptr = ptr->next )
		    addtostart(pl,ptr->p);
		clearpolyline(&list[j]);
		clearpolyline(new);
		joins++;
	    }
	    else if( samepoint(list[j].end->p, pl->end->p, gap) ) {
		new = reverseline(&list[j]);
		/* add points to end of pl */
		for(ptr=new->start->next; ptr!=NULL; ptr = ptr->next )
		    addtoend(pl, ptr->p);
		clearpolyline(&list[j]);
		clearpolyline(new);
		joins++;
	    }
	}
    }
    return(joins);
}


/** Postscript plot of polyline to file stream.
  * Assumes we are using SJW's PS header file.
  *
  * @param fp file stream. open and writable.
  * @param pl polyline.
  */
void
plotpolyline(FILE *fp, polyline pl)
{
    int n;
    linepoint *ptr;

    if( pl.np <= 1 )
	return;

    fprintf(fp,"%% Polyline, %ld points\n",pl.np);
    for(ptr=pl.start, n=0; ptr != NULL; ptr = ptr->next, n++ ) {
	if( n > 1200 ) {
	    fprintf(fp,"moredat\n");
	    n = 0;
	}
	fprintf(fp,"%e %e DP\n",ptr->p.x, ptr->p.y);
    }

    fprintf(fp,"enddat\n");
}


/** Remove spikes from a polyline.
  *
  * @param pl pointer to polyline.
  * @param tol maximum permitted jump distance.
  */
void
despike(polyline *pl, double tol)
{
    linepoint *before;
    linepoint *ptr;
    linepoint *after;
    dpoint p;

    if( pl->np < 3 )
	return;

    before = pl->start;
    for(ptr=pl->start->next; ptr != NULL; ptr = ptr->next ) {
        if( (after = ptr->next) == NULL ) break;
	if( samepoint(before->p,after->p,tol) ) {
	    before->next = after->next;
	    if( after->next == NULL )
		pl->end = before;
	    pl->np -= 2;
	    free(ptr);
	    free(after);
	    ptr = before;
	}
	before = ptr;
    }

    /* re-calculate bounding box */
    p = pl->start->p;
    pl->bound = def_rect(p.x,p.y,p.x,p.y);
    for(ptr=pl->start; ptr != NULL; ptr = ptr->next ) {
	p = ptr->p;
	if( p.x < pl->bound.botleft.x )
	    pl->bound.botleft.x = p.x;
	if( p.x > pl->bound.topright.x )
	    pl->bound.topright.x = p.x;
	if( p.y < pl->bound.botleft.y )
	    pl->bound.botleft.y = p.y;
	if( p.y > pl->bound.topright.y )
	    pl->bound.topright.y = p.y;
    }
}

/** Computes the area of a polygon.
  * The polyline is assumed closed - if not,
  * the area calculated corresponds to the area enclosed by joining
  * the start and end of the polyline
  *
  * @param pl polyline to compute area.
  * @param area of polygon (positive if polygon defined in
  * anticlockwise order).
  */
double
polyarea(polyline pl)
{
    linepoint *lp1;
    linepoint *lp2;
    double area = 0.0;

    for(lp1=pl.start; lp1!=NULL; lp1 = lp1->next) {
	lp2 = lp1->next;
	if( lp2 == NULL ) /* End of line so join back to start */
	    lp2 = pl.start;
	area -= (lp2->p.x - lp1->p.x)*(lp2->p.y + lp1->p.y)/2.0;
    }

    return(area);
}

/**
 * Calculate the Euclidian distance between two dpoints.
 */
double
dist(dpoint p, dpoint q) {
    return sqrt((p.x-q.x)*(p.x-q.x) + (p.y-q.y)*(p.y-q.y));
}

/* The implementation of polyclip below needs some auxiliary data
 * structures and functions; note that these are not declared in the
 * header file.
 */

typedef struct _waNode {
    dpoint          p;
    struct _waNode *next;
    struct _waNode *prev;
    struct _waNode *nextPoly;   /* pointer to the next polygon */
    struct _waNode *neighbor;   /* the coresponding intersection point */
    int             intersect;  /* 1 if an intersection point, 0 otherwise */
    int             entry;      /* 1 if an entry point, 0 otherwise */
    int             visited;    /* 1 if the node has been visited, 0 otherwise */
    double          alpha;      /* intersection point placemet */
} waNode;

void
wa_delete(waNode *p) {
    waNode *aux, *hold;
    if ((hold=p))                 /* note side-effecting assignment */
        do {
            aux = p;
            p = p->next;
            free(aux);
        } while(p && p!=hold);
}

/* Insert ins between first and last, using alpha to resolve ordering
 * of other waNodes inbetween */
void
wa_insert(waNode *ins, waNode *first, waNode *last) {
    waNode *aux = first;
    while (aux != last && aux->alpha < ins->alpha) aux = aux->next;
    ins->next = aux;
    ins->prev = aux->prev;
    ins->prev->next = ins;
    ins->next->prev = ins;
}

waNode *
wa_create(dpoint p, waNode *next, waNode *prev, waNode *nextPoly,
       waNode *neighbor, int intersect, int entry, int visited, double alpha) {
    waNode *new          = malloc(sizeof(waNode));
    new->p               = p;
    new->next            = next;
    new->prev            = prev;
    if (prev) prev->next = new;
    if (next) next->prev = new;
    new->nextPoly        = nextPoly;
    new->neighbor        = neighbor;
    new->intersect       = intersect;
    new->entry           = entry;
    new->visited         = visited;
    new->alpha           = alpha;
    return new;
}

/* Next non-intersecting node: */
waNode *
next_waNode(waNode *p) {
    waNode *aux = p;
    while (aux && aux->intersect) aux = aux->next;
    return aux;
}

/* Last node in the list */
waNode *
last_waNode(waNode *p) {
    waNode *aux = p;
    if (aux) while(aux->next) aux = aux->next;
    return aux;
}

/*
 * Find the starting node; first non-visited intersection (if there
 * are no intersections, it returns the first node)
 */
waNode *
wa_first(waNode *p) {
    waNode *aux = p;
    if (aux)
        do aux = aux->next;
        while (aux!=p && (!aux->intersect || (aux->intersect && aux->visited)));
    return aux;
}

void
wa_circle(waNode *p) {
    waNode *aux = last_waNode(p);
    aux->prev->next = p;
    p->prev = aux->prev;
    free(aux);
}

/* Calculate intersection point, and alpha coefficient for ordering.
 * Returns 1 if an intersection is found, 0 otherwise. */
int
wa_I(waNode *p1, waNode *p2, waNode *q1, waNode *q2,
     double *alpha_p, double *alpha_q, dpoint *pout) {
    double tp, tq, par;
    dpoint pi;

    par = (double) ((p2->p.x - p1->p.x)*(q2->p.y - q1->p.y) -
                    (p2->p.y - p1->p.y)*(q2->p.x - q1->p.x));

    if (!par) return 0;         /* parallel lines */

    tp = ((q1->p.x - p1->p.x)*(q2->p.y - q1->p.y) - (q1->p.y - p1->p.y)*(q2->p.x - q1->p.x))/par;
    tq = ((p2->p.y - p1->p.y)*(q1->p.x - p1->p.x) - (p2->p.x - p1->p.x)*(q1->p.y - p1->p.y))/par;

    if(tp<0 || tp>1 || tq<0 || tq>1) return 0;

    pi.x = p1->p.x + tp*(p2->p.x - p1->p.x);
    pi.y = p1->p.y + tp*(p2->p.y - p1->p.y);

    *alpha_p = dist(p1->p, pi) / dist(p1->p, p2->p);
    *alpha_q = dist(q1->p, pi) / dist(q1->p, q2->p);
    *pout = pi;

    return 1;
}

/* Entry or exit?  (This just toggles with each intersection node) */
int
wa_test(waNode *point, waNode *p) {
    waNode *aux, *left, i;
    dpoint pt = point->p;
    int type = 0;

    pt.x = 0;
    left = wa_create(pt, 0, 0, 0, 0, 0, 0, 0, 0.);
    for (aux=p; aux->next; aux=aux->next)
        if (wa_I(left, point, aux, aux->next, &i.alpha, &i.alpha, &i.p)) type++;
    return type%2;
}

polyline *
polylineFromWaNode(waNode *node) {
    polyline *pl  = createpolyline();
    waNode   *aux = node;
    if (aux) do {
            addtoend(pl, aux->p);
            aux = aux->next;
        } while (aux && aux != node);
    return pl;
}

/**
 * Clip a polygon using the given rectangle as a clipping window.
 * Returns the number of distinct polygons in the clip (or 0 if there
 * is no intersection), and an array of polylines constituting the
 * polygons making up the intersection (the caller is responsible for
 * freeing this array and its contents) as an output parameter.  Uses
 * the Weiler-Atherton algorithm which is a little more complex than
 * some of the alternatives, but will handle cases involve concave
 * polygons better.
 *
 * @param clipWindow The clipping rectangle
 * @param pl The polygon to clip
 * @param outputList An array of polyline pointers as an output
 *   parameter; memory will be allocated by the procedure and must be
 *   freed by the caller (both the array and its contents).
 * @return The number of polygons making up the clip region (that is,
 * the size of the outputList array after calling).
 */
int
polyclip(rect *clipWindow, polyline *pl, polyline ***outputList) {
    waNode *c, *s;
    waNode *auxs, *auxc, *is, *ic, *root = 0;
    dpoint pi;
    int e;
    double alpha_s, alpha_c;

    waNode *crt, *new, *old, *aux;
    linepoint *lp;
    int forward, i;

    int numPolys = 0;

    /* Convert clip and polyline into the richer structures needed for
     * the algorithm. */
    /* The clip window: */
    auxc = wa_create(clipWindow->botleft,                                       0,0,    0,0,0,0,0,0);
    aux  = wa_create(def_dpoint(clipWindow->botleft.x, clipWindow->topright.y), 0,auxc, 0,0,0,0,0,0);
    aux  = wa_create(clipWindow->topright,                                      0,aux,  0,0,0,0,0,0);
    aux  = wa_create(def_dpoint(clipWindow->topright.x, clipWindow->botleft.y), 0,aux,  0,0,0,0,0,0);
    c    = auxc;
    /* The polygon: */
    lp = pl->start;
    auxs = aux = wa_create(lp->p, 0,0, 0,0,0,0,0,0);
    for (lp = lp->next; lp && lp != pl->start; lp = lp->next)
        aux = wa_create(lp->p, 0,aux, 0,0,0,0,0,0);
    s = auxs;

    auxs = last_waNode(s);
    wa_create(s->p, 0, auxs, 0, 0, 0, 0, 0, 0.);
    auxc = last_waNode(c);
    wa_create(c->p, 0, auxc, 0, 0, 0, 0, 0, 0.);

    /* Calculate intersections: */
    for (auxs = s; auxs->next; auxs = auxs->next)
        if (!auxs->intersect)
            for (auxc = c; auxc->next; auxc = auxc->next)
                if (!auxc->intersect)
                    if (wa_I(auxs, next_waNode(auxs->next), auxc, next_waNode(auxc->next),
                             &alpha_s, &alpha_c, &pi)) {
                        is = wa_create(pi, 0, 0, 0, 0, 1, 0, 0, alpha_s);
                        ic = wa_create(pi, 0, 0, 0, 0, 1, 0, 0, alpha_c);
                        is->neighbor = ic;
                        ic->neighbor = is;
                        wa_insert(is, auxs, next_waNode(auxs->next));
                        wa_insert(ic, auxc, next_waNode(auxc->next));
                    }

    /* set entry/exit flags for the two polygons (These are toggled to
     * begin with; if you skip this step you will get the union, not
     * intersection): */
    e = 1-wa_test(s, c);
    for (auxs = s; auxs->next; auxs = auxs->next)
        if (auxs->intersect) {
            auxs->entry = e;
            e = 1-e;
        }

    e = 1-wa_test(c, s);
    for (auxc = c; auxc->next; auxc = auxc->next)
        if (auxc->intersect) {
            auxc->entry = e;
            e = 1-e;
        }

    /* close polygons: */
    wa_circle(s);
    wa_circle(c);

    /* The actual clipping: */
    while ((crt = wa_first(s)) != s) {
        /* Points to the chain of output polygons: */
        old = 0;
        for (; !crt->visited; crt = crt->neighbor) /* Jump to the neighbouring poly for traversal each time we hit an intersection */
            for (forward = crt->entry ;; ) {
                /* For each waNode visited, add it to the chain (only bother with forward links): */
                new = wa_create(crt->p, old, 0, 0, 0, 0, 0, 0, 0.);
                old = new;
                crt->visited = 1;
                /* always go clockwise; direction through list
                 * depends on whether we have entered or exited: */
                crt = forward ? crt->next : crt->prev;
                if (crt->intersect) {
                    crt->visited = 1;
                    break;
                }
            }

        old->nextPoly = root;
        root = old;
        numPolys++;
    }

    if (numPolys) {
        /* We've found at least one intersection polygon; package
         * it/them up into outputList and we're done. */
        *outputList = malloc((size_t)numPolys * sizeof(polyline *));
        for (i = 0, aux = root; i < numPolys; i++, aux = aux->nextPoly)
            (*outputList)[i] = polylineFromWaNode(aux);
    }
    else {
        /* If we have no insersection polygons there are three
         * possibilities: the polygon completely contains the clipping
         * window, the clipping window completely contains the
         * polygon, or they are completely distinct.  We check the
         * first two scenarios by single point-in-polygon tests; the
         * third then follows by elimination (note that with no
         * intersections, if one point is contained they all must be).
         * Note that we can just use the original rect/polyline
         * arguments for extra convenience, but must return copies so
         * the user doesn't free the same memory twice. */

        if (pointinpoly(clipWindow->botleft, *pl)) {
            /* The clip window is contained by the polygon; return the
             * clip region as a polygon: */
            *outputList  = malloc(sizeof(polyline *));
            **outputList = createpolyline();
            addtoend(**outputList, clipWindow->botleft);
            addtoend(**outputList, def_dpoint(clipWindow->botleft.x, clipWindow->topright.y));
            addtoend(**outputList, clipWindow->topright);
            addtoend(**outputList, def_dpoint(clipWindow->topright.x, clipWindow->botleft.y));
            numPolys = 1;
        }
        else if (pointinrect(pl->start->p, *clipWindow)) {
            /* The polygon is completely contained by the clip window;
             * return a copy of the polygon: */
            *outputList  = malloc(sizeof(polyline *));
            **outputList = createpolyline();
            for (lp = pl->start; lp; lp = lp->next)
                addtoend(**outputList, lp->p);
            numPolys = 1;
        }
        /* Otherwise, we just leave outputList as NULL and numPolys as 0. */
    }

    /* Clean up: */
    wa_delete(s);
    wa_delete(c);
    if (root) do {
            aux = root;
            root = root->nextPoly;
            wa_delete(aux);
        } while (root);

    return numPolys;
}
