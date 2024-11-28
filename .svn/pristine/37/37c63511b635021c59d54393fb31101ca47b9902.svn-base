/*********************************************************************

    File:           polyline.h

    Created:        Mon Aug 24 13:53:08 EST 1992

    Author:         Stephen Walker
                    CSIRO Division of Oceanography

    Purpose:        Include file for polyline data structures

    Revisions:      none

    $Id: polyline.h 2761 2011-07-01 04:35:25Z gor171 $

*********************************************************************/

typedef struct
{
	double x;
	double y;
} dpoint;

typedef struct
{
	dpoint botleft;
	dpoint topright;
} rect;

typedef struct linepoint
{
	dpoint	p;
	struct linepoint	*next;
} linepoint;

typedef struct
{
	long	np;	/* number of points */
	rect    bound;  /* bounding rectangle */
	linepoint	*start;	/* starting point */
	linepoint	*end;	/* ending point */
} polyline;

/* Prototypes */
dpoint 	*createdpoint(double x, double y);
linepoint 	*createlinepoint(dpoint p);
polyline *createpolyline(void);
void freepolyline(polyline *pl);
rect 	*createrect(double x1, double y1, double x2, double y2);
dpoint 	def_dpoint(double x, double y);
rect	def_rect(double left, double bot, double right, double top);
int	samepoint(dpoint p1, dpoint p2, double dist);
int	pointinrect(dpoint p, rect r);
int	pointinpoly(dpoint p, polyline pl);
int	pointonpolyline(dpoint p, polyline pl);
int	sectrect(rect r1, rect r2);
void	clearpolyline(polyline *pl);
void	printpolyline(FILE *fp, polyline *pl);
long	readpolyline(FILE *fp, polyline *pl);
polyline *resamplepolyline(polyline *pl, double resolution);
polyline *reverseline(polyline *pl);
polyline *smoothpolyline(polyline *pl, long ns);
void	addtostart(polyline *pl, dpoint p);
void	addtoend(polyline *pl, dpoint p);
void	deletefromline(polyline *pl, dpoint p);
long	joinsegments(polyline list[], long n, double gap);
void	plotpolyline(FILE *fp, polyline pl);
void	despike(polyline *pl, double tol);
double  polyarea(polyline pl);
double  dist(dpoint p, dpoint q);
int     polyclip(rect *clipWindow, polyline *pl, polyline ***outputList);
