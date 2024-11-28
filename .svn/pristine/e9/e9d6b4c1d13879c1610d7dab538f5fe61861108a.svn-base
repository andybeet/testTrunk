/*

    File:           xytoij.c
    
    Created:        Fri Feb 19 11:36:24 EST 1992
    
    Author:         Daniel Delbourgo/Stephen Walker
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
    
    Revisions:      none

    $Id: xytoij.c 2761 2011-07-01 04:35:25Z gor171 $

*/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <sjwlib.h>

#define EPS 1e-6
#define iszero(x) (((x) > -EPS) && ((x) < EPS))


leaf
*createleaf(polyline *pl, long int i1, long int i2, long int j1, long int j2)
{
	leaf *petal;

	if ((petal=(leaf *)malloc(sizeof(leaf))) == NULL)
		quit("createleaf: memory allocation failure\n");
        memset(petal, 0, sizeof(leaf));
	petal->boundary = pl;
	petal->mini = i1;
	petal->maxi = i2;
	petal->minj = j1;
	petal->maxj = j2;
	petal->half1 = NULL;
	petal->half2 = NULL;

	return(petal);
}

tree
*createtree(leaf *startleaf, polyline *outline)
{
	tree *conifer;

	if ((conifer=(tree *)malloc(sizeof(tree))) == NULL)
		quit("createtree: memory allocation failure\n");
        memset(conifer, 0, sizeof(tree));
	conifer->leaves = 1;
	conifer->outline = outline;
	conifer->trunk = startleaf;
	conifer->lasti = -1;
	conifer->lastj = -1;
	conifer->lastcell = NULL;

	return(conifer);
}

void
addtoleaf(leaf *currleaf, leaf *leaf1, leaf *leaf2)
{
	currleaf->half1 = leaf1;
	currleaf->half2 = leaf2;
}

void
addtotree(tree *conifer, leaf *currleaf, leaf *leaf1, leaf *leaf2)
{
	addtoleaf(currleaf,leaf1,leaf2);
	if (leaf1 != NULL) 
		++(conifer->leaves);
	if (leaf2 != NULL)
		++(conifer->leaves);
}

long
intavrge(long int n1, long int n2)
{
	return((n1+n2)/2);
}

polyline
*formbndry(long int i1, long int i2, long int j1, long int j2, double **gx, double **gy)
{
	polyline *pl;
	dpoint point;
	long count;

	pl = createpolyline();
	for(count = i1; count <= i2; ++count) {
		point = def_dpoint(gx[j1][count],gy[j1][count]);
		addtoend(pl,point);
	}
	for(count = j1+1; count <= j2; ++count) {
		point = def_dpoint(gx[count][i2],gy[count][i2]);
		addtoend(pl,point);
	}
	for(count = i2-1; count >= i1; --count) {
		point = def_dpoint(gx[j2][count],gy[j2][count]);
		addtoend(pl,point);
	}
	for(count = j2-1; count >= j1; --count) {
		point = def_dpoint(gx[count][i1],gy[count][i1]);
		addtoend(pl,point);
	}

	return(pl);
}

void
divideleaf(leaf *currleaf, leaf **leaf1, leaf **leaf2, double **gx, double **gy)
{
	polyline *pl1;
	polyline *pl2;
	long i1;
	long i2;
	long i3;
	long i4;
	long j1;
	long j2;
	long j3;
	long j4;

	if ((currleaf->maxi==currleaf->mini+1) && (currleaf->maxj==currleaf->minj+1)) {
		*leaf1 = *leaf2 = NULL;
		return;
	}

	if (currleaf->maxi-currleaf->mini > currleaf->maxj-currleaf->minj) {
		j1 = j3 = currleaf->minj;
		j2 = j4 = currleaf->maxj;
		i1 = currleaf->mini;
		i2 = i3 = intavrge(currleaf->mini,currleaf->maxi);
		i4 = currleaf->maxi;
	}
	else {
		i1 = i3 = currleaf->mini;
		i2 = i4 = currleaf->maxi;
		j1 = currleaf->minj;
		j2 = j3 = intavrge(currleaf->minj,currleaf->maxj);
		j4 = currleaf->maxj;
	}
	pl1 = formbndry(i1,i2,j1,j2,gx,gy);
	pl2 = formbndry(i3,i4,j3,j4,gx,gy);

	*leaf1 = createleaf(pl1,i1,i2,j1,j2);
	*leaf2 = createleaf(pl2,i3,i4,j3,j4);
}

void
subdivide(tree *partition, leaf *currleaf, double **gx, double **gy)
{
	leaf *leaf1;
	leaf *leaf2;

	divideleaf(currleaf,&leaf1,&leaf2,gx,gy);
	addtotree(partition,currleaf,leaf1,leaf2);
	if (leaf1 != NULL)
		subdivide(partition,leaf1,gx,gy);
	if (leaf2 != NULL)
		subdivide(partition,leaf2,gx,gy);
}


/** Initialises a tree structure, given a pair of 2 dimensional
  * regular coordinates, to facilitate conversion from coordinate
  * to inidice space.
  *
  * @param gx array of X coordinates (of size (nce1+1)*(nce2+1)).
  * @param gy array of Y coordinates (of size (nce1+1)*(nce2+1)).
  * @param nce1 number of cells in e1 direction.
  * @param nce2 number of cells in e2 direction.
  * @return a partition tree to be use by xytoij.
  */
tree *
xytoij_init(double **gx, double **gy, int nce1, int nce2)
{
tree *partition;

	leaf *firstleaf;
	polyline *firstpl;

	firstpl = formbndry(0,nce1,0,nce2,gx,gy);
	firstleaf = createleaf(firstpl,0,nce1,0,nce2);
	partition = createtree(firstleaf, firstpl);
	subdivide(partition, firstleaf,gx,gy);

	/* The values are not copies only referenced */
	partition->gridx = gx;
	partition->gridy = gy;
	partition->nce1 = nce1;
	partition->nce2 = nce2;

	return partition;
}


/** calculates the indices (i,j) of a topologically rectangular
  * grid cell containing the point (x,y).
  *
  * @param partition a tree structure returned from xytoij_init.
  * @param x X coordinate.
  * @param y Y coordinate.
  * @param ival pointer to returned I indice value.
  * @param jval pointer to returned J indice value.
  * @return non-zero if successful.
  */
int
xytoij(tree *partition, double x, double y, int *ival, int *jval)
{
	dpoint point;
	leaf *currleaf;

/* Check to see that initialisation has been done */
	if( partition == NULL )
		quit("xytoij: initialisation of partition not done\n");

	if (ival == NULL || jval == NULL)
		return(0);

	point = def_dpoint(x,y);

/* check if point is in same cell as last point found */
	if ((partition->lasti > 0)
	   && (pointinpoly(point,*partition->lastcell))) {
		*ival = partition->lasti;
		*jval = partition->lastj;
		return(1);
	}

/* check if point is in grid outline */
	currleaf = partition->trunk;
	if( !pointinpoly(point,*currleaf->boundary) )
		return(0);

/* do the full search */
	while (currleaf->half1 != NULL) {
		if (pointinpoly(point,*(currleaf->half1->boundary)))
			currleaf = currleaf->half1;
		else
			currleaf = currleaf->half2;
	}

	*ival = (int)currleaf->mini;
	*jval = (int)currleaf->minj;
	partition->lastcell = currleaf->boundary;
	return(1);
}


/** Calculates the XY coordinate for the specified fractional (i,j)
  * within a topologically rectangular grid.
  *
  * The tranformation used to compute the coords is a forward
  * tetragonal bilinear texture mapping.
  *
  * @param partition a tree structure returned from xytoij_init.
  * @param ival I indice value.
  * @param jval J indice value.
  * @param x Pointer to returned X coordinate.
  * @param y Pointer to returned Y coordinate.
  * @return non-zero if successful.
  */
int fijtoxy(tree *partition, double ival, double jval, double *x, double *y)
{
   int status = 1;
   int i, j;
   double u, v;
   double **gx = partition->gridx;
   double **gy = partition->gridy;
   double a, b, c, d, e, f, g, h;

   /* Trim I to range 0 to nce1 */
   if (ival < 0) {
     ival = 0;
     status = 0;
   }

   if (ival > partition->nce1) {
     ival = partition->nce1-EPS;
     status = 0;
   }

   /* Trim J to range 0 to nce2 */
   if (jval < 0) {
     jval = 0;
     status = 0;
   }

   if (jval > partition->nce2) {
     jval = partition->nce2-EPS;
     status = 0;
   }

   i = (int)ival;
   j = (int)jval;
   u = ival - i;
   v = jval - j;

   a = gx[j][i] - gx[j][i+1] - gx[j+1][i] + gx[j+1][i+1];
   b = gx[j][i+1] - gx[j][i];
   c = gx[j+1][i] - gx[j][i];
   d = gx[j][i];
   e = gy[j][i] - gy[j][i+1] - gy[j+1][i] + gy[j+1][i+1];
   f = gy[j][i+1] - gy[j][i];
   g = gy[j+1][i] - gy[j][i];
   h = gy[j][i];

   *x = a*u*v + b*u + c*v + d;
   *y = e*u*v + f*u + g*v + h;

   return status;
}


/** Calculates the XY coordinate for the specified (i,j) within a
  * topologically rectangular grid.
  *
  * @param partition a tree structure returned from xytoij_init.
  * @param ival I indice value.
  * @param jval J indice value.
  * @param x Pointer to returned X coordinate.
  * @param y Pointer to returned Y coordinate.
  * @return non-zero if successful.
  */
int ijtoxy(tree *partition, int ival, int jval, double *x, double *y)
{
   int status = 1;

   /* Trim I to range 0 to nce1 */
   if (ival < 0) {
     ival = 0;
     status = 0;
   }

   if (ival > partition->nce1) {
     ival = partition->nce1;
     status = 0;
   }

   /* Trim J to range 0 to nce2 */
   if (jval < 0) {
     jval = 0;
     status = 0;
   }

   if (jval > partition->nce2) {
     jval = partition->nce2;
     status = 0;
   }

   *x = partition->gridx[jval][ival];
   *y = partition->gridy[jval][ival];

   return status;
}

/** Calculates the XY coordinate for the specified fractional (i,j)
  * within a topologically rectangular grid.
  *
  * The tranformation used to compute the indices is an inverse
  * tetragonal bilinear texture mapping.
  *
  * @param partition a tree structure returned from xytoij_init.
  * @param ival I indice value.
  * @param jval J indice value.
  * @param x Pointer to returned X coordinate.
  * @param y Pointer to returned Y coordinate.
  * @return non-zero if successful.
  */
int xytofij(tree *partition, double x, double y, double *ival, double *jval)
{
   int i, j;
   double u, v, v_denom;
   double **gx = partition->gridx;
   double **gy = partition->gridy;
   int status = xytoij(partition, x, y, &i, &j);

   if (status > 0) {
      double a = gx[j][i] - gx[j][i+1] - gx[j+1][i] + gx[j+1][i+1];
      double b = gx[j][i+1] - gx[j][i];
      double c = gx[j+1][i] - gx[j][i];
      double d = gx[j][i];
      double e = gy[j][i] - gy[j][i+1] - gy[j+1][i] + gy[j+1][i+1];
      double f = gy[j][i+1] - gy[j][i];
      double g = gy[j+1][i] - gy[j][i];
      double h = gy[j][i];

      double A = a*f -b*e;
      double B = e*x - a*y + a*h - d*e + c*f - b*g;
      double C = g*x - c*y + c*h -d*g;

      u = (iszero(A)) ? -C/B : (-B - sqrt(B*B -4*A*C))/(2*A);
      v_denom = a*u + c;
      v = (iszero(v_denom)) ? (y - h)/g : (x - b*u - d)/v_denom;

      *ival = i + u;
      *jval = j + v;
   }

   return status;
}
