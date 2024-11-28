/*******************************************************************\\**
\file
\brief Routines to find various indices in Box Model.
    File:           atindex.c

    Created:        Wed Dec 14 11:10:38 EST 1994

    Author:         Stephen Walker
                    CSIRO Division of Oceanography

    Purpose:        Routines to find various indices in Box Model.

    Arguments:      See routines below

    Returns:        int index value if index is found, or -1 if not found

    Revisions:      8/8/2004 EA Fulton
					Ported across the code from the southeast (sephys) model

					17/11/2004 EA Fulton
					Converted original routine definitions from
					void
					routine_name(blah,blahblah)
					int blah;
					double blahblah;

					to

					void routine_name(int blah, double blahblah)
*********************************************************************/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <sjwlib.h>
#include <atlantisboxmodel.h>


/*******************************************************************//**
Routines to allow external libraries to access indexing routines
without knowing details of box structure etc.
*********************************************************************/

/* Static pointer to box model */
/*TODO: Remove global box reference */
/*TODO: Functions should return once index is found , not complete the loop. Check with beth that
you don't need to find the last matching value. Thats teh only reason why you wouldn't return once match found.
If this is the case a comment should be added to indicate that not returning is a decision not a bug. */
MSEBoxModel *this_bm = NULL;

void setIndexBM(MSEBoxModel *bm)
{
    this_bm = bm;
}

int ext_trIndex(void *dummy, char *name)
{
    if( !this_bm )
		quit("ext_trIndex: NULL box model pointer\n");
    return( trIndex(this_bm->tinfo, this_bm->ntracer, name) );
}

int ext_xyztoijk(void *dummy, double x, double y, double z, int *i, int *j, int *k)
{
    int b = 0;
    int kk = 0;

    if( !this_bm )
		quit("ext_xyztoijk: NULL box model pointer\n");

    b = xytobox(x, y, this_bm);

    kk = wczbtok(z, &this_bm->boxes[b]);

	if( b < 0 || kk < 0 ){
		return(-1);
	}else {
		*i = b;
		*j = b;
		*k = kk;
    }
    return(0);
}

/*******************************************************************//**
Routine to find tracer index
    @param tinfo - array of TracerInfo structures
    @param n     - number of elements in array
    @param name  - name of variable to search for
*********************************************************************/
int trIndex(TracerInfo *tinfo, int n, char *name)
{
    int i = 0;
    int index = 0;

    index = -1;
    for(i=0; i<n && index<0; i++) {
    	if( strcmp(name,tinfo[i].name) == 0 ){
		  index = tinfo[i].n;
		  return index;
		}
    }

    return(index);
}


/*******************************************************************//**
Routine to find epi variable index
    einfo - array of EpiInfo structures
    n     - number of elements in array
    name  - name of variable to search for
*********************************************************************/
int epiIndex(EpiInfo *einfo, int n, char *name)
{
    int i = 0;
    int index = 0;

    index = -1;
    for(i=0; i<n && index<0; i++) {
		if( strcmp(name,einfo[i].name) == 0 )
		   index = einfo[i].n;
    }

    return(index);
}

/*******************************************************************//**
Routine to find land variable index
    einfo - array of landInfo structures
    n     - number of elements in array
    name  - name of variable to search for
*********************************************************************/
int landIndex(EpiInfo *landinfo, int n, char *name)
{
    int i = 0;
    int index = 0;

    index = -1;
    for(i=0; i<n && index<0; i++) {
		if( strcmp(name,landinfo[i].name) == 0 )
		   index = landinfo[i].n;
    }

    return(index);
}

/*******************************************************************//**
Routine to find diagnostic variable index
    @param dinfo - array of DiagInfo structures
    @param n     - number of elements in array
    @param name  - name of variable to search for
*********************************************************************/
int diagIndex(DiagInfo *dinfo, int n, char *name)
{
    int i = 0;
    int index = 0;

    index = -1;
    for(i=0; i<n && index<0; i++) {
		if( strcmp(name,dinfo[i].name) == 0 )
		  index = dinfo[i].n;
    }

    return(index);
}

/*******************************************************************//**
Routine to find diagnostic variable index
    @param finfo - array of FstatInfo structures
    @param n     - number of elements in array
    @param name  - name of variable to search for
*********************************************************************/
int fstatIndex(FstatInfo *finfo, int n, char *name)
{
    int i = 0;
    int index = 0;

    index = -1;
    for(i=0; i<n && index<0; i++) {
		if( strcmp(name,finfo[i].name) == 0 )
		    index = finfo[i].n;
    }

    return(index);
}
/*******************************************************************//**
Routine to find box index from point location
    @param x    - x coordinate of point
    @param y    - y coordinate of point
*********************************************************************/
int xytobox(double x, double y, MSEBoxModel *bm)
{
    int i = 0;
    int index = 0;
    dpoint p;

    p.x = x;
    p.y = y;
    index = -1;
    for(i=0; (i<bm->nbox) && (index<0); i++) {
		if( pointinpoly(p,*bm->boxes[i].bnd)  ){
		    index = i;
		}
    }
    return(index);
}


/*******************************************************************//**
Routine to find water column layer index from z coordinate
and box number.

    @param z    - z coordinate of point
*********************************************************************/
int wczbtok(double z, Box *bp)
{
    int k = 0;
    int index = -1;

    /* Truncate values below bottom or above top - this
     * is necessary so that point sources don't end up out
     * of the water during a 'hot start'.
     */
    if( z < bp->gridz[0] )
		return(0);
    if( z >= bp->gridz[bp->nz] ){
		return(bp->nz - 1);
    }

    /* Loop through layers to find index */
    for(k=0; (k<bp->nz) && (index<0); k++) {
		if( (z>=bp->gridz[k]) && (z<bp->gridz[k+1]) ){
		    index = k;
		}
    }
    return(index);
}


/*******************************************************************//**
 Routine that returns the wc surface for a specified box.

*********************************************************************/
double wcsurface(Box *bp)
{
    return(bp->gridz[0]);
}


/*******************************************************************//**
 Routine that returns the wc bottom for a specified box.

*********************************************************************/
double wcbottom(Box *bp)
{
    return(bp->gridz[bp->nz]);
}


/*******************************************************************//**
Routine to find sediment column layer index from z coordinate
and box number.

   @param z    - z coordinate of point
*********************************************************************/
int sedzbtok(double z, Box *bp)
{
    int k = 0;
    int index = 0;

    index = -1;
    for(k=bp->sm.topk; (k<bp->sm.nz) && (index<0); k++) {
		if( (z<=bp->sm.gridz[k]) && (z>bp->sm.gridz[k+1]) )
		    index = k;
    }
    return(index);
}


/*******************************************************************//**
 Routine that returns the sediment surface for a specified box.

*********************************************************************/
double sedsurface(Box *bp)
{
    return(bp->sm.gridz[bp->sm.topk]);
}


/*******************************************************************//**
 Routine that returns the sediment bottom for a specified box.

*********************************************************************/
double sedbottom(Box *bp)
{
    return(bp->sm.gridz[bp->sm.nz]);
}
