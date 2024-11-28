/*********************************************************************

    File:           colourtable.h
    
    Created:        Mon May 16 15:24:46 EST 1994
    
    Author:         Stephen Walker
                    CSIRO Division of Oceanography
    
    Purpose:        Include file for colourtable data structures
    
    Revisions:      none

    $Id: colourtable.h 2761 2011-07-01 04:35:25Z gor171 $

*********************************************************************/

typedef struct
{
	int n;
	double *v;
	double *r;
	double *g;
	double *b;
} colourtable;

/* Prototypes */
colourtable	*readCT(char *fname);
void		valcolour(double v, colourtable ct,
			  double *r, double *g, double *b);
