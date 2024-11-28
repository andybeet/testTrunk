/*

    File:           colourtable.c

    Created:        Mon May 16 13:31:25 EST 1994

    Author:         Stephen Walker
                    CSIRO Division of Oceanography

    Purpose:        Colour table routines

    Arguments:      See routines below

    Returns:

    Revisions:

    $Id: colourtable.c 3369 2012-08-29 06:16:46Z gor171 $

  */

#include <stdio.h>
#include <stdlib.h>
#include "sjwlib.h"



/** Opens and reads a colour table file, and returns an
  * instance of a colourtable structure.
  * @param fname colourtable file name.
  * @return pointer to a colourtable structure.
  */
colourtable *readCT(char *fname)
{
    FILE *fp;
    int i;
    int ntvals;
    char line[MAXLINELEN];
    colourtable *p;

    if( (fp=fopen(fname,"r"))==NULL )
        return(NULL);

    /* Count the number of table entries */
    ntvals = 0;
    while( fgets(line,MAXLINELEN,fp) != NULL )
       if(line[0]!='#' && line[0]!=0)
	   ntvals++;
    if( ntvals < 1 ) {
    	fclose(fp);
	return(NULL);
    }

    /* Allocate memory for table */
    if( (p=(colourtable *)malloc(sizeof(colourtable)))==NULL ) {
    	fclose(fp);
        return(NULL);
    }
    p->v = d_alloc1d(ntvals);
    p->r = d_alloc1d(ntvals);
    p->g = d_alloc1d(ntvals);
    p->b = d_alloc1d(ntvals);

    p->n = ntvals;

    /* Read table */
    fseek(fp,0L,0);
    for(i=0; fgets(line,MAXLINELEN,fp) != NULL; ) {
        double v, r, g, b;

        if(line[0]!='#' && line[0]!=0) {
	    if( sscanf(line,"%lf %lf %lf %lf",&v,&r,&g,&b)!= 4 )
		quit("readCT: Can't read table values\n");
            p->v[i] = v;
            p->r[i] = r;
            p->g[i] = g;
            p->b[i] = b;
	    i++;
	}
    }

    /* Check values are monotonic increasing */
    for(i=1; i<ntvals; i++)
        if( p->v[i] <= p->v[i-1] )
	    quit("readCT: Table out of order\n");

    /* Close file and return colour table pointer */
    fclose(fp);
    return(p);
}

void writeCT(void)
{
    /* NOT IMPLEMENTED YET */
}


/** Get an RGB colour from the colourtable.
  *
  * Get the colour (r, g, b values) corresponding to
  * a particular value from the colour table.
  *
  * @param v A value within the colourtable range.
  * @param ct Colour table.
  * @param r Red corresponding to specified value.
  * @param g Green corresponding to specified value.
  * @param b Blue corresponding to specified value.
  */
void
valcolour(double v, colourtable ct, double *r, double *g, double *b)
{
    double frac;
    int ilow;
    int imid;
    int ihigh;

    /* first check whether v is within the table range */
    if( v <= ct.v[0] ) {
    	*r = ct.r[0];
    	*g = ct.g[0];
    	*b = ct.b[0];
    	return;
    }
    if( v >= ct.v[ct.n-1] ) {
    	*r = ct.r[ct.n-1];
    	*g = ct.g[ct.n-1];
    	*b = ct.b[ct.n-1];
    	return;
    }

    /* perform binary chop to determine values either side of v */
    ilow = 0;
    ihigh = ct.n-1;
    while( ihigh-ilow > 1 ) {
	    imid = (ilow+ihigh)/2;
	    if( v >= ct.v[imid] )
		    ilow = imid;
	    else
		    ihigh = imid;
    }

    /* Calculate fractional position */
    frac = (v-ct.v[ilow])/(ct.v[ihigh]-ct.v[ilow]);

    if( ct.n >= 100 ) {
	/* High-resolution colour table - return the nearest colour.
	 * This allows an application to use a known set of colours
	 */
	int index = (frac <= 0.5) ? ilow : ihigh;
	*r = ct.r[index];
	*g = ct.g[index];
	*b = ct.b[index];
    }
    else {
	/* Interpolate red, green and blue values */
	*r = ct.r[ilow]*(1.0-frac) + ct.r[ihigh]*frac;
	*g = ct.g[ilow]*(1.0-frac) + ct.g[ihigh]*frac;
	*b = ct.b[ilow]*(1.0-frac) + ct.b[ihigh]*frac;
    }

    return;
}
