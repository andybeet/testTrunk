/*

    File:           ppbfetch.c
    
    Created:        Thu Jul 18 20:10:01 EST 1996
    
    Author:         Stephen Walker
                    CSIRO Division of Oceanography
    
    Revisions:      none

    $Id: ppbfetch.c 2761 2011-07-01 04:35:25Z gor171 $

*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sjwlib.h>

static char fetchname[] = "/home/walker/ppb/waves/fetch.nc";

/** Look up fetch data for Port Phillip Bay, as calculated by
  * Chris Sherwood.
  *
  * @param useuv 0 -> use dir value as passed,
  * 		 1 -> calculate dir from u and v.
  * @param dir direction of wind (FROM)
  * @param u East cpt of wind velocity
  * @param v North cpt of wind velocity
  * @param x Easting of point where fetch is wanted
  * @param y Northing of point where fetch is wanted
  * @return fetch value in metres.
  */
double
ppbfetch(char *folderPath, int useuv, double dir, double u, double v, double x, double y)
{
    static TimeSeries fts;
    static int fid = -1;
    double f;

    if( fid == -1 ) {
	/* Read fetch data */
	tsRead(folderPath, fetchname,&fts);
	if( (fid=tsIndex(&fts,"fetch")) < 0 )
	    quit("ppbfetch: Data file %s doesn't contain fetch variable\n",
			    fetchname);
    }

    if( useuv )
	/* Calculate direction from u and v */
	dir = 270 - atan2(v,u)*180/M_PI;
    
    dir = ((int)(dir)%360);
    
    f = tsEvalXY(&fts,fid,dir,x,y);
    return( f );
}
