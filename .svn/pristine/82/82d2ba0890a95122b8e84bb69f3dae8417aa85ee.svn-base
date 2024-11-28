/*********************************************************************

    File:           atprofile.c
    
    Created:        Thu May 18 16:48:03 EST 1995
    
    Author:         Stephen Walker
                    CSIRO Division of Oceanography
    
    Purpose:        This file contains routines to
					calculate various functional forms.
					They are used for example, to give values which
					represent how bioturbation declines with
					depth in the sediments.

    Arguments:      See below

    Returns:        void
    
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
#include <math.h>
#include <sjwlib.h>
#include <atlantisboxmodel.h>

double bioprofile(MSEBoxModel *bm, double z, double scale)
/* z          depth at which to evaluate function
   scale      length scale for function */
{
    switch(bm->biosedprofile) {
	case 'c': 
		return( w_tophat(z,scale) );
	case 'l': 
		return( w_linear(z,scale) );
	case 'p': 
		return( w_parabolic(z,scale) );
	case 'g': 
		return( w_gaussian(z,scale) );
	default : quit("bioprofile: Unknown profile (%c)\n",bm->biosedprofile);
    }

    /* NOTREACHED */
    return(0.0);
}
