/*********************************************************************

    File:           atsedprops.c

    Created:        Wed Jun 14 12:40:34 EST 1995

    Author:         Stephen Walker
                    CSIRO Division of Oceanography

    Purpose:        This file contains routines which
					calculate physical sediment properties

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

/* Prototypes for routines used here */
double calc_por(MSEBoxModel *bm, SedModel *sm, int k);
double calc_css(MSEBoxModel *bm, SedModel *sm, int k);
double calc_er(MSEBoxModel *bm, SedModel *sm, int k);

void sedprops(MSEBoxModel *bm)
{
    int b = 0;
    int k = 0;

    /* Loop over each box */
    for(b=0; b<bm->nbox; b++) {
		SedModel *sm = &bm->boxes[b].sm;

		/* Loop over empty sediment layers */
		for( k=0; k<sm->topk; k++) {
			sm->porosity[k] = 1.0;
			sm->css[k] = 0.0;
		}

		/* Loop over each sediment layer */
		for( k=sm->topk; k<sm->nz; k++) {
			/* FIX - find tracer with max concentration */
			/* FIX - maybe split into fine and coarse, and
			 * base calculations on this split
			 */
			/* FIX - find average layer density */
			/* FIX - base calculations below on those properties */
			/* Porosity */
			sm->porosity[k] = calc_por(bm,sm,k);

			/* Critical shear stress */
			sm->css[k] = calc_css(bm,sm,k);
			/* Erosion rate constant - must come after css */
			sm->er[k] = calc_er(bm,sm,k);


		}
    }
}

double calc_por(MSEBoxModel *bm, SedModel *sm, int k)
{
    int n = 0;
    double solidvol = 0.0;

    /* Loop over each tracer */
    for(n=0; n<bm->ntracer; n++) {
	/* Skip tracers which are not particulate or
	 * not in sediment or which have weirdo values
	 */
	if( !bm->tinfo[n].partic ||
	    !bm->tinfo[n].insed ||
	    bm->tinfo[n].b_dens <= 0 ||
	    isnan(sm->tr[k][n]) )
	    continue;

		/* Add volume of particles to solid volume */
		solidvol += sm->tr[k][n]*sm->volume[k]/bm->tinfo[n].b_dens;

    }

    if( sm->volume[k] <= 0.0 || solidvol == 0.0 )
	return(1.0);

    return( (sm->volume[k]-solidvol)/sm->volume[k]);

}

/* Routine to calculate critical shear stress.
 * The units of the value returned here are N m-2.
 * FIX - For the  moment,  we return a constant, which is a typical
 * value for fine sediments consolidated for about 24 hours
 * - see Teisson J Hydraul. Res. 29 1992
 */
double calc_css(MSEBoxModel *bm, SedModel *sm, int k)
{
    return(0.2);
}

/* Routine to calculate erosion rate constant.
 * The units of the value returned here are m s-1.
 * This value is similar to the empirical constant M
 * described in Teisson. However, M has units of kg m-2 s-1,
 * whereas I want m s-1 (thickness of sediment eroded per second)
 * This is obtained by dividing M by the sediment density
 *
 * The value returned here is calculated according to Delo 1988,
 * as cited in Uncles and Stephens (1989) JGR, 94 C10 14395-14405
 *
 * FIX - For the moment, the value returned assumes a deposit
 * density of 1000 kg m-3.
 *
 */
double calc_er(MSEBoxModel *bm, SedModel *sm, int k)
{
    return(0.002*sm->css[k]/1000.0);
}
