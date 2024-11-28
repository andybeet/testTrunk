/*********************************************************************

    File:           atboundary.c
    
    Created:        Wed Mar  1 09:06:49 EST 1995
    
    Author:         Stephen Walker
                    CSIRO Division of Oceanography
    
    Purpose:        Boundary routines for general box model
    
    Arguments:      See each routine below

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
#include <stdlib.h>
#include <math.h>
#include <sjwlib.h>
#include <atlantisboxmodel.h>

/* For now, all that the boundary code does is preserve the
 * values in the boundary boxes. This means that the boundary
 * values for all tracers stays constant at the value that
 * was read at the start of the model run.
 */
void boundaries(MSEBoxModel *bm, double ***newwctr, double ***newsedtr,  double ***newicetr, double **newlandtr,  FILE *llogfp)
{
    int b = 0;
    int k = 0;
    int n = 0;

    /* Loop over each box */
    for(b=0; b<bm->nbox; b++) {
	/* Check for boundary flag */
		if( bm->boxes[b].type == BOUNDARY ) {
			Box *bp = &bm->boxes[b];
	
			/* FIX - this code may be inadequate, as I need to worry about
			* volume, dz, sediment properties etc
			*/

			/* Loop over each water column layer and tracer */
			for(k=0; k<bp->nz; k++){
				for(n=0; n<bm->ntracer; n++){
					newwctr[b][k][n] = bp->tr[k][n];
				}
			}
			
			/* Loop over each sediment layer and tracer */
			for(k=0; k<bp->sm.nz; k++){
				for(n=0; n<bm->ntracer; n++){
					newsedtr[b][k][n] = bp->sm.tr[k][n];
				}
			}

			if(bm->ice_on == TRUE){
				/* Loop over each ice layer and tracer */
				for(k=0; k<bp->ice.nz; k++){

					for(n=0; n<bm->ntracer; n++){
						newicetr[b][k][n] = bp->ice.tr[k][n];
					}
				}
			}

			if(bm->terrestrial_on){
				/* Loop over each land tracer */
				for(n=0; n<bm->nland; n++){
					newlandtr[b][n] = bp->land.tr[n];
				}
			}
		}

	}
}
