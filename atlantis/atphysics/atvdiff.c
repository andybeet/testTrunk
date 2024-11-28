/*********************************************************************

    File:           atvdiff.c
    
    Created:        Tue Dec 6 11:24:32 EST 1994
    
    Author:         Stephen Walker
                    CSIRO Division of Oceanography
    
    Purpose:        This file contains routines which perform vertical diffusion.
					vdiffBMwc      - water column mixing.
    
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
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <sjwlib.h>
#include <atlantisboxmodel.h>
#include <atUtilLib.h>

/* Vertical mixing in the water column */
void vdiffBMwc(MSEBoxModel *bm, double ***newval)
{
    int b = 0;
    int k = 0;
    int n = 0;
    double dt = bm->dt;
    double *c = alloc1d(bm->wcnz);

	if( verbose )
		fprintf(stderr,"Entering vdiffBMwc\n");

    /* Loop over each box */
    for(b=0; b<bm->nbox; b++) {
        Box *bp = &bm->boxes[b];

	if( bp->type != BOUNDARY && bp->type != LAND ) {
	    /* Loop over each tracer */
	    for(n=0; n<bm->ntracer; n++) {
		/* Skip tracers not in water column */
		if( !bm->tinfo[n].inwc )
		    continue;

		/* Don't mix water variable! */
		if( n==bm->waterid )
		    continue;

		/* Collect tracer values for each layer */
		for(k=0; k<bp->nz; k++)
		    c[k] = newval[b][k][n];

		/* Set fluxes at top and bottom to zero */
		bp->kz[0] = bp->kz[bp->nz] = 0.0;

		/* Calculate diffusion */
		diffusion1d(bp->nz,c,bp->cellz,bp->kz,bp->gridz,dt,bm->a_wc);

		/**
		if(!strcmp(bm->tinfo[n].name, "Oxygen")){
			for(k=0; k<bp->nz; k++){
				fprintf(bm->logFile, "Time: %e, b%d-k%d %s start %e end %e (with kz: %e)\n",
					bm->dayt, b, k, bm->tinfo[n].name, newval[b][k][n], c[k], bp->kz[k]);
			}
		}
		**/
		/* If mixing with deep ocean allow replenishment lowest water column box */
		if(bm->mix_deep && (bm->mix_deep_depth > bm->boxes[b].botz)){

			if(!strcmp(bm->tinfo[n].name, "NH3"))
				c[0] = bp->bottNH;
			else if(!strcmp(bm->tinfo[n].name, "NO3"))
				c[0] = bp->bottNO;
            else if((!strcmp(bm->tinfo[n].name, "Oxygen")) && bm->mix_deep_O2)
                c[0] = bp->bottO2;
            else if(!strcmp(bm->tinfo[n].name, "Si"))
				c[0] = bp->bottSi;
			else if(!strcmp(bm->tinfo[n].name, "MicroNut"))
				c[0] = bp->bottFe;
			else if (bm->track_atomic_ratio && (strcmp(bm->tinfo[n].name, "Phosphorus") == 0))
				c[0] = bp->bottP;
			else if (bm->track_atomic_ratio && (strcmp(bm->tinfo[n].name, "Carbon") == 0))
				c[0] = bp->bottC;

		}

		/* Store new values */
		for(k=0; k<bp->nz; k++)
		    newval[b][k][n] = c[k];
	    }
	}
    }

    free1d(c);
}
