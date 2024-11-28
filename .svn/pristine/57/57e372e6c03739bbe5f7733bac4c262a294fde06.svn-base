/*********************************************************************

    File:           atresuspension.c
    
    Created:        Wed Mar  8 12:06:36 EST 1995
    
    Author:         Stephen Walker
                    CSIRO Division of Oceanography
    
    Purpose:        Implement simple resuspension for pilot Box Model
    
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

void resuspendBM(MSEBoxModel *bm, double ***newwc, double ***newsed, FILE *llogfp)
{
    int b = 0;

    /* Loop over each box */
    for(b=0; b<bm->nbox; b++) {
	Box *bp = &bm->boxes[b];
	SedModel *sm = &bm->boxes[b].sm;
	double tleft = bm->dt;
	double nxsbs = 0;
	double er_limited;
    double ask_t;

	if( verbose )
		fprintf(stderr,"Entering resuspendBM\n");

	/* Skip boundary boxes */
	if( bp->type == BOUNDARY || bp->type == LAND)
	    continue;

	/* Check that there is some sediment present */
	if( sm->topk >= sm->nz-1 && sm->dz[sm->topk] <= 0.0 )
	    /* No - skip this box */
	    continue;

	if(!bm->supplied_stress){
		/* Read the normalised excess bottom stress value.
		 * Here we assume that this represents the value
		 * appropriately integrated over the time step, and including
		 * wave effects. These calculations do not depend on the transport
		 * model as such (except for time varying critical shear stresses,
		 * making the integral over long time steps problematic, and,
		 * perhaps, for time varying bottom roughness?).
		 * For efficiency, they should be done by the Interface
		 * between the Hydrodynamic model and the Transport model.
		 * The calculations are described in Grant and Madsen 1986 
		 * Ann. Rev. Fluid Mech. 18 265-305. An implementation
		 * for a simple 1 column model was part of the work that
		 * Reg Uncles did during his visit to 
		 * CSIRO Division of Oceanography in 1994.
		 */
        
        if (bm->ts_on_hydro_time)
            ask_t = bm->hd.t;
        else
            ask_t = bm->t;

		nxsbs = tsEvalXY(bm->stress,bm->stress_id,ask_t,bp->inside.x,bp->inside.y);
		/*Store new stress values*/
		bm->boxes[b].stress = nxsbs;
	}
	bp->erosion_rate = bm->boxes[b].stress * sm->er[sm->topk];
	er_limited = bp->erosion_rate;

	/* Limit maximum erosion depth in any one time step */
	if( er_limited*tleft > bm->max_erosion )
	    er_limited = bm->max_erosion/tleft;

	/* Loop while erosion will occur and more time remains for this step */
	while( nxsbs > 0.0 && tleft > 0.0 ) {
	    int n;
	    /* Get sediment surface index */
	    int tk = sm->topk;
	    /* Erosion rate for this sediment layer (m/s) */
	    /* double er = sm->er[tk]*nxsbs; FIX - removed layer dependence SJW 12/11/96 */
	    double er = er_limited;
	    /* Amount to erode from this layer */
	    double erdz = min(er*tleft, sm->dz[tk]);
	    /* Volume eroded */
	    double ervol = erdz*bp->area;
	    /* Pore water volume */
	    double pwv = ervol*sm->porosity[tk];
	    /* New bottom water layer volume */
	    double newwcvol = bp->volume[0]+pwv;

	    /* Subtract erosion time for this layer */
		tleft -= erdz/er;

	    /* Erode tracers and calculate new water column concentrations */
	    for(n=0; n<sm->ntr; n++) {
		if( bm->tinfo[n].insed && bm->tinfo[n].partic && bm->tinfo[n].can_be_moved) {

		    /* Particulate  tracers */
		    newwc[b][0][n] = (newwc[b][0][n]*bp->volume[0] +
				      newsed[b][tk][n]*ervol)/newwcvol;

		} else if( bm->tinfo[n].insed && bm->tinfo[n].dissol && bm->tinfo[n].can_be_moved)
		    /* Pore water */
		    newwc[b][0][n] = (newwc[b][0][n]*bp->volume[0] +
				      newsed[b][tk][n]*pwv)/newwcvol;
	    }

	    /* Adjust water column volume */
	    bp->volume[0] += pwv;
	    bp->dz[0] = bp->volume[0]/bp->area;

	    /* Adjust sediment thickness */
	    sm->dz[tk] -= erdz;
	    sm->volume[tk] = sm->dz[tk]*bp->area;

	    /* If layer fully eroded, shift to next layer */
	    if( sm->dz[tk] <= 0.0 ) {
		/* Clear this layer */
		sm->dz[tk] = sm->volume[tk] = 0.0;
		for(n=0; n<bm->ntracer; n++)
		    newsed[b][tk][n] = 0.0;

		/* Are more layers present? */
		if( tk < sm->nz-1 )
		    /* Yes - change top index */
		    sm->topk += 1;
		else {
		    /* No - can't erode any more */
		    tleft = 0.0;
		    if( verbose )
			warn("resuspendBM: Sediment totally eroded, box %d, time %g\n",b,bm->t);
		}
	    }
	}
    }

    /* Re-calculate layer coordinates */
    layer_coords(bm,llogfp);
}
