/*********************************************************************

    File:           athdiff.c
    
    Created:        Tue Jan 16 11:23:32 EST 2006
    
    Author:         Beth Fulton
                    CSIRO Division of Marine and Atmospheric Research
    
    Purpose:        This file contains routines which perform horizontal diffusion.
    
    Arguments:      See below

    Returns:        void
    
    Revisions:      
*********************************************************************/

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <sjwlib.h>
#include <atlantisboxmodel.h>

/* Horizontal mixing in the water column */
void hdiffBMwc(MSEBoxModel *bm, double ***newval, FILE *llogfp)
{
    int b = 0;
	int nb = 0;
	int bb = 0;
    int k = 0;
    int n = 0;
	int nnconn = 0;
    double dt = bm->dt;

	if( verbose )
		fprintf(stderr,"Entering hdiffBMwc\n");

    /* Loop over each box */
    for(b=0; b<bm->nbox; b++) {
        Box *bp = &bm->boxes[b];

	if( bp->type != BOUNDARY && bp->type != LAND) {
	    /* Loop over each tracer */
        
	    for(n=0; n<bm->ntracer; n++) {
		/* Skip tracers not in water column */
		if( !bm->tinfo[n].inwc )
		    continue;

		/* Don't mix water variable! */
		if( n==bm->waterid )
		    continue;

		/* Collect tracer values for doing diffusion for each layer in each adjoining box */
        nnconn = bm->boxes[b].nconn;
		for(k=0; k<bp->nz; k++){
			for(nb=0; nb<bm->boxes[b].nconn; nb++){
				bb = bm->boxes[b].ibox[nb];
                bm->hc[nb] = newval[bb][k][n];
                bm->hx[nb] = bm->boxes[bb].inside.x;
			}

			/* Calculate diffusion */
			diffusion1d(nnconn,bm->hc,bm->hx,bm->hk,bm->hx,dt,bm->a_wc);

			/* Store new values */
			for(nb=0; nb<bm->boxes[b].nconn; nb++){
				bb = bm->boxes[b].ibox[nb];
                newval[bb][k][n] = bm->hc[nb];
			}
		}
		}
    }
	}
}

/* Horizontal mixing in the water column done as simple weighted averaging */
void filler_hdiffBMwc(MSEBoxModel *bm, double ***newval, FILE *llogfp)
{
    int b = 0;
	int nb = 0;
	int bb = 0;
    int k = 0;
    int n = 0;
	int nnconn = 0;
    double dt = bm->dt;
    double flushin_days, flushout_days, totvol, tot_exvol, new_avg_val, added_tr, removed_tr, hmix;
    double orig_newval;
    
	if( verbose )
		fprintf(stderr,"Entering filler_hdiffBMwc\n");
    
    /* Loop over each box */
    for(b=0; b<bm->nbox; b++) {
        Box *bp = &bm->boxes[b];
        
        if( bp->type != BOUNDARY && bp->type != LAND) {
            /* Loop over each tracer */
            
            for(n=0; n<bm->ntracer; n++) {
                /* Skip tracers not in water column */
                if( !bm->tinfo[n].inwc )
                    continue;
                
                if ((bm->do_availflag && !bm->tinfo[n].flagid) || (bm->tinfo[n].partic && !bm->tinfo[n].passive)) {
                    continue;  // As not allow to move it
                } else
                    // Can do these so let the code flow
                    
                /* Don't mix water variable! */
                    if( n==bm->waterid )
                        continue;
                
                /* Collect tracer values for doing diffusion for each layer in each adjoining box */
                for(k=0; k<bp->nz; k++){
                    for(nb=0; nb<bm->max_nconn; nb++){
                        bm->hc[nb] = 0.0;
                    }
                    for(nb=0; nb<bm->nbox; nb++){
                        bm->hdone[nb] = 0;
                    }
                    
                    flushin_days = bm->dayt - bp->flush_in[k];
                    flushout_days = bm->dayt - bp->flush_out[k];
                    if (bm->fill_zero_exchange && ((flushin_days < bm->flush_threshold) && (flushout_days < bm->flush_threshold)))
                        continue;  // Don't do it if flushing frequent enough using the transport model
                    
                    orig_newval = newval[b][k][n];
                    /* Calculate exchanges */
                    nnconn = 0;
                    tot_exvol = 0.0;
                    for(nb=0; nb<bm->boxes[b].nconn; nb++){
                        bb = bm->boxes[b].ibox[nb];
                        if (!bm->hdone[bb] && (bm->boxes[bb].nz > k)) {
                            bm->hc[nnconn] = newval[bb][k][n];
                            
                            if(bm->use_fill_horizmix) {
                                hmix = bm->boxes[bb].horizmix;
                            } else {
                                hmix = 1.0;
                            }
                            
                            bm->hx[nnconn] = bm->boxes[bb].volume[k] * bm->wc_kz * hmix * dt;
                            tot_exvol += bm->boxes[bb].volume[k] * bm->wc_kz * hmix * dt;
                            nnconn++;
                            
                        }
                        bm->hdone[bb] = 1;
                    }
                    // Add on local box
                    bm->hc[nnconn] = newval[b][k][n];
                    
                    if(bm->use_fill_horizmix) {
                        hmix = bm->boxes[b].horizmix;
                    } else {
                        hmix = 1.0;
                    }
                    bm->hx[nnconn] = bm->boxes[b].volume[k] * bm->wc_kz * hmix * dt;
                    totvol = bm->boxes[b].volume[k] * hmix * bm->wc_kz * dt;  // As dictated by this box not all boxes or volumes get swamped
                    nnconn++;
                    
					/* Rescale the exchange volumes (to proportions) as final flow dictated
                     by the recieving box so don't blow out the volumes */
					for(nb=0; nb<nnconn - 1; nb++){
						bm->hx[nb] /= tot_exvol;
					}
					/* Calculate diffusion */
                    added_tr = 0.0;
                    for(nb=0; nb<nnconn - 1; nb++){  // Flows into system
                        added_tr += bm->hc[nb] * bm->hx[nb] * totvol;
                    }
                    removed_tr = bm->hc[nb] * bm->hx[nb];  // Final flow is out of the system
                    
                    /* Store new values */
                    for(nb=0; nb<bm->nbox; nb++){
                        bm->hdone[nb] = 0;
                    }
                    nnconn = 0;
                    for(nb=0; nb<bm->boxes[b].nconn; nb++){
                        bb = bm->boxes[b].ibox[nb];
                        if (!bm->hdone[bb] && (bm->boxes[bb].nz > k)) {
                            
							// New value is current value  + amount entering - amount exiting
                            newval[bb][k][n] = (newval[bb][k][n] * bm->boxes[bb].volume[k] + bm->hx[nnconn] * totvol * (newval[b][k][n] - bm->hc[nnconn])) / bm->boxes[bb].volume[k];
                            
                            /*
                             if(bm->checkbox == b)
                             fprintf(bm->logFile, "Time: %e n: %s b: %d bb: %d layer: %d nnconn: %d newval: %e volume: %e hx: %e totvol: %e newvalb: %e hc: %e\n", bm->dayt, bm->tinfo[n].name, b, bb, k, nnconn, newval[bb][k][n], bm->boxes[bb].volume[k], bm->hx[nnconn], totvol, newval[b][k][n], bm->hc[nnconn]);
                             */
                            
                            nnconn++;
                            
                        }
                        bm->hdone[bb] = 1;
                    }
					//orig_newval = newval[b][k][n];
                    new_avg_val = (newval[b][k][n] * bm->boxes[b].volume[k] + added_tr - removed_tr) / bm->boxes[b].volume[k];
                    
                    /**
                     if(bm->checkbox == b)
                     fprintf(bm->logFile, "Time: %e box: %d k: %d n: %s old_newval= %e, new_avg_val= %e, added_tr= %e, removed_tr= %e vol= %e\n", bm->dayt, b, k, bm->tinfo[n].name, orig_newval, new_avg_val, added_tr, removed_tr, bm->boxes[b].volume[k]);
                     
                     **/

                    
                    newval[b][k][n] = new_avg_val;  // Store local value
                    
                    if (newval[b][k][n] < 0) {
                        warn("hdiff transport: %s negative (box: %d k: %d orig_newval= %e, new_avg_val= %e, added_tr= %e, removed_tr= %e vol= %e)\n",
                             bm->tinfo[n].name, b, k, orig_newval, new_avg_val, added_tr, removed_tr, bm->boxes[b].volume[k]);
                       /* fprintf(llogfp,
                                "hdiff transport: %s negative (box: %d k: %d orig_newval= %e, new_avg_val= %e, added_tr= %e, removed_tr= %e vol= %e)\n",
                                bm->tinfo[n].name, b, k, orig_newval, new_avg_val, added_tr, removed_tr, bm->boxes[b].volume[k]);*/
                        
                        
                        newval[b][k][n] = 0;
                    }
                }
            }
        }
	}
}

