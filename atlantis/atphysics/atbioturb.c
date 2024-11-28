/*********************************************************************

    File:           atbioturb.c

    Created:        Wed May 17 13:45:02 EST 1995

    Author:         Stephen Walker
                    CSIRO Division of Oceanography

    Purpose:        This file contains routines which represent the bioturbation processes.
		            These processes primarily move particulate tracers within the sediments.

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
#include <sjwlib.h>
#include <atlantisboxmodel.h>
#include <atUtilLib.h>

/* Prototypes for local routines */
void partic_diff(MSEBoxModel *bm, Box *bp, SedModel *sm, double **newwc, double **newsed);
void expul(MSEBoxModel *bm, Box *bp, SedModel *sm, double **newwc, double **newsed);
void mix_in(MSEBoxModel *bm, Box *bp, SedModel *sm, double **newwc, double **newsed);

/* Overall routine called by physics() */
void bioturb(MSEBoxModel *bm, double ***newwc, double ***newsed)
{
    int b = 0;

	if( verbose )
		fprintf(stderr,"Entering bioturb\n");

    /* Loop over each box */
    for(b=0; b<bm->nbox; b++) {
        Box *bp = &bm->boxes[b];
	SedModel *sm = &bm->boxes[b].sm;

	if( bp->type != BOUNDARY && bp->type != LAND) {
	    /* Diffusion - based on diffusion in bioirrigation */
	    partic_diff(bm,bp,sm,newwc[b],newsed[b]);

	    /* Expulsion */
        expul(bm,bp,sm,newwc[b],newsed[b]);

	    /* Expulsion and intake */
        mix_in(bm,bp,sm,newwc[b],newsed[b]);
	}
    }

}

/* This routine implements particulate diffusion within the sediments, for one
sediment model (1 box) */

void partic_diff(MSEBoxModel *bm, Box *bp, SedModel *sm, double **newwc, double **newsed)
{
	int k = 0;
	int n = 0;
	int tk = sm->topk;			        /* Index of top sediment layer */
	double *thick = alloc1d(sm->nz);		/* Temporary storage - thickness (amount) of
											sed to transfer*/

	/* Nothing to do if less than 2 layers */
	if( tk >= sm->nz-1){
		free1d(thick);
		return;
	}


	/* Calculate the amount to 'diffuse' (move) between layers */
	for(k=tk+1; k<sm->nz; k++) {
		sm->partic_kz[k] = bm->bt_partic_kz * sm->biodens * sm->turbenh
			* bioprofile(bm,sm->gridz[k],sm->biodepth);
	/* To convert this rate to a thickness multiply by timestep (as then
	converted to seconds ensure that timestep is always in days) and divide
	by cell thickness. This is because the diffusion coefficient is in m2 s-1
	thus to get it to m must multiply by s m-1.
	*/
		thick[k] = sm->partic_kz[k] * bm->dt * 86400.0 / sm->dz[k];
	}

	/* Loop over tracers */
	for(n=0; n<bm->ntracer; n++) {

	/* Skip tracers which are dissolved or not allowed in the sediment or
	macrobenthos (indicated by particle size).
	*/
		if( !(bm->tinfo[n].partic) || !(bm->tinfo[n].insed) || (bm->tinfo[n].psize >= 1.0) || !bm->tinfo[n].can_be_moved)
			continue;

	/* Do diffusion - using simplification of simply moving amounts in
	the following way: x-1 <-> x <-> x+1 (where <-> means movement goes both ways
	to and from x etc.). This is a simplification, but its OK wrt to the
	small amounts being moved between layers in multi-layered cases. */
		for(k=tk; k<sm->nz-1; k++){
			if( k<sm->nz )
				transfer_sed(sm,k,k+1,thick[k]);
			if( k>tk )
				transfer_sed(sm,k,k-1,thick[k]);
		}
	}

    free1d(thick);

}

void expul(MSEBoxModel *bm, Box *bp, SedModel *sm, double **newwc, double **newsed)
{
	int k = 0;
	int n = 0;
	int tk = sm->topk;			        /* Index of top sediment layer */
	double *thick = alloc1d(sm->nz);		/* Temporary storage - thickness (amount) of
											sed to transfer*/

	/* Nothing to do if less than 2 layers */
	if( tk >= sm->nz-1){
		free1d(thick);
		return;
	}

	/* Calculate the amount to 'expel' (move to surface) */
	for(k=tk+1; k<sm->nz; k++) {
		thick[k] = bm->bt_expulsion * sm->biodens * sm->turbenh
			* bioprofile(bm,sm->gridz[k],sm->biodepth) * bm->dt * 86400.0;
	/* Note that as expulsion coefficient has units m s-1 then must multiply
	by the timestep (converted from days to seconds) to get the flux as a
	thickness not a rate.*/
	}

	/* Loop over tracers */
	for(n=0; n<bm->ntracer; n++) {

	/* Skip tracers which are dissolved or not allowed in the sediment or
	macrobenthos (indicated by particle size).
	*/
		if( !(bm->tinfo[n].partic) || !(bm->tinfo[n].insed) || (bm->tinfo[n].psize >= 1.0) || !bm->tinfo[n].can_be_moved )
			continue;

	/* Do expulsion - using simplification of simply moving amounts in
	the following way: x -> 0 (where -> means movement goes one way,
	to the top most layer in this case). This is a simplification, but
	its OK wrt to the small amounts being moved between layers in
	multi-layered cases. */
		for(k=tk; k<sm->nz; k++)
			if( k>tk )
				transfer_sed(sm,k,tk,thick[k]);
	}
    free1d(thick);

}

void mix_in(MSEBoxModel *bm, Box *bp, SedModel *sm, double **newwc, double **newsed)
{

	int k = 0;
	int n = 0;
	int tk = sm->topk;			        /* Index of top sediment layer */
	double *thick = alloc1d(sm->nz);		/* Temporary storage - thickness (amount) of
											sed to transfer*/

	/* Nothing to do if less than 2 layers */
	if( tk >= sm->nz-1){
		free1d(thick);
		return;
	}

	/* Calculate the amount to 'exchange' (move to depth from surface and vice versa) */
	for(k=tk+1; k<sm->nz; k++) {
		thick[k] = bm->bt_exchange * sm->biodens * sm->turbenh
			* bioprofile(bm,sm->gridz[k],sm->biodepth) * bm->dt * 86400.0;
	/* Note that as exchange coefficient has units m s-1 then must multiply
	by the timestep (converted from days to seconds) to get the flux as a
	thickness not a rate.*/
	}

	/* Loop over tracers */
	for(n=0; n<bm->ntracer; n++) {

	/* Skip tracers which are dissolved or not allowed in the sediment or
	macrobenthos (indicated by particle size).
	*/
		if( !(bm->tinfo[n].partic) || !(bm->tinfo[n].insed) || (bm->tinfo[n].psize >= 1.0) || !bm->tinfo[n].can_be_moved )
			continue;

	/* Do exhange with surface (mixing) - using simplification of simply moving amounts in
	the following way: x <-> 0 (where <-> means movement goes both ways
	to and from x etc.). This is a simplification, but its OK wrt to the
	small amounts being moved between layers in multi-layered cases. */
		for(k=tk; k<sm->nz; k++){
			if( k>tk ) {
				transfer_sed(sm,k,tk,thick[k]);
				transfer_sed(sm,tk,k,thick[k]);
			}
		}

	}
    free1d(thick);

}
