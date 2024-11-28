/*********************************************************************

 File:           atgas.c

 Created:        Tue Mar 28 11:18:09 EST 1995

 Author:         Stephen Walker
 CSIRO Division of Oceanography

 Purpose:        Routines to implement gas exchange at the water surface.

 It is  assumed that the gas is freely available
 in the atmosphere as an 'infinite' source.
 The surface flux F is given by

 (G-S)
 F = -K -----
 dz
 where K is the diffusion coefficient (molecular diffusion), G is the gas
 concentration in the water, S is the saturation concentration, and dz
 is a small stagnant layer thickness. Both K and S will in general vary with
 temperature. dz varies with wind speed, oil slicks etc.

 For the time being, these values are held constant.

 The rate of change of concentration in the surface layer is then

 dG        (G-S)  1
 --  =  -K -----  -
 dt         dz    D

 where D is the thickness of the surface layer.
 This equation implies that G exponentially approaches S (it is
 analogous to a decay equation in (G-S), assuming that S remains
 constant over a time step).


 between

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

/* Indices for gas tracers */
//static int Oxygen_i;

/* Data structure for gases */
struct {
	char *name; /* Name of gas tracer */
	int *index; /* pointer to gas tracer index */
	double K; /* Diffusion cfft */
	double sat; /* Saturation value */
	double dz; /* stagnant layer thickness */

} gaslist[] =
	{
		{ "Oxygen", &Oxygen_i, 2e-9, 8000.0, 2e-5 },
		{ "CO2", &CO2_i, 2e-9, 0.64, 2e-5 },
		{ NULL, NULL, 0.0, 0.0, 0.0 } };

/* Prototypes for routines used here */

void gasExchange(MSEBoxModel *bm, double ***newwc) {
	int b = 0;
	static int first = 1;

	if (verbose)
		fprintf(stderr, "Entering gasExchange\n");

	/* Initialisation - Loop to find tracer indices for each gas */
	if (first) {
		int i = 0;
		int id = 0;

		for (i = 0; gaslist[i].name != NULL; i++) {
			id = trIndex(bm->tinfo, bm->ntracer, gaslist[i].name);
			if (id < 0)
				warn("gasExchange: Can't find variable %s\n", gaslist[i].name);
			*(gaslist[i].index) = id;
		}
		first = 0;
	}

	/* Loop over each box */
	for (b = 0; b < bm->nbox; b++) {
		Box *bp = &bm->boxes[b];
		double layerdz;
		int i;

		if (bp->type != BOUNDARY && bp->type != LAND) {
			/* Calculate surface layer depth */
			layerdz = bp->gridz[bp->nz] - bp->gridz[bp->nz - 1];

			/* Loop over each gas */
			for (i = 0; gaslist[i].name != NULL; i++) {
				int n = *(gaslist[i].index);
				double sat;
				double dz;
				double diff;
				double r;

				if (n >= 0) {
					/* Calculate saturation values etc
					 * FIX - these are constants for now
					 */
					sat = gaslist[i].sat;
					dz = gaslist[i].dz;
					diff = gaslist[i].K;
					r = diff / (dz * layerdz + small_num);

					/* Calculate new gas concentration */
					newwc[b][bp->nz - 1][n] = sat + decay_exact(newwc[b][bp->nz - 1][n] - sat, r, bm->dt);
				}
			}
		}
	}
}
