/*********************************************************************

 File:           atsaturation.c

 Created:        Fri Dec 5 2:23:34 EST 2003

 Author:         Beth Fulton
 CSIRO Division of Marine Research

 Purpose:        This file contains routines which check
 that gases and nutrients haven't surpassed saturation values.

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
#include <string.h>
#include <math.h>
#include <sjwlib.h>
#include <atlantisboxmodel.h>

/* Checking for saturation of oxygen and nitrogen in the model layers */

void Saturation_Check(MSEBoxModel *bm, double ***newwc, double ***newsed) {
	int b, n, k, do_break;
	double DIN, NH = 0, NO = 0;

	if (verbose)
		fprintf(stderr, "Entering Saturation_Check\n");

	/* Loop over each box */
	for (b = 0; b < bm->nbox; b++) {
		Box *bp = &bm->boxes[b];
		SedModel *sm = &bm->boxes[b].sm;

		if (bp->type != BOUNDARY && bp->type != LAND) {
			/* Check for saturation of oxygen levels */
			/* Loop over each tracer */
			for (n = 0; n < bm->ntracer; n++) {
				/* Find oxygen */
				if (strcmp(bm->tinfo[n].name, "Oxygen") == 0) {
					/* Loop through layers and check for saturation */
					for (k = 0; k < bp->nz; k++) {
						if (newwc[b][k][n] > 8000.0)
							newwc[b][k][n] = 8000.0;
					}
					for (k = 0; k < sm->nz; k++) {
						if (newsed[b][k][n] > 8000.0)
							newsed[b][k][n] = 8000.0;
					}
					break;
				}
			}

			/* Check Nutrients in the watercolumn */
			for (k = 0; k < bp->nz; k++) {
				DIN = 0.0;
				do_break = 0;
				NH = 0;
				NO = 0;
				for (n = 0; n < bm->ntracer && do_break < 2; n++) {
					/* Find ammonium and nitrate */
					if (strcmp(bm->tinfo[n].name, "NH3") == 0) {
						NH = newwc[b][k][n];
						do_break++;
					} else if (strcmp(bm->tinfo[n].name, "NO3") == 0) {
						NO = newwc[b][k][n];
						do_break++;
					}
				}
				DIN = NH + NO;

				if (DIN > 12500.0) {
					//					printf("day: %e, box: %d-%d start DIN: %e (NH: %e NO: %e)", bm->dayt, b, k, DIN, NH, NO);
					do_break = 0;
					for (n = 0; n < bm->ntracer && do_break < 2; n++) {
						/* Correct ammonium and nitrate */
						if (strcmp(bm->tinfo[n].name, "NH3") == 0) {
							newwc[b][k][n] = 12500.0 * NH / DIN;
							NH = newwc[b][k][n];
							do_break++;
						} else if (strcmp(bm->tinfo[n].name, "NO3") == 0) {
							newwc[b][k][n] = 12500.0 * NO / DIN;
							NO = newwc[b][k][n];
							do_break++;
						}
					}
					DIN = NH + NO;
					//					printf(" end DIN: %e (NH: %e NO: %e)\n", DIN, NH, NO);
				}
			}

			/* Check Nutrients in the sediments */
			for (k = 0; k < sm->nz; k++) {
				DIN = 0.0;
				do_break = 0;
				NH = 0;
				NO = 0;
				for (n = 0; n < bm->ntracer && do_break < 2; n++) {
					/* Find ammonium and nitrate */
					if (strcmp(bm->tinfo[n].name, "NH3") == 0) {
						NH = newsed[b][k][n];
						do_break++;
					} else if (strcmp(bm->tinfo[n].name, "NO3") == 0) {
						NO = newsed[b][k][n];
						do_break++;
					}
				}
				DIN = NH + NO;

				if (DIN > 12500.0) {
					do_break = 0;
					for (n = 0; n < bm->ntracer && do_break < 2; n++) {
						/* Correct ammonium and nitrate */
						if (strcmp(bm->tinfo[n].name, "NH3") == 0) {
							newsed[b][k][n] = 12500.0 * NH / DIN;
							do_break++;
						} else if (strcmp(bm->tinfo[n].name, "NO3") == 0) {
							newsed[b][k][n] = 12500.0 * NO / DIN;
							do_break++;
						}
					}
				}
			}
		}
	}

	return;
}

