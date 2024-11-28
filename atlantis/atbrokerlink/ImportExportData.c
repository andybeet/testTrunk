/*
 * ImportExportData.c
 *
 *  Created on: Jun 27, 2012
 *      Author: bec
 */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <sjwlib.h>
#include <netcdf.h>
#include <atlantisboxmodel.h>
#include <atUtilLib.h>
#include <LinkageInterface.h>

/**
 * Set the model time.
 *
 *
 */
int Link_SetTime(MSEBoxModel *bm, double t) {

	bm->t = t;
	bm->nt++;
	return TRUE;
}

/**
 * Set the model time.
 *
 *
 */
int Link_RunNextTimeStep(MSEBoxModel *bm) {

	/* Run the next time step */
	return runNextTimeStep(bm);
}

/**
 * Set the biomass of the given groupCode group for the given boxes to the given values.
 *
 * Will work out what proportion of existing biomass is in each layer within the boxes of interest
 * and then divide up the new biomass based on these proportions.
 *
 */
int Link_SetBiomass(MSEBoxModel *bm, char *groupCode, int numBoxes, int *boxIndexValues, double *boxValues) {

	double tolerance = 1.0;
	int cohort, guild, i, b, k, bid;
	double epi, difference;
	double *tr = Util_Alloc_Init_1D_Double(bm->wcnz, 0.0);
	double *smtr = Util_Alloc_Init_1D_Double(bm->sednz, 0.0);
	double **columnFractions = Util_Alloc_Init_2D_Double(bm->wcnz, bm->nbox, 0.0);
	double **sedColumnFractions = Util_Alloc_Init_2D_Double(bm->sednz, bm->nbox, 0.0);
	double **totalValue = Util_Alloc_Init_2D_Double(bm->K_num_max_cohort * bm->K_num_max_genetypes, bm->nbox, 0.0);

	/* Grab the group code */
	guild = Util_Get_FG_Index(bm, groupCode);
	if (guild == -1) {
		fprintf(stderr, "Linkage Link_SetBiomass - group code %s not recognised\n", groupCode);
		quit("");
	}

	fprintf(bm->logFile, "Link_SetBiomass \n");

	/* Initialise column fractions outside main logic, because
	 * they are need for all types of calculations: we may
	 * need to initialse the array containing the original
	 * fraction of biomass for that particular level in the
	 * water column, for those whose habitat is vertically
	 * distributed, to redistribute the imported value in the
	 * same proportions.  For safety, if there is no biomass
	 * left in that box (ie, a total of 0) we initialise to 0
	 * to avoid NaN or div-by-zero errors. */
	for (cohort = 0; cohort < FunctGroupArray[guild].numCohortsXnumGenes; cohort++) {
		bid = FunctGroupArray[guild].totNTracers[cohort];

		for (i = 0; i < numBoxes; i++) {
			b = boxIndexValues[i];

			totalValue[b][cohort] = 0.0;

			if (FunctGroupArray[guild].habitatType == EPIFAUNA) {
				totalValue[b][cohort] = bm->boxes[b].epi[bid];
				epi = 0;
			} else {

				if (FunctGroupArray[guild].habitatCoeffs[WC] > 0) {
					for (k = 0; k < bm->boxes[b].nz; k++)
						totalValue[b][cohort] += bm->boxes[b].tr[k][bid];
				}
				if (FunctGroupArray[guild].habitatCoeffs[SED] > 0) {
					for (k = 0; k < bm->boxes[b].sm.nz; k++)
						totalValue[b][cohort] += bm->boxes[b].sm.tr[k][bid];
				}

				/* Now work out the fractions */
				if (FunctGroupArray[guild].habitatCoeffs[WC] > 0) {
					for (k = 0; k < bm->boxes[b].nz; k++)
						columnFractions[b][k] = totalValue[b][cohort] ? (bm->boxes[b].tr[k][bid] / totalValue[b][cohort]) : 0.0;
				}

				if (FunctGroupArray[guild].habitatCoeffs[SED] > 0) {
					for (k = 0; k < bm->boxes[b].sm.nz; k++)
						sedColumnFractions[b][k] = totalValue[b][cohort] ? (bm->boxes[b].sm.tr[k][bid] / totalValue[b][cohort]) : 0.0;
				}
			}
			fprintf(bm->logFile, "totalValue[%d][%d] = %e\n", b, cohort, totalValue[b][cohort]);
		}
	}

	/* Now divide up the biomass given into each column based on columnFraction values */
	for (cohort = 0; cohort < FunctGroupArray[guild].numCohortsXnumGenes; cohort++) {
		bid = FunctGroupArray[guild].totNTracers[cohort];

		for (i = 0; i < numBoxes; i++) {
			b = boxIndexValues[i];

			Util_Init_1D_Double(tr, bm->boxes[b].nz, 0.0);
			Util_Init_1D_Double(smtr, bm->boxes[b].sm.nz, 0.0);
			if (FunctGroupArray[guild].habitatType == EPIFAUNA) {
				/* bm->boxes[b].epi[bid] += */
				epi += boxValues[i];
			} else {
				if (FunctGroupArray[guild].habitatCoeffs[WC] > 0) {
					for (k = 0; k < bm->boxes[b].nz; k++) {
						tr[k] += boxValues[i] * columnFractions[b][k];
					}
				}
				if (FunctGroupArray[guild].habitatCoeffs[SED] > 0) {
					for (k = 0; k < bm->boxes[b].sm.nz; k++) {
						smtr[k] += boxValues[i] * sedColumnFractions[b][k];
					}
				}
			}

			/* Now do the comparisons, and see if we want
			 * to update, change slightly, or quit, etc
			 */
			if (FunctGroupArray[guild].habitatType == EPIFAUNA) {
				difference = epi - bm->boxes[b].epi[bid];
				if (fabs(difference) / (bm->boxes[b].epi[bid] + small_num) > tolerance) {
					warn("box %2d:group %s cohort %d:EPIFAUNA: Imported value %lf and current value %lf are too far apart; pruning.\n", b, groupCode, cohort,
							epi, bm->boxes[b].epi[bid]);
					/* note if difference is negative this subtracts from the existing value, which is what we want: */
					bm->boxes[b].epi[bid] += tolerance * difference;
				} else
					bm->boxes[b].epi[bid] = epi;
			} else {
				for (k = 0; k < bm->boxes[b].nz; k++) {
					difference = tr[k] - bm->boxes[b].tr[k][bid];
					if (fabs(difference) / (bm->boxes[b].tr[k][bid] + small_num) > tolerance) {
						warn("box %2d:group %s cohort %d: WC, layer %d: Imported value %lf and current value %lf are too far apart; pruning.\n", b, groupCode,
								cohort, k, tr[k], bm->boxes[b].tr[k][bid]);

						/* note if difference is negative this subtracts from the existing value, which is what we want: */
						bm->boxes[b].tr[k][bid] += tolerance * difference;
					} else
						bm->boxes[b].tr[k][bid] = tr[k];
				}

				if (FunctGroupArray[guild].habitatCoeffs[SED] > 0) {
					for (k = 0; k < bm->boxes[b].sm.nz; k++) {
						difference = smtr[k] - bm->boxes[b].sm.tr[k][bid];
						if (fabs(difference) / (bm->boxes[b].sm.tr[k][bid] + small_num) > tolerance) {
							warn("box %2d:group %s cohort %d: SED layer %d: Imported value %lf and current value %lf are too far apart; pruning\n", b,
									groupCode, cohort, k, smtr[k], bm->boxes[b].sm.tr[k][bid]);
							/* note if difference is negative this subtracts from the existing value, which is what we want: */
							bm->boxes[b].sm.tr[k][bid] += tolerance * difference;
						} else
							bm->boxes[b].sm.tr[k][bid] = smtr[k];
					}
				}
			}
		}
	}
	d_free1d(tr);
	d_free1d(smtr);
	d_free2d(columnFractions);
	d_free2d(sedColumnFractions);
	d_free2d(totalValue);

	return TRUE;
}

/**
 * Get the total biomass in each specified box of interest for the given groupCode.
 *
 * Assumes that boxValues has already be allocated and will be freed once used.
 *
 */
int Link_GetBiomass(MSEBoxModel *bm, char *groupCode, int numBoxes, int *boxIndexValues, double *boxValues) {

	int guild, i, k, b;
	int cohort, bid;

	/* Grab the group code */
	guild = Util_Get_FG_Index(bm, groupCode);
	if (guild == -1) {
		fprintf(stderr, "Linkage Link_GetBiomass - group code %s not recognised\n", groupCode);
		quit("");
	}

	if (boxValues == NULL) {
		quit("Link_GetBiomass boxValues array given is NULL\n");
	}

	for (cohort = 0; cohort < FunctGroupArray[guild].numCohortsXnumGenes; cohort++) {
		bid = FunctGroupArray[guild].totNTracers[cohort];
		for (i = 0; i < numBoxes; i++) {
			b = boxIndexValues[i];

			/* Multiply by area in this case as want to have amounts in tonnes
			 and then divide by area of EwE grid cell later on to get t/km2
			 */
			if (FunctGroupArray[guild].habitatType == EPIFAUNA) {
				boxValues[i] += bm->boxes[b].epi[bid] * mg_2_tonne;
			} else {
				if (FunctGroupArray[guild].habitatCoeffs[WC] > 0) {
					for (k = 0; k < bm->boxes[b].nz; k++) {
						boxValues[i] += bm->boxes[b].tr[k][bid] * mg_2_tonne * bm->boxes[b].dz[k];
					}
				}
				if (FunctGroupArray[guild].habitatCoeffs[SED] > 0) {
					for (k = 0; k < bm->boxes[b].sm.nz; k++) {
						boxValues[i] += bm->boxes[b].sm.tr[k][bid] * mg_2_tonne * bm->boxes[b].sm.dz[k];
					}
				}
			}
		}
	}
	return TRUE;
}

/**
 * Set the biomass of the given groupCode group for the given boxes to the given values.
 *
 * Will just call the setBiomass function.
 */
int Link_SetDetritus(MSEBoxModel *bm, char *groupCode, int numBoxes, int *boxIndexValues, double *boxValues) {

	int guild, i, b, k, bid;
	double **wcColumnFractions = Util_Alloc_Init_2D_Double(bm->wcnz, bm->nbox, 0.0);
	double **sedColumnFractions = Util_Alloc_Init_2D_Double(bm->sednz, bm->nbox, 0.0);
	double *totalValue = Util_Alloc_Init_1D_Double(bm->nbox, 0.0);

	/* Grab the group code */
	guild = Util_Get_FG_Index(bm, groupCode);
	if (guild == -1) {
		fprintf(stderr, "Linkage Link_SetDetritus - group code %s not recognised\n", groupCode);
		quit("");
	}

	/* Check that the group is actually a detritus group */
	if (!FunctGroupArray[guild].isDetritus) {
		fprintf(stderr, "Linkage Link_SetDetritus - group %s not a detritus group.\n", groupCode);
		quit("");
	}

	fprintf(bm->logFile, "Link_SetBiomass \n");

	/* Initialise column fractions outside main logic, because
	 * they are need for all types of calculations: we may
	 * need to initialse the array containing the original
	 * fraction of biomass for that particular level in the
	 * water column, for those whose habitat is vertically
	 * distributed, to redistribute the imported value in the
	 * same proportions.  For safety, if there is no biomass
	 * left in that box (ie, a total of 0) we initialise to 0
	 * to avoid NaN or div-by-zero errors. */
	bid = FunctGroupArray[guild].totNTracers[0];

	for (i = 0; i < numBoxes; i++) {
		b = boxIndexValues[i];

		totalValue[b] = 0.0;

		if (FunctGroupArray[guild].habitatType != EPIFAUNA) {

			if (FunctGroupArray[guild].habitatCoeffs[WC] > 0) {
				for (k = 0; k < bm->boxes[b].nz; k++)
					totalValue[b] += bm->boxes[b].tr[k][bid];
			}
			if (FunctGroupArray[guild].habitatCoeffs[SED] > 0) {
				for (k = 0; k < bm->boxes[b].sm.nz; k++)
					totalValue[b] += bm->boxes[b].sm.tr[k][bid];
			}

			/* Now work out the fractions */
			if (FunctGroupArray[guild].habitatCoeffs[WC] > 0) {
				for (k = 0; k < bm->boxes[b].nz; k++)
					wcColumnFractions[b][k] = totalValue[b] ? (bm->boxes[b].tr[k][bid] / totalValue[b]) : 0.0;
			}

			if (FunctGroupArray[guild].habitatCoeffs[SED] > 0) {
				for (k = 0; k < bm->boxes[b].sm.nz; k++)
					sedColumnFractions[b][k] = totalValue[b] ? (bm->boxes[b].sm.tr[k][bid] / totalValue[b]) : 0.0;
			}
		}
		fprintf(bm->logFile, "totalValue[%d] = %e\n", b, totalValue[b]);
	}

	for (i = 0; i < numBoxes; i++) {
		b = boxIndexValues[i];

		if (FunctGroupArray[guild].habitatType != EPIFAUNA) {
			if (FunctGroupArray[guild].habitatCoeffs[WC] > 0) {
				for (k = 0; k < bm->boxes[b].nz; k++) {
					bm->linkageInterface->linkageWCDetritusFlux[b][k][guild] = boxValues[b] / 100 * wcColumnFractions[b][k] / bm->boxes[b].dz[k]; /* Bec change 18th Oct */

				}
			}
			if (FunctGroupArray[guild].habitatCoeffs[SED] > 0) {
				for (k = 0; k < bm->boxes[b].sm.nz; k++) {
					bm->linkageInterface->linkageSEDDetritusFlux[b][k][guild] = boxValues[b] / 100 * sedColumnFractions[b][k] / bm->boxes[b].sm.dz[k]; /* Bec change 18th Oct */
				}
			}
		}
		/* Shouldn't be any in the EPI layer, but
		 * we can assume the functional groups are
		 * set up correctly I think */

	}
	for (b = 0; b < bm->nbox; b++) {
		for (k = 0; k < bm->boxes[b].nz; k++) {
			if (bm->linkageInterface->linkageWCDetritusFlux[b][k][guild] > 0)
				fprintf(bm->logFile, "bm->atEcologyModule->linkageWCDetritusFlux[%d][%d][%s] = %e\n", b, k, FunctGroupArray[guild].groupCode,
						bm->linkageInterface->linkageWCDetritusFlux[b][k][guild]);
		}
	}

	d_free2d(wcColumnFractions);
	d_free2d(sedColumnFractions);
	d_free1d(totalValue);
	return TRUE;

}
/**
 * Get the total biomass in each specified box of interest for the given groupCode.
 *
 * Will just call the getBiomass function.
 */
int Link_GetDetritus(MSEBoxModel *bm, char *groupCode, int numBoxes, int *boxIndexValues, double *boxValues) {
	return Link_GetBiomass(bm, groupCode, numBoxes, boxIndexValues, boxValues);

}

int Link_SetMortality(MSEBoxModel *bm, char *groupCode, int numBoxes, int *boxIndexValues, double *boxValues) {
	int b, i, cohort, guild, bid, k;

	/*
	 http://wiki.csiro.au/confluence/display/Atlantis/Interaction+with+atecology

	 These values are stored in the following new arrays:


	 double ****linkageWCMortality;
	 double ****linkageSEDMortality;
	 double ***linkageEPIMortality;

	 These are allocated as:

	 linkageEPIMortality = AllocInit3DArrayDouble(bm->K_num_cohorts, bm->K_num_tot_sp, bm->nbox, 0.0);
	 linkageWCMortality = AllocInit4DArrayDouble(bm->K_num_cohorts, bm->K_num_tot_sp, bm->wcnz, bm->nbox, 0.0);
	 linkageSEDMortality = AllocInit4DArrayDouble(bm->K_num_cohorts, bm->K_num_tot_sp, bm->sednz, bm->nbox, 0.0);

	 The appropriate mortality value is then used in each
	 rate of change (ROC) equation for all functional
	 groups. These values are subtracted similar to how the
	 mortality values are indicating a reduction in
	 biomass.

	 The units of the mortality values should be mgN /m^3
	 for water column and sediment and mgN/m^2 for the
	 epibenthic layer (2D).

	 So the standard ROC equation is:
	 flux = growth - naturalmortality - mortalityDueToPredator;

	 So we have now included an additional term in these equations:
	 flux = growth - naturalmortality - mortalityDueToPredator - linkageMortality.
	 */

	/* Grab the group code */
	guild = Util_Get_FG_Index(bm, groupCode);
	if (guild == -1) {
		fprintf(stderr, "Linkage Link_GetBiomass - group code %s not recognised\n", groupCode);
		quit("");
	}

	for (cohort = 0; cohort < FunctGroupArray[guild].numCohortsXnumGenes; cohort++) {
		bid = FunctGroupArray[guild].totNTracers[cohort];

		for (i = 0; i < numBoxes; i++) {
			b = boxIndexValues[i];

			if (FunctGroupArray[guild].habitatCoeffs[EPIFAUNA] > 0) {
				bm->linkageInterface->linkageEPIMortality[b][guild][cohort] = (boxValues[i] / 100) * bm->boxes[b].epi[bid];
			} else {
				if (FunctGroupArray[guild].habitatCoeffs[WC] > 0) {
					for (k = 0; k < bm->boxes[b].nz; k++) {
						bm->linkageInterface->linkageWCMortality[b][k][guild][cohort] = (boxValues[i] / 100) * bm->boxes[b].tr[k][bid];
					}
				}
				if (FunctGroupArray[guild].habitatCoeffs[SED] > 0) {
					for (k = 0; k < bm->boxes[b].sm.nz; k++) {
						bm->linkageInterface->linkageSEDMortality[b][k][guild][cohort] = (boxValues[i] / 100) * bm->boxes[b].sm.tr[k][bid];
					}
				}
			}
		}

	}
	return TRUE;
}

int Link_GetMortality(MSEBoxModel *bm, char *groupCode, int numBoxes, int *boxIndexValues, double *boxValues) {
	quit("Linkage - Link_GetMortality - Function not yet supported\n");

	return TRUE;
}

