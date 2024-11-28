/**
\file
\brief The code to get variables from EwE to pass to Atlantis
\ingroup atLink

    File:		atlinkimport.c
	Created:	4/9/2008
	Author:		Beth Fulton,   CSIRO Marine Laboratories,  2008
	Purpose:	C file for pulling EwE variables to format into
				Atlantis via link server message
    Note:

    Revisions:

*************/
#include <stdio.h>
#include <stdlib.h>
#include <signal.h>
#include <math.h>
#include <ctype.h>
#include <time.h>
#include <string.h>
#include <sjwlib.h>
#include <netcdf.h>
#include <atlantisboxmodel.h>
#include <atUtilLib.h>
#include <atlink.h>
#include <atComms.h>


/************************************************************************************
Routine to unpack EwE array to Atlantis state variables
*/

void UnpackArrays(MSEBoxModel *bm)
{
	int fgIndex, layer=0, b=0, nr=0, nc=0, row=0, col=0, k=0;
	double area=0, totval=0, origval;
	int cohort, bid;

	/* Update the array */
	switch(bm->eweLinkModule->convertPRMS.spacehandling){
		case WEIGHTED_AVERAGE:
			for(fgIndex = 0; fgIndex < bm->K_num_tot_sp; fgIndex++){
			/* Note only running for biomass pools for now.
			   More consultation needed with Villy for age structured groups as will
			   need to know about stanzas etc in EwE
			*/
				if(FunctGroupArray[fgIndex].isImported == TRUE){
					if(FunctGroupArray[fgIndex].groupAgeType == BIOMASS){
						switch(bm->eweLinkModule->ImportArray[layer].dataType){
							case TYPE_BIOMASS:  /* Create biomass array */
								for(cohort=0; cohort < FunctGroupArray[fgIndex].numCohortsXnumGenes; cohort++){
									bid = FunctGroupArray[fgIndex].totNTracers[cohort];
									for(b=0; b<bm->nbox; b++){
										for(nr=0; nr<bm->boxes[b].export->numrows; nr++){
											for(nc=0; nc<bm->boxes[b].export->numcols; nc++){
												row = bm->boxes[b].export->rows[nr];
												col	= bm->boxes[b].export->cols[nc];
												area = bm->boxes[b].export->bcellarea[nr][nc];  // area of box that contributes to cell [nr][nc] - calculated on startup
												if(FunctGroupArray[fgIndex].habitatType == EPIFAUNA){
													bm->boxes[b].epi[bid] = bm->eweLinkModule->ImportArray[layer].values[row][col] * tonne_2_mg * (area /  bm->eweLinkModule->ExportArray[layer].area[nr][nc]) * per_m2_from_per_km2;
												} else {
													if(FunctGroupArray[fgIndex].habitatType == WC){
														totval = 0;
														for(k=0; k<bm->boxes[b].nz; k++){
															totval += bm->boxes[b].tr[k][bid] * bm->boxes[b].dz[k];
														}
														for(k=0; k<bm->boxes[b].nz; k++){
															origval = bm->boxes[b].tr[k][bid] * bm->boxes[b].dz[k];
															bm->boxes[b].tr[k][bid] = bm->eweLinkModule->ImportArray[layer].values[row][col] * tonne_2_mg * (area / bm->eweLinkModule->ExportArray[layer].area[nr][nc]) * per_m2_from_per_km2 * (origval / totval) / bm->boxes[b].dz[k];
														}
													}
													if(FunctGroupArray[fgIndex].habitatType == SED){
														totval = 0;
														for(k=0; k<bm->boxes[b].sm.nz; k++){
															totval += bm->boxes[b].sm.tr[k][bid] * bm->boxes[b].sm.dz[k];
														}
														for(k=0; k<bm->boxes[b].sm.nz; k++){
															origval = bm->boxes[b].sm.tr[k][bid] * bm->boxes[b].sm.dz[k];
															bm->boxes[b].sm.tr[k][bid] = bm->eweLinkModule->ImportArray[layer].values[row][col] * tonne_2_mg * (area / bm->eweLinkModule->ExportArray[layer].area[nr][nc]) * per_m2_from_per_km2 * (origval / totval) / bm->boxes[b].sm.dz[k];
														}
													}
												}
											}
										}
									}
								}
								layer++;
								break;
							case TYPE_PRODUCTION:  /* Create production array - actually P/B */
								quit("UnpackArrays: This kind of data unpacking is not defined yet\n");
								break;
							case TYPE_MORTALITY:  /* Create mortality array */
								quit("UnpackArrays: This kind of data unpacking is not defined yet\n");
								break;
							default:
								quit("UnpackArrays: This kind of data unpacking is not defined yet\n");
								break;
						}
					} else {
						quit("UnpackArrays: This code doesn't exist yet\n");
					}
				}
			}

			/* Assuming not importing nutrient or physical properties from EwE as yet - may need to rethink later on */

			break;
		case FIXED_POINT:
			quit("UnpackArrays: No code for centre_point unpack as yet\n");
			break;
		default:
			quit("UnpackArrays: No code for this form of unpack as yet\n");

	}

	return;
}


/************************************************************************************
Routine where Atlantis calls EwE and gets state variables from EwE
*/

int GetArrays(MSEBoxModel *bm)
{
	int ans = 0, l=0;

	/* Get data from EwE */
	for(l=0; l<bm->eweLinkModule->convertPRMS.ingroup; l++){

		/* Call to EwE */
		switch(bm->eweLinkModule->ImportArray[l].dataType){
			case TYPE_BIOMASS:
				ans += getBiomass(bm->eweLinkModule->ImportArray[l].parameterName, bm->eweLinkModule->ImportArray[l].values, bm->eweLinkModule->convertPRMS.io_ncol, bm->eweLinkModule->convertPRMS.io_nrow);
				break;
			case TYPE_PRODUCTION:
				ans += getProduction(bm->eweLinkModule->ImportArray[l].parameterName, bm->eweLinkModule->ImportArray[l].values, bm->eweLinkModule->convertPRMS.io_ncol, bm->eweLinkModule->convertPRMS.io_nrow);
                break;
			case TYPE_MORTALITY:
				ans += getMortality(bm->eweLinkModule->ImportArray[l].parameterName, bm->eweLinkModule->ImportArray[l].values, bm->eweLinkModule->convertPRMS.io_ncol, bm->eweLinkModule->convertPRMS.io_nrow);
                break;
			case TYPE_NUTRIENTS:
				ans += getNutrients(bm->eweLinkModule->ImportArray[l].parameterName, bm->eweLinkModule->ImportArray[l].values, bm->eweLinkModule->convertPRMS.io_ncol, bm->eweLinkModule->convertPRMS.io_nrow);
                break;
			case TYPE_PHYSCHEM:
				ans += getProperties(bm->eweLinkModule->ImportArray[l].parameterName, bm->eweLinkModule->ImportArray[l].values, bm->eweLinkModule->convertPRMS.io_ncol, bm->eweLinkModule->convertPRMS.io_nrow);
                break;
			default:
				quit("GetArrays: This data type does not exist\n");
				break;
		}
	}

	UnpackArrays(bm);

	if(ans < bm->eweLinkModule->convertPRMS.ingroup)
		return 0;
	else
		return 1;
}
