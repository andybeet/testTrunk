/**
\file
\brief The conversion code to take Atlantis variables to pass to EwE
\ingroup atLink

    File:		atlinkconversion.c
	Created:	4/9/2008
	Author:		Beth Fulton,   CSIRO Marine Laboratories,  2008
	Purpose:	C file for converting atlantis variables to format and group
				structure for passing to EwE model via link server message
    Note:

    Revisions:	02-02-2009 Beth Fulton
    			Changed the variable called small to small_num.

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

/************************************************************************************
Routine to handle the data transfer to and from EwE
*/

void LinkToModel(MSEBoxModel *bm)
{
	int l=0;

	/* Initialise the array if necessary */
	if(	bm->eweLinkModule->convertPRMS.refresh == TRUE){
		for(l=0; l<bm->eweLinkModule->convertPRMS.ingroup; l++){
			Util_Init_2D_Double(bm->eweLinkModule->ExportArray[l].values, eweLinkModule->convertPRMS.io_nrow, bm->eweLinkModule->convertPRMS.io_ncol, 0.0);
		}
	}

	/* Create the array */
	bm->eweLinkModule->convertPRMS.cumday += bm->dt;

	switch(bm->eweLinkModule->convertPRMS.timehandling){
		case WEIGHTED_AVERAGE:
			WeightedAvg(bm);
			break;
		case FIXED_POINT:
			quit("CreateConvertedArray: No code for snapshot form of conversion as yet\n");
			break;
		default:
			quit("CreateConvertedArray: No code for this form of conversion as yet\n");
	}

	/* Send data to EwE */
	if(bm->eweLinkModule->convertPRMS.cumday > bm->eweLinkModule->convertPRMS.EwEdt){
		if(SendArrays(bm) == FALSE)
			quit("Array send through link failed\n");

	/* Call back results from EwE */
		if(GetArrays(bm) == FALSE)
			quit("Array recieve through link failed\n");

		bm->eweLinkModule->convertPRMS.refresh = FALSE;
	}

	return;
}


/************************************************************************************
Routine to convert Atlantis state variables to the array to send to EwE
*/

void CreateConvertedArray(MSEBoxModel *bm)
{
	/* Create the array */
	bm->eweLinkModule->convertPRMS.cumday += bm->dt;

	switch(bm->eweLinkModule->convertPRMS.timehandling){
		case WEIGHTED_AVERAGE:
			WeightedAvg(bm);
			break;
		case FIXED_POINT:
			quit("CreateConvertedArray: No code for snapshot form of conversion as yet\n");
			break;
		default:
			quit("CreateConvertedArray: No code for this form of conversion as yet\n");
	}

	return;
}

/************************************************************************************
Routine to convert Atlantis state variables to the array to send to EwE
*/

void WeightedAvg(MSEBoxModel *bm)
{
	int layer=0, fgIndex=0, b=0, nr=0, nc=0, col=0, row=0, k=0, l;
	double area=0.0;
	int cohort, bid, pid, nid;
	int c, r;
	double dz = 1.0; //REMOVE

	/* Update the array */
	switch(bm->eweLinkModule->convertPRMS.spacehandling){
		case WEIGHTED_AVERAGE:
			for(fgIndex = 0; fgIndex < bm->K_num_tot_sp; fgIndex++){

			/* Note only running for biomass pools for now.
			   More consultation needed with Villy for age structured groups as will
			   need to know about stanzas etc in EwE
			*/
				if(FunctGroupArray[fgIndex].isExported == TRUE){
					if(FunctGroupArray[fgIndex].groupAgeType == BIOMASS){
						switch(FunctGroupArray[fgIndex].exportType){
							case TYPE_BIOMASS:  /* Create biomass array */
								for(cohort=0; cohort < FunctGroupArray[fgIndex].numCohortsXnumGenes; cohort++){

									bid = FunctGroupArray[fgIndex].totNTracers[cohort];
									for(b=0; b<bm->nbox; b++){
										for(nr=0; nr<bm->boxes[b].export->numrows; nr++){
											for(nc=0; nc<bm->boxes[b].export->numcols; nc++){
												row = bm->boxes[b].export->rows[nr];
												col	= bm->boxes[b].export->cols[nc];
												area = bm->boxes[b].export->bcellarea[nr][nc];  // area of box that contributes to cell [nr][nc] - calculated on startup
												/* Multiply by area in this case as want to have amounts in tonnes
												   and then divide by area of EwE grid cell later on to get t/km2
											     */
												if(FunctGroupArray[fgIndex].habitatType == EPIFAUNA){
													bm->eweLinkModule->ExportArray[layer].values[row][col] += bm->boxes[b].epi[bid] * area * mg_2_tonne;
													bm->eweLinkModule->ExportArray[layer].dataType = TYPE_BIOMASS;
												} else {
													if(FunctGroupArray[fgIndex].habitatType == WC){
														for(k=0; k<bm->boxes[b].nz; k++){
															bm->eweLinkModule->ExportArray[layer].values[row][col] += bm->boxes[b].tr[k][bid] * area * mg_2_tonne * bm->boxes[b].dz[k];
															bm->eweLinkModule->ExportArray[layer].dataType = TYPE_BIOMASS;
														}
													}
													if(FunctGroupArray[fgIndex].habitatType == SED){
														for(k=0; k<bm->boxes[b].sm.nz; k++){
															bm->eweLinkModule->ExportArray[layer].values[row][col] += bm->boxes[b].sm.tr[k][bid] * area * mg_2_tonne * bm->boxes[b].sm.dz[k];
															bm->eweLinkModule->ExportArray[layer].dataType = TYPE_BIOMASS;
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
								for(cohort=0; cohort < FunctGroupArray[fgIndex].numCohortsXnumGenes; cohort++){
									bid = FunctGroupArray[fgIndex].totNTracers[cohort];
									pid = FunctGroupArray[fgIndex].prodnTracers[cohort];

									for(b=0; b<bm->nbox; b++){
										for(nr=0; nr<bm->boxes[b].export->numrows; nr++){
											for(nc=0; nc<bm->boxes[b].export->numcols; nc++){
												row = bm->boxes[b].export->rows[nr];
												col	= bm->boxes[b].export->cols[nc];
												area = bm->boxes[b].export->bcellarea[nr][nc];  // area of box that contributes to cell [nr][nc] - calculated on startup
												/* Multiply by area in this case so get weighted contribution
												   from Atlantis cell to EwE cell based on relative overlap
												   of the cells
											     */
												if(FunctGroupArray[fgIndex].habitatType == EPIFAUNA){
													bm->eweLinkModule->ExportArray[layer].values[row][col] += area * s_2_yr * (bm->boxes[b].epi[pid] / (bm->boxes[b].epi[bid] + small_num));
													bm->eweLinkModule->ExportArray[layer].dataType = TYPE_PRODUCTION;
												} else {
													if(FunctGroupArray[fgIndex].habitatType == WC){
														for(k=0; k<bm->boxes[b].nz; k++){
															bm->eweLinkModule->ExportArray[layer].values[row][col] += area * s_2_yr * (bm->boxes[b].tr[k][pid] / (bm->boxes[b].tr[k][bid] + small_num));
															bm->eweLinkModule->ExportArray[layer].dataType = TYPE_PRODUCTION;
														}
													}
													if(FunctGroupArray[fgIndex].habitatType == SED){
														for(k=0; k<bm->boxes[b].sm.nz; k++){
															bm->eweLinkModule->ExportArray[layer].values[row][col] += area * s_2_yr * (bm->boxes[b].sm.tr[k][pid] / (bm->boxes[b].tr[k][bid] + small_num));
															bm->eweLinkModule->ExportArray[layer].dataType = TYPE_PRODUCTION;
														}
													}
												}
											}
										}
									}
								}
								layer++;
								break;
							case TYPE_MORTALITY:  /* Create mortality array */
								quit("WeightedAvg: This kind of data packing is not defined yet\n");
								break;
							default:
								quit("WeightedAvg: This kind of data packing is not defined yet\n");
								break;
						}
					} else {
						quit("WeightedAvg: This code doesn't exist yet\n");
					}
				}
			}

			for(fgIndex = 0; fgIndex < bm->K_num_physiochem; fgIndex++){
			/* For physical and chemical properties	*/
				if(PhysioChemArray[fgIndex].isExported == TRUE){
					nid = *PhysioChemArray[fgIndex].tracerIndex; /* Get the tracer index pointed to by the totTracers pointer */
					switch(PhysioChemArray[fgIndex].exportType){
						case TYPE_NUTRIENTS:  /* Create nutrient array */
							for(b=0; b<bm->nbox; b++){
								for(nr=0; nr<bm->boxes[b].export->numrows; nr++){
									for(nc=0; nc<bm->boxes[b].export->numcols; nc++){
										row = bm->boxes[b].export->rows[nr];
										col	= bm->boxes[b].export->cols[nc];
										area = bm->boxes[b].export->bcellarea[nr][nc];  // area of box that contributes to cell [nr][nc] - calculated on startup
										if(PhysioChemArray[fgIndex].habitat_type == WC){
											for(k=0; k<bm->boxes[b].nz; k++){
												bm->eweLinkModule->ExportArray[layer].values[row][col] += bm->boxes[b].tr[k][bid] * area * mg_2_tonne * dz;
												bm->eweLinkModule->ExportArray[layer].dataType = TYPE_NUTRIENTS;
											}
										}
										if(PhysioChemArray[fgIndex].habitat_type == SED){
											for(k=0; k<bm->boxes[b].sm.nz; k++){
												bm->eweLinkModule->ExportArray[layer].values[row][col] += bm->boxes[b].sm.tr[k][bid] * area * mg_2_tonne * dz;
												bm->eweLinkModule->ExportArray[layer].dataType = TYPE_NUTRIENTS;
											}
										}
									}
								}
							}
							layer++;
							break;
						case TYPE_PHYSCHEM: /* Create physical properties array */
							for(b=0; b<bm->nbox; b++){
								k = bm->boxes[b].nz-1;
								for(nr=0; nr<bm->boxes[b].export->numrows; nr++){
									for(nc=0; nc<bm->boxes[b].export->numcols; nc++){
										row = bm->boxes[b].export->rows[nr];
										col	= bm->boxes[b].export->cols[nc];
										area = bm->boxes[b].export->bcellarea[nr][nc];  // area of box that contributes to cell [nr][nc] - calculated on startup

										bm->eweLinkModule->ExportArray[layer].values[row][col] += bm->boxes[b].tr[k][bid] * area;
										bm->eweLinkModule->ExportArray[layer].dataType = TYPE_PHYSCHEM;
									}
								}
							}
							layer++;
							break;
						default:
							quit("WeightedAvg: This kind of physiochem data packing is not defined yet\n");
							break;
					}
				}
			}


			/* Convert biomass, production, mortality and the like from tonnes in cell to t/km2 in cell */
			for(l=0; l<bm->eweLinkModule->convertPRMS.ingroup; l++){
				for(r=0; r<bm->eweLinkModule->convertPRMS.io_nrow; r++){
					for(c=0; c<bm->eweLinkModule->convertPRMS.io_ncol; c++){
						bm->eweLinkModule->ExportArray[l].values[r][c] /= bm->eweLinkModule->ExportArray[l].area[r][c]; // area created on readin
					}
				}
			}
			break;
		case FIXED_POINT:
			quit("CreateConvertedArray: No code for centre_point conversion as yet\n");
			break;
		default:
			quit("CreateConvertedArray: No code for this form of conversion as yet\n");

	}

	return;
}
