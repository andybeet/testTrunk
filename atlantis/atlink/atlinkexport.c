/**
\file
\brief The code for Atlantis to call to EwE - and do iniitalisation of params for that
\ingroup atLink

    File:		atlinkexport.c
	Created:	4/9/2008
	Author:		Beth Fulton,   CSIRO Marine Laboratories,  2008
	Purpose:	C file for sending atlantis variables to EwE
				model via link server message
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
#include "atComms.h"


/**
 * \brief Initialised the Link Module
 *
 */
void initialiseLinkModule(MSEBoxModel *bm)
{
	int group;
	bm->eweLinkModule = (EwELinkStruct *)malloc(sizeof(EwELinkStruct));

	/* The transport arrays */

	/* Set up Export and Import arrays */
	//ExportArray = Util_Alloc_Init_4D_Double(bm->eweLinkModule->convertPRMS.io_ncol,bm->eweLinkModule->convertPRMS.io_nrow,2, bm->eweLinkModule->convertPRMS.ingroup, 0.0);
	//ImportArray = Util_Alloc_Init_4D_Double(bm->eweLinkModule->convertPRMS.io_ncol,bm->eweLinkModule->convertPRMS.io_nrow,2, bm->eweLinkModule->convertPRMS.outgroup, 0.0);

	bm->eweLinkModule->ImportArray = (LinkTransportStruct *)malloc(sizeof(LinkTransportStruct) * bm->eweLinkModule->convertPRMS.ingroup);
	bm->eweLinkModule->ExportArray = (LinkTransportStruct *)malloc(sizeof(LinkTransportStruct) * bm->eweLinkModule->convertPRMS.outgroup);

	for(group=0; group<bm->eweLinkModule->convertPRMS.ingroup; group++){
		bm->eweLinkModule->ImportArray[group].values = Util_Alloc_Init_2D_Double(bm->eweLinkModule->convertPRMS.io_ncol,bm->eweLinkModule->convertPRMS.io_nrow, 0.0);
		bm->eweLinkModule->ExportArray[group].values = Util_Alloc_Init_2D_Double(bm->eweLinkModule->convertPRMS.io_ncol,bm->eweLinkModule->convertPRMS.io_nrow, 0.0);

		bm->eweLinkModule->ImportArray[group].area = Util_Alloc_Init_2D_Double(bm->eweLinkModule->convertPRMS.io_ncol,bm->eweLinkModule->convertPRMS.io_nrow, 0.0);
		bm->eweLinkModule->ExportArray[group].area = Util_Alloc_Init_2D_Double(bm->eweLinkModule->convertPRMS.io_ncol,bm->eweLinkModule->convertPRMS.io_nrow, 0.0);
	}
}

/**
 * \brief Free the Link Module
 *
 */
void freeLinkModule(MSEBoxModel *bm)
{
	int group;

	for(group=0; group<bm->eweLinkModule->convertPRMS.ingroup; group++){
		free(bm->eweLinkModule->ImportArray[group].values);
		free(bm->eweLinkModule->ExportArray[group].values);

		free(bm->eweLinkModule->ImportArray[group].area);
		free(bm->eweLinkModule->ExportArray[group].area);

	}
	free(bm->eweLinkModule->ImportArray);
	free(bm->eweLinkModule->ExportArray);

	free(eweLinkModule);
}

/**
 * \brief Load the Atlantis to EwE conversion values.
 *
 */
void LoadCoversionParameters(MSEBoxModel *bm)
{
	int b=0, r=0, c=0, nr=0, nc=0;
	int l;
	double area;
	double xxx = 1.0, yyy = 1.0; //REMOVE

	initialiseLinkModule(bm);


	/* Do loady stuff here */

	/* Hard wired code here */
	bm->eweLinkModule->convertPRMS.refresh = TRUE;

	bm->eweLinkModule->convertPRMS.numdt = 0.0;
	bm->eweLinkModule->convertPRMS.ingroup = 1;
	bm->eweLinkModule->convertPRMS.outgroup = 1;
	bm->eweLinkModule->convertPRMS.io_nrow = 72;	// Get rid of this hardwiring for Cathy's ASAP
	bm->eweLinkModule->convertPRMS.io_ncol = 60;    // Get rid of this hardwiring for Cathy's ASAP
	bm->eweLinkModule->convertPRMS.cumday = 0.0;
	bm->eweLinkModule->convertPRMS.EwEdt = 30.0;    // How deal with fact not always 30 days a month??????
	bm->eweLinkModule->convertPRMS.timehandling = 0;
	bm->eweLinkModule->convertPRMS.spacehandling = 0;
	bm->eweLinkModule->convertPRMS.ulx = 1.0;
	bm->eweLinkModule->convertPRMS.uly = 1.0;
	bm->eweLinkModule->convertPRMS.ullat = 1.0;
	bm->eweLinkModule->convertPRMS.ullong = 1.0;
	bm->eweLinkModule->convertPRMS.ulref = 0.0;		// Spatial referencing method (0=(x,y) or 1=(lat,long))
	bm->eweLinkModule->convertPRMS.io_area = 5000;

	/* Get number of cells per box being overlapped */
	for(b=0; b<bm->nbox; b++){
		/* Get polygon minx, maxx, miny, maxy */


		/* Get number of rows overlapped - TO DO */
		nr = 4;

		/* Get number of cols overlapped - TO DO */
        nc = 4;

		/* Create memory for the export structure */
		bm->boxes[b].export = (ExportInfo *)malloc(sizeof(ExportInfo));

		/* Create memory allocation for cell ids */
		bm->boxes[b].export->rows = (int *)i_alloc1d(nr);
		bm->boxes[b].export->cols = (int *)i_alloc1d(nc);

		/* Create memory allocation for areas */
		bm->boxes[b].export->bcellarea = (double **)alloc2d(nc,nr);

		/** Now load in data **/

		/* Readin numrows and numcols */
		bm->boxes[b].export->numrows = nr;
		bm->boxes[b].export->numrows = nc;

		/* Readin in cell ids */
		for(r=0; r<nr; r++){
			bm->boxes[b].export->rows[r] = (int)xxx;
		}
		for(c=0; c<nc; c++){
			bm->boxes[b].export->cols[c] = (int)yyy;
		}

		/* Readin cell areas */
		area = bm->eweLinkModule->convertPRMS.io_area;  // For now just assume 100% overlap of cells - TODO NEEDS TO BE FIXED
		for(r=0; r<nr; r++){
			for(c=0; c<nc; c++){
				bm->boxes[b].export->bcellarea[nr][nc] = area;
			}
		}
	}

	/* Set up Export and Import arrays */

	//ExportName = define array here and then load in the names
	/* The ararys are already initialised when they are allocated in initialiseLinkModule*/
	for(l=0; l<bm->eweLinkModule->convertPRMS.ingroup; l++){
		bm->eweLinkModule->ImportArray[l].dataType = TYPE_BIOMASS;
		bm->eweLinkModule->ExportArray[l].dataType = -1;


//		for(r=0; r<bm->eweLinkModule->convertPRMS.io_nrow; r++){
//			for(c=0; c<bm->eweLinkModule->convertPRMS.io_ncol; c++){
//				ExportArray[l][value_id][r][c] = 0.0;
//				ExportArray[l][type_id][r][c] = -1;
//
//				ImportArray[l][value_id][r][c] = 0.0;
//
//				// For now all import array types are biomass, but that needs to be generalised
//				ImportArray[l][type_id][0][0] = TYPE_BIOMASS;
//			}
//		}
	}

	return;
}

/************************************************************************************
Routine where Atlantis calls EwE and sends state variables to EwE
*/

int SendArrays(MSEBoxModel *bm)
{
	int ans = 0, l=0;

	/* Intialise the array - if need be */
	for(l=0; l<bm->eweLinkModule->convertPRMS.ingroup; l++){

//		bm->eweLinkModule->ExportArray[l].dataType = -1;
		// Not sure why you would want to do this.
//		for(r=0; r<bm->eweLinkModule->convertPRMS.io_nrow; r++){
//			for(c=0; c<bm->eweLinkModule->convertPRMS.io_ncol; c++){
//				bm->eweLinkModule->ExportArray[l].values[r][c] = 0.0;
//				ExportArray[l][value_id][r][c] = 0.0;
//				ExportArray[l][type_id][r][c] = -1;
//			}
//		}

		/* Call to EwE */
		switch(bm->eweLinkModule->ExportArray[l].dataType){
			case TYPE_BIOMASS:
				setBiomass(bm->eweLinkModule->ExportArray[l].parameterName, bm->eweLinkModule->ExportArray[l].values,
				bm->eweLinkModule->convertPRMS.io_ncol, bm->eweLinkModule->convertPRMS.io_nrow);
                ans++;
                break;
			case TYPE_PRODUCTION:
				setProduction(bm->eweLinkModule->ExportArray[l].parameterName, bm->eweLinkModule->ExportArray[l].values, bm->eweLinkModule->convertPRMS.io_ncol, bm->eweLinkModule->convertPRMS.io_nrow);
                ans++;
                break;
			case TYPE_MORTALITY:
				setMortality(bm->eweLinkModule->ExportArray[l].parameterName, bm->eweLinkModule->ExportArray[l].values, bm->eweLinkModule->convertPRMS.io_ncol, bm->eweLinkModule->convertPRMS.io_nrow);
                ans++;
                break;
			case TYPE_NUTRIENTS:
				setNutrients(bm->eweLinkModule->ExportArray[l].parameterName, bm->eweLinkModule->ExportArray[l].values, bm->eweLinkModule->convertPRMS.io_ncol, bm->eweLinkModule->convertPRMS.io_nrow);
                ans++;
                break;
			case TYPE_PHYSCHEM:
				setProperties(bm->eweLinkModule->ExportArray[l].parameterName, bm->eweLinkModule->ExportArray[l].values, bm->eweLinkModule->convertPRMS.io_ncol, bm->eweLinkModule->convertPRMS.io_nrow);
                ans++;
                break;
			default:
				quit("WeightedAvg: This data type does not exist\n");
				break;
		}

	}

	/* Identify that array reinitialisation is needed */
	if(ans < bm->eweLinkModule->convertPRMS.ingroup){
		return 0;
	}

	bm->eweLinkModule->convertPRMS.refresh = TRUE;
	return 1;
}
