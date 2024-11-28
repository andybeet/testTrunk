/**************************************************************************//**
\file
\brief The code to interface with the CLAM model.
\ingroup atCLAMLink

    File:		atCLAMLink.c
	Created:	4/9/2008
	Author:		Bec Gorton,   CSIRO Marine Laboratories,  2009
	Purpose:	The code to link to the CLAM model.
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
#include <atlantisboxmodel.h>
#include <atUtilLib.h>
#include "atCLAMLinkPrivate.h"

static void Calculate_Max_Output_States(MSEBoxModel *bm);

/**
 * \brief Initialise the CLAM-Atlantis linkage module.
 *
 * This function will read in the linkage parameters from the bm->CLAMLinkageprmIfname file.
 * The necessary arrays will be allocated, the CLAM ini file will be generated (this only needs to happen at the start
 * of the model run) and the maximum number of output states will be calculated.
 *
 */
int CLAM_Init(MSEBoxModel *bm){

	int i, index;

	for(i = 0; i < numClamIndicators; i++){
		for(index= 0; index < CLAMOutputOptions[i].numStates; index++){
			fprintf(bm->logFile, "CLAMOutputOptions[%d].OutputStates[index] = %e\n", i, CLAMOutputOptions[i].OutputStates[index]);
		}
	}


    Read_CLAM_Link_Params(bm, bm->CLAMLinkageprmIfname);

	Allocate_CLAM_Arrays(bm);

    /* Create the CLAM ini file */
    //Create_CLAM_Ini_File(bm);

    /*Find the maximum number of output states */
    Calculate_Max_Output_States(bm);

    return TRUE;
}


/**
 * \brief Allocate the arrays used in the code that determines the 'best' management scenario to implement
 * within Atlantis.
 */
void Allocate_CLAM_Arrays(MSEBoxModel *bm){

    outcomeValueArray = Util_Alloc_Init_4D_Double(MAXNUMOUTPUTSTATES, maxStrategyOptions, numScenarioOptions, numClamIndicators, 0.0);
    expectedValue = Util_Alloc_Init_3D_Double(numClamIndicators, maxStrategyOptions, numScenarioOptions, 0.0);
    variance = Util_Alloc_Init_3D_Double(numClamIndicators, maxStrategyOptions, numScenarioOptions, 0.0);
    distance = Util_Alloc_Init_3D_Double(numClamIndicators, maxStrategyOptions, numScenarioOptions, 0.0);
    euclidDistance = Util_Alloc_Init_2D_Double(maxStrategyOptions, numScenarioOptions, 0.0);
    indicatorWeightings = Util_Alloc_Init_1D_Double(numClamIndicators, 0.0);
}
/**
 * \brief Free the arrays used in the code that determines the 'best' management scenario to implement
 * within Atlantis.
 */
void Free_CLAM_Arrays(MSEBoxModel *bm){

    free4d(outcomeValueArray);
    free3d(expectedValue);
    free3d(variance);
    free3d(distance);
    free2d(euclidDistance);
    free1d(indicatorWeightings);
    free(indicators);
    free(CLAMScenarioOptions);
}

void CLAM_Free(MSEBoxModel *bm){
	Free_CLAM_Arrays(bm);
}
/*
 * \brief Find the max number of output states from the CLAMOutputOptions array.
 */
static void Calculate_Max_Output_States(MSEBoxModel *bm){

    int i;
    int count = 0;
    int maxCount = 0;

    /* Find the max number of output states. */
    for(i = 0; i < numClamIndicators; i++){
        count = CLAMOutputOptions[i].numStates;
        maxCount = max(maxCount, count);
    }
    maxNumOutputStates = maxCount;
}
