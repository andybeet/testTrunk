/**
 * atCLAMImplement.c
 *
 *	\brief This file contains the functions used to implement CLAM strategies in Atlantis.
 *
 *  Created on: 15/04/2010
 *      Author: bec
 */
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
#include <atEcologyLib.h>
#include "atHarvestLib.h"
#include <atCLAMLinkPrivate.h>

static void Scale_PointSource(MSEBoxModel *bm, char *name, int location, double start, double period, double mult, FILE *llogfp);
//void Scale_Habitat(MSEBoxModel *bm, char *name, int box, double start, double period, double mult, FILE *llogfp);
//void Scale_Linear_Mortality(MSEBoxModel *bm, int speciesIndex, int cohort, int box, double start, double period, double mult, FILE *llogfp);

void Implement_Strategy(MSEBoxModel *bm, int scenarioIndex, int strategyIndex, FILE *llogfp) {
	int leverID;
	CLAMLeverStructure *lever;

	fprintf(llogfp, "CLAM: Implementing scenario '%s' strategy '%s'\n", CLAMScenarioOptions[scenarioIndex].scenarioOptionName, CLAMScenarioOptions[scenarioIndex].strategies[strategyIndex].strategyName);
	for (leverID = 0; leverID < CLAMScenarioOptions[scenarioIndex].strategies[strategyIndex].numLevers; leverID++) {
		lever = &CLAMScenarioOptions[scenarioIndex].strategies[strategyIndex].strategyLevers[leverID];

		switch (lever->leverType) {
		case PSS_CHANGE:
			Scale_PointSource(bm, lever->targetName, lever->location, lever->start, lever->period, lever->mult, llogfp);
			break;
		case Q_CHANGE:
			Harvest_Scale_Q(bm, lever->groupCode, lever->location, lever->start, lever->period, lever->mult, llogfp);
			break;
		case HABITAT_CHANGE:
			Scale_Habitat(bm, lever->targetName, lever->location, lever->start, lever->period, lever->mult, llogfp);
			break;
		case LINEAR_MORT_CHANGE:
			Scale_Group_Linear_Mortality(bm, lever->groupCode, lever->cohort, lever->location, lever->start, lever->period, lever->mult, llogfp);
			break;
		default:
			quit("atCLAMImplement - Lever type %d is not supported\n", lever->leverType);
			break;
		}
	}

	/* Check to see if we have a lock out period for this strategy */
	if(CLAMScenarioOptions[scenarioIndex].strategies[strategyIndex].lockoutPeriod > 0){
		CLAMScenarioOptions[scenarioIndex].strategies[strategyIndex].timeLockOutStarted = bm->t;
		CLAMScenarioOptions[scenarioIndex].strategies[strategyIndex].isLocked = TRUE;
		fprintf(llogfp, "Locking scenario %d, strategy %d, t = %e\n", scenarioIndex, strategyIndex, bm->t);
	}
}

void SwapDouble(double ***x, double ***y) {
	double **tempArray;

	tempArray = *x;
	*x = *y;
	free2d(tempArray);
}

void SwapInt(int ***x, int ***y) {
	int **tempArray;

	tempArray = *x;
	*x = *y;
	i_free2d(tempArray);
}


/**
 * \Scale the point source.
 *
 *
 * This is done in a slightly different way to the other change arrays.
 * The point source scaling data is not sorted and more than one scaling can take place at any time.
 *
 */
void Scale_PointSource(MSEBoxModel *bm, char *varName, int pssIndex, double start, double period, double mult, FILE *llogfp) {

	double **multArray;
	int **startArray;
	int **periodArray;
	int **changeFlagArray;
	int i, j, variableIndex;

	fprintf(llogfp, "Scale_PointSource: Scaling %s: %s by mult %e from %e over %e\n", bm->pss[pssIndex].name, varName, mult,
			start, period);

	bm->pss[pssIndex].numPssChanges++;
	/* Now increase the array sizes */
	changeFlagArray = Util_Alloc_Init_2D_Int(bm->pss[pssIndex].ts.nv - 1, bm->pss[pssIndex].numPssChanges, FALSE);
	multArray = Util_Alloc_Init_2D_Double(bm->pss[pssIndex].ts.nv - 1, bm->pss[pssIndex].numPssChanges, 1.0);
	startArray = Util_Alloc_Init_2D_Int(bm->pss[pssIndex].ts.nv - 1, bm->pss[pssIndex].numPssChanges, 0);
	periodArray = Util_Alloc_Init_2D_Int(bm->pss[pssIndex].ts.nv - 1, bm->pss[pssIndex].numPssChanges, 0);

	/* Now copy the values across */
	for (i = 0; i < bm->pss[pssIndex].numPssChanges - 1; i++) {
		for (j = 0; j < bm->pss[pssIndex].ts.nv - 1; j++) {
			changeFlagArray[i][j] = bm->pss[pssIndex].pssChange[i][j];
			multArray[i][j] = bm->pss[pssIndex].pssMult[i][j];
			startArray[i][j] = bm->pss[pssIndex].pssStart[i][j];
			periodArray[i][j] = bm->pss[pssIndex].pssPeriod[i][j];
		}
	}

	/* Find the index of the variables in the source list of variables */
	variableIndex = tsIndex(&bm->pss[pssIndex].ts, varName) - 1;
	if (variableIndex == -1)
		quit("Scale_PointSource - The variable %s is not in the point source %s\n", varName, bm->pss[pssIndex].name);

	//fprintf(llogfp, "variableIndex = %d\n", variableIndex);
	/* Now add this change to the end of the array - the rest of the variables will have been set 'off' when the
	 * new array is allocated.
	 */
	changeFlagArray[bm->pss[pssIndex].numPssChanges - 1][variableIndex] = TRUE;
	multArray[bm->pss[pssIndex].numPssChanges - 1][variableIndex] = mult;
	startArray[bm->pss[pssIndex].numPssChanges - 1][variableIndex] = bm->dayt + start;
	periodArray[bm->pss[pssIndex].numPssChanges - 1][variableIndex] = period;

	/* Now change the pointers and free up the old arrays */
	SwapDouble(&bm->pss[pssIndex].pssMult, &multArray);
	SwapInt(&bm->pss[pssIndex].pssStart, &startArray);
	SwapInt(&bm->pss[pssIndex].pssPeriod, &periodArray);
	SwapInt(&bm->pss[pssIndex].pssChange, &changeFlagArray);

	/* After scaling */
	fprintf(llogfp, "After scaling the point source scaling will be:\n");
	for (i = 0; i < bm->pss[pssIndex].numPssChanges; i++) {
		for (j = 0; j < bm->pss[pssIndex].ts.nv - 1; j++) {
			fprintf(llogfp, "bm->pss[pssIndex].pssMult = %e\n", bm->pss[pssIndex].pssMult[i][j]);
			fprintf(llogfp, "bm->pss[pssIndex].psschange = %d\n", bm->pss[pssIndex].pssChange[i][j]);
			fprintf(llogfp, "bm->pss[pssIndex].pssStart = %d\n", bm->pss[pssIndex].pssStart[i][j]);
			fprintf(llogfp, "bm->pss[pssIndex].pssPeriod = %d\n", bm->pss[pssIndex].pssPeriod[i][j]);
		}
	}

}
//
///**
// * Add code to scale linear mortality for each group in each box.
// *
// *	Max_Num_Mort_Changes
// * MortScale[sp][ageclass][b][index][start_id]
// * MortScale[sp][ageclass][b][index][period_id]
// * MortScale[sp][ageclass][b][index][end_id]
// */
//void Scale_Linear_Mortality(MSEBoxModel *bm, int speciesIndex, int cohort, int box, double start, double period, double mult, FILE *llogfp) {
//
//	double *****tempArray;
//	double *****tempArrayPointer;
//	int i;
//
//	fprintf(llogfp, "Scale_Linear_Mortality: Scaling %s:%d in box %d by mult %e from %e over %e\n", FunctGroupArray[speciesIndex].groupCode, cohort, box, mult,
//			start, period);
//	/* Increase the number of scaling events */
//	numMortChanges[speciesIndex][cohort][box]++;
//
//	printf("maxMortChange = %d\n", maxMortChange);
//	if (((int) numMortChanges[speciesIndex][cohort][box]) > maxMortChange) {
//		maxMortChange = numMortChanges[speciesIndex][cohort][box];
//	}
//	/* Allocate a new array */
//	tempArray = Util_Alloc_Init_5D_Double(3, maxMortChange, bm->nbox, 2, bm->K_num_tot_sp, 0);
//
//	/* Copy all existing values across */
////	memcpy(tempArray, LinearMortChange, sizeof(LinearMortChange));
//
//	/* Now insert the new value in the correct place */
//	Util_Copy_Change_Values(bm, LinearMortChange[speciesIndex][cohort][box],
//			tempArray[speciesIndex][cohort][box],
//			numMortChanges[speciesIndex][cohort][box],
//			start, period, mult);
//
//	/* Now do some pointer magic and free up the 'old' array*/
//	tempArrayPointer = LinearMortChange;
//	LinearMortChange = tempArray;
//	free5d(tempArrayPointer);
//
//	for(i = 0; i < ((int) numMortChanges[speciesIndex][cohort][box]); i++){
//		fprintf(llogfp, "LinearMortChange[speciesIndex][cohort][box][%d][start_id = %e\n", i, LinearMortChange[speciesIndex][cohort][box][i][start_id]);
//		fprintf(llogfp, "LinearMortChange[speciesIndex][cohort][box][%d][period_id] = %e\n", i, LinearMortChange[speciesIndex][cohort][box][i][period_id]);
//		fprintf(llogfp, "LinearMortChange[speciesIndex][cohort][box][%d][mult_id] = %e\n", i, LinearMortChange[speciesIndex][cohort][box][i][mult_id]);
//	}
//
//}

//
//void Scale_Temperature(MSEBoxModel *bm, double start, double period, double mult) {
//
//	double **TchangeArray;
//	int i, j, arrayIndex;
//	int beenInserted = FALSE;
//	int clayer;
//	Tchange_max_num++;
//
//	/* Now increase the array sizes */
//	TchangeArray = (double **) alloc2d(3, Tchange_max_num);
//	/* increase the size of the
//	VchangeArray = (double **) alloc2d(3, Tchange_max_num);
//
//	/* Now copy the values across inserting the new value in the correct place.*/
//	arrayIndex = 0;
//	for (i = 0; i < Tchange_max_num - 1; i++) {
//		if (beenInserted == FALSE) {
//			if (Tchange[i][start_id] > start) {
//				/* insert the new values into the array here */
//				TchangeArray[arrayIndex][start_id] = start;
//				TchangeArray[arrayIndex][period_id] = period;
//				TchangeArray[arrayIndex][mult_id] = mult;
//				for(clayer = 0; clayer <  bm->wcnz; clayer++)
//					Vchange[arrayIndex][clayer][tempscale_id] = 1.0;
//				arrayIndex++;
//				beenInserted = TRUE;
//			}
//		}
//		TchangeArray[arrayIndex][start_id] = Tchange[i][start_id];
//		TchangeArray[arrayIndex][period_id] = Tchange[i][period_id];
//		TchangeArray[arrayIndex][mult_id] = Tchange[i][mult_id];
//		arrayIndex++;
//
//	}
//
//	Vmult = Vchange[now_change][clayer][tempscale_id];
//
//	/* Find the index of the variables in the source list of variables */
//	variableIndex = tsIndex(&bm->pss[pssIndex].ts, varName);
//
//	/* Now add this change to the end of the array - the rest of the variables will have been set 'off' when the
//	 * new array is allocated.
//	 */
//	changeFlagArray[bm->pss[pssIndex].numPssChanges - 1][variableIndex] = TRUE;
//	multArray[bm->pss[pssIndex].numPssChanges - 1][variableIndex] = mult;
//	startArray[bm->pss[pssIndex].numPssChanges - 1][variableIndex] = bm->dayt + start;
//	periodArray[bm->pss[pssIndex].numPssChanges - 1][variableIndex] = period;
//
//	/* Now change the pointers and free up the old arrays */
//	SwapDouble(&bm->pss[pssIndex].pssMult, &multArray);
//
//	SwapInt(&bm->pss[pssIndex].pssStart, &startArray);
//	SwapInt(&bm->pss[pssIndex].pssPeriod, &periodArray);
//	SwapInt(&bm->pss[pssIndex].pssChange, &changeFlagArray);
//
//}

void Scale_Salinity(MSEBoxModel *bm, double start, double period, double mult) {

}
//
///**
// * \brief Scale the specified habitat.
// */
//void Scale_Habitat(MSEBoxModel *bm, char *name, int box, double start, double period, double mult, FILE *llogfp) {
//
//	int habitatIndex;
//	double ***tempArray;
//	double ***tempArrayPointer;
//	int i;
//
//	fprintf(llogfp, "Scale_Habitat: Scaling %s\n", name);
//
//	/* First - work out which habitat we are changing */
//	if (strcmp(name, "REEF") == 0) {
//		habitatIndex = reef_id;
//	} else if (strcmp(name, "FLAT") == 0) {
//		habitatIndex = flat_id;
//	} else if (strcmp(name, "SOFT") == 0) {
//		habitatIndex = soft_id;
//	} else {
//		quit("Scale_Habitat: Habitat %s is not recognised. Must be either REEF, FLAT or SOFT\n");
//	}
//
//	/* Turn on habitat degradation in this box */
//	Box_degradedi[box] = TRUE;
//	bm->flagdegrade = TRUE;
//
//	/* Check to see if there is just a single dummy value in the array - if so just overwrite it */
//	if(BEDchange_max_num == 1 && BEDchange[habitatIndex][0][start_id] == 0){
//		BEDchange[habitatIndex][0][start_id] = start + bm->dayt;
//		BEDchange[habitatIndex][0][period_id] = period;
//		BEDchange[habitatIndex][0][mult_id] = mult;
//		return;
//	}
//	/* Otherwise increase the size and insert it into the correct location */
//	BEDchange_max_num++;
//
//	fprintf(llogfp, "Scale_Habitat: Scaling %s:%d in box %d by mult %e from %e over %e\n", name, habitatIndex, box, mult,
//				start, period);
//
//	/* Now allocate a temp array to store the updated values */
//	tempArray = Util_Alloc_Init_3D_Double(3, BEDchange_max_num, K_num_sed_type, 0);
//
//	/* Now copy the values across inserting the new value in the correct place.*/
//	Util_Copy_Change_Values(bm, BEDchange[habitatIndex], tempArray[habitatIndex],
//			BEDchange_max_num, start, period, mult);
//
//	/* Now do some pointer magic and free up the 'old' array*/
//	tempArrayPointer = BEDchange;
//	BEDchange = tempArray;
//	free3d(tempArrayPointer);
//
//	for(i = 0; i < BEDchange_max_num; i++){
//		fprintf(llogfp, "BEDchange[habitatIndex][%d][start_id = %e\n", i, BEDchange[habitatIndex][i][start_id]);
//		fprintf(llogfp, "BEDchange[habitatIndex][%d][period_id] = %e\n", i, BEDchange[habitatIndex][i][period_id]);
//		fprintf(llogfp, "BEDchange[habitatIndex][%d][mult_id] = %e\n", i, BEDchange[habitatIndex][i][mult_id]);
//	}
//
//}
