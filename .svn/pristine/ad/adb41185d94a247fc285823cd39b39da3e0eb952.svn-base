/**************************************************************************//**
 \file
 \brief The code to interface with the CLAM model.
 \ingroup atCLAMLink

 File:       atCLAMManage.c
 Created:    20/11/2009
 Author:     Bec Gorton,   CSIRO Marine Laboratories,  2009
 Purpose:    The code to determine the 'best' management option based on the output from the CLAM.
 Note:

 The algorithm for determining the 'best' management scenario is in a spreadsheet in the code repository.
 This was developed by Rich Little, Scott Condie, Beth Fulton and Bec Gorton.

 Documentation of the algorithm is on the Atlantis wiki at:

 Revisions:

 07-05-2010 Bec Gorton
 Changed the code that calculates the corrected distance to multiple by the normalised variance. The 'best' scenario selection
 was not quite working - it was tending to select options with the highest variance.
 Updated spreadsheet in the project repository for this project and linked from the wiki.

 06-08-2010 Bec Gorton
 Changed the code to consider all indicators when selecting the best strategy.
 Added a lockout period for each strategy.
 Added a cost for each strategy.

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
#include <atCLAMLinkPrivate.h>
#include <atAssessLib.h>

static void Run_CLAM(MSEBoxModel *bm, FILE *llogfp);
static double Get_Expected_Value(MSEBoxModel *bm, int indicatorIndex, int scenarioIndex, int strategyIndex);
static double Get_Variance(MSEBoxModel *bm, int indicatorIndex, int scenarioIndex, int strategyIndex, double E);
static double Get_Indicator_SD(MSEBoxModel *bm, int indicatorID, int b);
static double Get_Box_Indicator_Value(MSEBoxModel *bm, int indicatorIndex, int b);
static void Determine_Indicator_Weightings(MSEBoxModel *bm);

/**
 *
 * \brief Check the CLAM indicators.
 *
 * If any of the indicators is outside of the allowed range the we need to call the CLAM
 * for assistance in making a decision about which course of action to take to remedy the problem.
 *
 * At the moment we have the limitation of really only allowing a single region defined by the intense sampling sites.
 * We should be looping over a number of regions within the model.
 *
 */
void CLAM_Check_Indicators(MSEBoxModel *bm, FILE *llogfp) {
	int i, b, bb;
	double indicatorValue, change, requiredChange;
	double tempIndicatorValue;
	int runCLAM = FALSE; /* Assume that all the indicators are within the target range */
	double SD;
	//double maxScenarioImpact;
	//int numOutputStates;

	//maxScenarioImpact = -1;
	for (i = 0; i < numIndicators; i++) {

		/* By default we are not running the scenarios attached to this indicator */
		indicators[i].runScenarios = FALSE;
		//numOutputStates = CLAMOutputOptions[indicators[i].CLAMOutputIndex].numStates;
		requiredChange = 0;

		/* How we calculate the required change depends on the type of indicator.
		 * If its a standard dev indicator then we calculate the standard deviation over all boxes.
		 * If its a percent change indicator then we calculate the indicator value in a number of special boxes.
		 */
		/* For each site of interest - get the required change. Just look at the intense sites.*/
		for (bb = 0; bb < bm->nsbox_i; bb++) {
			b = bm->nsboxes[bb];
			tempIndicatorValue = Get_Box_Indicator_Value(bm, i, b);

			fprintf(bm->logFile, "CLAM - Indicator %s value is %e in box %d\n", indicators[i].CLAMOutputName, tempIndicatorValue, b);

			/* Check that the indicator is inside the allowed range */
			if (tempIndicatorValue < indicators[i].minValue || tempIndicatorValue > indicators[i].maxValue) {

				/* get the required change */
				change = indicators[i].idealValue - tempIndicatorValue;

				/* Save if its the max required change - might need other ways of doing this - might want to do weighted average over sites */
				if (change > requiredChange || requiredChange == 0) {
					requiredChange = change;
				}
			}
		}
		if (requiredChange != 0) {
			runCLAM = TRUE;

		} else {
			fprintf(bm->logFile, "CLAM run not required\n");
		}
	}

	/* Ok one of our indicators is outside its acceptable range - now take all of the other indicators into account
	 * so we don't select a strategy that will have nasty effects on the rest of them
	 */
	if (runCLAM) {
		printf("numIndicators = %d\n", numIndicators);
		for (i = 0; i < numIndicators; i++) {

			/* By default we are not running the scenarios attached to this indicator */
			indicators[i].runScenarios = TRUE;
			//numOutputStates = CLAMOutputOptions[indicators[i].CLAMOutputIndex].numStates;
			requiredChange = 0;

			/* How we calculate the required change depends on the type of indicator.
			 * If its a standard dev indicator then we calculate the standard deviation over all boxes.
			 * If its a percent change indicator then we calculate the indicator value in a number of special boxes.
			 */
			/* For each site of interest - get the required change. Just look at the intense sites.*/
			for (bb = 0; bb < bm->nsbox_i; bb++) {
				b = bm->nsboxes[bb];
				tempIndicatorValue = Get_Box_Indicator_Value(bm, i, b);

				printf("name = %s, tempIndicatorValue = %e\n", indicators[i].CLAMOutputName, tempIndicatorValue);
				/* Check that the indicator is inside the allowed range */
				//if (tempIndicatorValue < indicators[i].minValue || tempIndicatorValue > indicators[i].maxValue) {

				/* get the required change */
				change = indicators[i].idealValue - tempIndicatorValue;

				/* Save if its the max required change - might need other ways of doing this - might want to do weighted average over sites */
				if (change > requiredChange || requiredChange == 0) {
					requiredChange = change;
					indicatorValue = tempIndicatorValue;

					/* if this indicator values are returned from the CLAM in SD then get the SD of this indicator so we can
					 * work out the required change as a prop of the SD.
					 */
					if (CLAMOutputOptions[indicators[i].CLAMOutputIndex].outputType == STANDARD_DEV) {

						SD = Get_Indicator_SD(bm, i, b);
					}
				}
				//}
			}
			if (requiredChange != 0) {
				/* Assign the impact values for this indicator */
				// allocatedImpactWeightings(bm, i, requiredChange);
				if (CLAMOutputOptions[indicators[i].CLAMOutputIndex].outputType == PERCENT) {
					/* Calculate the percentage change required */
					requiredChange = requiredChange / indicatorValue * 100;
					indicators[i].targetChange = requiredChange;
				} else {

					requiredChange = requiredChange / SD;
					indicators[i].targetChange = requiredChange;
					indicators[i].currentValue = indicatorValue;
				}
				/* log the info to the log file */
				fprintf(llogfp, "CLAM: Indicator %s, indicatorValue = %e, RequiredChange = %e, SD =  %e\n", indicators[i].CLAMOutputName, indicatorValue,
						requiredChange, SD);

				/* Add the scenarios for this indicator to the list of scenarios to run */
				indicators[i].runScenarios = TRUE;
				runCLAM = TRUE;
			}
		}
	}
	/* Now run the CLAM and implement the 'best' strategy at this site.*/
	if (runCLAM == TRUE) {
		Run_CLAM(bm, llogfp);
	}

	/* Now work out when we should run the CLAM again */
	bm->clamDay = (bm->clamDay + CLAMFreq);
}

/**
 *
 */
static void Run_CLAM(MSEBoxModel *bm, FILE *llogfp) {

	int linkID, strategyIndex;
	int i;
	int clamOutputIndicatorIndex;
	int scenarioID, scenarioIndex;
	int indicatorIndex;
	double targetChange, minDistance, sum, averageDistance;
	int count;
	double varSum;

	/* Clear the previous outcome values - just in case we have some strategies locked. */
	Util_Init_4D_Double(outcomeValueArray, numClamIndicators, numScenarioOptions, maxStrategyOptions, MAXNUMOUTPUTSTATES, 0.0);
	Util_Init_3D_Double(expectedValue, numScenarioOptions, maxStrategyOptions, numClamIndicators, 0.0);
	Util_Init_3D_Double(variance, numScenarioOptions, maxStrategyOptions, numClamIndicators, 0.0);
	Util_Init_3D_Double(distance, numScenarioOptions, maxStrategyOptions, numClamIndicators, 0.0);
	Util_Init_2D_Double(euclidDistance, numScenarioOptions, maxStrategyOptions, 0.0);

	/* log the info to the log file */
	fprintf(llogfp, "CLAM: Running the CLAM at %.10g.\n", bm->t / 86400.0);

	/* Allocate weightings to the indicators - this depends on the manager type*/
	Determine_Indicator_Weightings(bm);

	/* Build an array of all of the scenarios that we are interested in running
	 *
	 * Select all of the output indicators that we might potentially be interested in
	 */
	/* Turn off all CLAM output options */
	for (indicatorIndex = 0; indicatorIndex < numClamIndicators; indicatorIndex++) {
		CLAMOutputOptions[indicatorIndex].optionSelected = 0;
	}
	/* then 'turn on' all of the CLAM indicators that we are interested in */
	for (i = 0; i < numIndicators; i++) {
		if (indicators[i].runScenarios == TRUE) {
			clamOutputIndicatorIndex = indicators[i].CLAMOutputIndex;
			CLAMOutputOptions[clamOutputIndicatorIndex].optionSelected = 1;
		}
	}

	/* Initialise the scenarios */
	for (i = 0; i < numScenarioOptions; i++) {
		/* Flag all scenarios as 'off' */
		CLAMScenarioOptions[i].doRun = FALSE;
		/* initialise the scenario strategies to 0*/
		CLAMScenarioOptions[i].optionSelected = 0;
	}
	/* Now work out what scenarios we are interested in running */
	for (i = 0; i < numIndicators; i++) {

		/* If we are interested in changing this indicator 'turn on' all of the scenarios that are linked to the
		 * indicators - these are currently read in from the input file
		 */
		if (indicators[i].runScenarios == TRUE) {
			for (linkID = 0; linkID < indicators[i].numScenarioLinks; linkID++) {
				scenarioID = indicators[i].scenarioLinks[linkID];

				/* Flag this scenario as 'on' */
				CLAMScenarioOptions[scenarioID].doRun = TRUE;
				//printf("CLAMScenarioOptions[%d].doRun \n", scenarioID);
			}
		}
	}

	/* Now run each strategy within each relevant scenario and store the outcome values for each indicator
	 * we are interested in
	 */
	for (scenarioIndex = 0; scenarioIndex < numScenarioOptions; scenarioIndex++) {
		if (CLAMScenarioOptions[scenarioIndex].doRun == TRUE) {
			//printf("Trying scenario %s\n", CLAMScenarioOptions[scenarioIndex].scenarioOptionName);
			for (strategyIndex = 0; strategyIndex < CLAMScenarioOptions[scenarioIndex].maxOptionValue; strategyIndex++) {

				/* check that this option is enabled */
				if (CLAMScenarioOptions[scenarioIndex].strategies[strategyIndex].isEnabled > 0) {
					/* Check that this strategy hasen't been locked out */
					if (CLAMScenarioOptions[scenarioIndex].strategies[strategyIndex].isLocked == TRUE) {
						/* if its been locked then don't bother trying this strategy */
						if (bm->t < (CLAMScenarioOptions[scenarioIndex].strategies[strategyIndex].timeLockOutStarted
								+ CLAMScenarioOptions[scenarioIndex].strategies[strategyIndex].lockoutPeriod)) {
							continue;
						}
						/* else the lock is finished and we can consider this strategy again */
						CLAMScenarioOptions[scenarioIndex].strategies[strategyIndex].isLocked = FALSE;
						fprintf(llogfp, "UnLocking scenario %d, strategy %d, t = %e\n", scenarioIndex, strategyIndex, bm->t);
					}

					/* Select this strategy - this is done so its flagged in the CLAM input file */
					CLAMScenarioOptions[scenarioIndex].optionSelected = CLAMScenarioOptions[scenarioIndex].strategies[strategyIndex].CLAMStrategyID;

					/* Start the CLAM - will write out the required files and fire off the process */
					Start_CLAM(bm, inputFileName);

					/* Parse the outcome values for each indicator */
					Parse_CLAM_OutputFile(bm, scenarioIndex, strategyIndex);
				}
			}
		}
	}

	fprintf(llogfp, "CLAM: CLAM run - determining the 'best' strategy\n");

	/* Now we have got all the data we need from the CLAM
	 *
	 * - We just need to select the 'best' strategy
	 */

	/* Calculate the expected value and variance of each scenario-strategy against each indicator */

	/* Now work out what scenarios we are interested in running */
	for (indicatorIndex = 0; indicatorIndex < numIndicators; indicatorIndex++) {
		if (indicators[indicatorIndex].runScenarios == TRUE) {
			targetChange = indicators[indicatorIndex].targetChange;
			sum = 0; /* The sum of the non-trivial distances */
			count = 0; /* The count of the non-trivial distances */

			fprintf(bm->logFile, "CLAM output states = '");
			for (i = 0; i < CLAMOutputOptions[indicators[indicatorIndex].CLAMOutputIndex].numStates; i++) {
				fprintf(bm->logFile, "%e ", CLAMOutputOptions[indicators[indicatorIndex].CLAMOutputIndex].OutputStates[i]);
			}
			fprintf(bm->logFile, "'\n");

			for (scenarioIndex = 0; scenarioIndex < numScenarioOptions; scenarioIndex++) {
				if (CLAMScenarioOptions[scenarioIndex].doRun == TRUE) {
					for (strategyIndex = 0; strategyIndex < CLAMScenarioOptions[scenarioIndex].maxOptionValue; strategyIndex++) {
						if (CLAMScenarioOptions[scenarioIndex].strategies[strategyIndex].isLocked == FALSE
								&& CLAMScenarioOptions[scenarioIndex].strategies[strategyIndex].isEnabled > 0) {

							/* Calculate the estimated value */
							expectedValue[scenarioIndex][strategyIndex][indicatorIndex] = Get_Expected_Value(bm, indicatorIndex, scenarioIndex, strategyIndex);

							variance[scenarioIndex][strategyIndex][indicatorIndex] = Get_Variance(bm, indicatorIndex, scenarioIndex, strategyIndex,
									expectedValue[scenarioIndex][strategyIndex][indicatorIndex]);
						}
					}

					/* get the total to the variance so we can normal it */
					varSum = 0.0;
					for (strategyIndex = 0; strategyIndex < CLAMScenarioOptions[scenarioIndex].maxOptionValue; strategyIndex++) {
						if (CLAMScenarioOptions[scenarioIndex].strategies[strategyIndex].isLocked == FALSE
								&& CLAMScenarioOptions[scenarioIndex].strategies[strategyIndex].isEnabled > 0) {
							varSum += variance[scenarioIndex][strategyIndex][indicatorIndex];
						}
					}
					fprintf(llogfp, "variance sum = %e\n", varSum);
					for (strategyIndex = 0; strategyIndex < CLAMScenarioOptions[scenarioIndex].maxOptionValue; strategyIndex++) {
						if (CLAMScenarioOptions[scenarioIndex].strategies[strategyIndex].isLocked == FALSE
								&& CLAMScenarioOptions[scenarioIndex].strategies[strategyIndex].isEnabled > 0) {

							double normaliseVar = variance[scenarioIndex][strategyIndex][indicatorIndex] / varSum;

							fprintf(llogfp, "normaliseVar = %e\n", normaliseVar);
							/* Calculate the distance to the required value */
							if (expectedValue[scenarioIndex][strategyIndex][indicatorIndex] == 0 && variance[scenarioIndex][strategyIndex][indicatorIndex] == 0) {
								distance[scenarioIndex][strategyIndex][indicatorIndex] = 0.0;
							} else if (variance[scenarioIndex][strategyIndex][indicatorIndex] == 0) {
								distance[scenarioIndex][strategyIndex][indicatorIndex] = fabs(targetChange
										- expectedValue[scenarioIndex][strategyIndex][indicatorIndex]);
								sum = sum + distance[scenarioIndex][strategyIndex][indicatorIndex];
								count++;
							} else {
								distance[scenarioIndex][strategyIndex][indicatorIndex] = //variance[scenarioIndex][strategyIndex][indicatorIndex] / varSum *
										fabs(targetChange - expectedValue[scenarioIndex][strategyIndex][indicatorIndex]) / (fabs(targetChange) * (1
												- variance[scenarioIndex][strategyIndex][indicatorIndex]));
								//variance[scenarioIndex][strategyIndex][indicatorIndex] / fabs(targetChange));
								//fprintf(llogfp, "fabs(targetChange - expectedValue[scenarioIndex][strategyIndex][indicatorIndex]) = %e\n", fabs(targetChange - expectedValue[scenarioIndex][strategyIndex][indicatorIndex]));
								//fprintf(llogfp, "(fabs(targetChange) * (1 + normaliseVar)) = %e\n", (fabs(targetChange) * (1  - variance[scenarioIndex][strategyIndex][indicatorIndex])));
								sum = sum + distance[scenarioIndex][strategyIndex][indicatorIndex];
								count++;
							}

							fflush(llogfp);
							fprintf(llogfp, "CLAM: Indicator %s, Scenario %s, strategy %d, target change = %e, E = %e, Var = %e, D = %e\n",
									indicators[indicatorIndex].CLAMOutputName, CLAMScenarioOptions[scenarioIndex].scenarioOptionName, strategyIndex,
									targetChange, expectedValue[scenarioIndex][strategyIndex][indicatorIndex],
									variance[scenarioIndex][strategyIndex][indicatorIndex], distance[scenarioIndex][strategyIndex][indicatorIndex]);

						}
					}
				}
			}

			/* calculate the average non-trivial distance */
			if (count > 0)
				averageDistance = sum / count;
			else
				averageDistance = 0;

			fprintf(llogfp, "CLAM: Indicator %s, Average Distance = %e\n", indicators[indicatorIndex].CLAMOutputName, averageDistance);
			//printf("CLAM: Indicator %s, Average Distance = %e\n", indicators[indicatorIndex].CLAMOutputName, averageDistance);

			/* Now assign the trivial case (where E == 0 and V == 0) the average non-trivial distance value */
			/* This calculates the corrected distance */
			for (scenarioIndex = 0; scenarioIndex < numScenarioOptions; scenarioIndex++) {
				if (CLAMScenarioOptions[scenarioIndex].doRun == TRUE) {
					for (strategyIndex = 0; strategyIndex < CLAMScenarioOptions[scenarioIndex].maxOptionValue; strategyIndex++) {
						if (CLAMScenarioOptions[scenarioIndex].strategies[strategyIndex].isLocked == FALSE
								&& CLAMScenarioOptions[scenarioIndex].strategies[strategyIndex].isEnabled > 0) {
							if (expectedValue[scenarioIndex][strategyIndex][indicatorIndex] == 0 && variance[scenarioIndex][strategyIndex][indicatorIndex] == 0) {
								distance[scenarioIndex][strategyIndex][indicatorIndex] = averageDistance;

								fprintf(llogfp, "CLAM: Indicator %s, Scenario %s, strategy %d, CorrectedDistance = %e\n",
										indicators[indicatorIndex].CLAMOutputName, CLAMScenarioOptions[scenarioIndex].scenarioOptionName, strategyIndex,
										distance[scenarioIndex][strategyIndex][indicatorIndex]);
								printf("CLAM: Indicator %s, Scenario %s, strategy %d, CorrectedDistance = %e\n", indicators[indicatorIndex].CLAMOutputName,
										CLAMScenarioOptions[scenarioIndex].scenarioOptionName, strategyIndex,
										distance[scenarioIndex][strategyIndex][indicatorIndex]);
							}
						}
					}
				}
			}

		} /* end of the runScenarios == TRUE if statement */
	}/* End of the for loop */

	/* Now combine the values for each indicator based on their weightings and find the min value - this is the strategy we should run*/
	fprintf(bm->logFile, "Combine indicator values based on weightings, scale by strategy weightings and find min value\n");
	minDistance = HUGE_VAL;
	bestScenarioIndex = -1;
	bestStrategyIndex = -1;
	for (scenarioIndex = 0; scenarioIndex < numScenarioOptions; scenarioIndex++) {
		if (CLAMScenarioOptions[scenarioIndex].doRun == TRUE) {
			for (strategyIndex = 0; strategyIndex < CLAMScenarioOptions[scenarioIndex].maxOptionValue; strategyIndex++) {
				if (CLAMScenarioOptions[scenarioIndex].strategies[strategyIndex].isLocked == FALSE
						&& CLAMScenarioOptions[scenarioIndex].strategies[strategyIndex].isEnabled > 0) {
					sum = 0;
					for (indicatorIndex = 0; indicatorIndex < numIndicators; indicatorIndex++) {
						//						fprintf(llogfp, "indicatorWeightings[%d] = %e\n", indicatorIndex, indicatorWeightings[indicatorIndex]);
						//						fprintf(llogfp, "distance[scenarioIndex][strategyIndex][indicatorIndex] = %e\n", distance[scenarioIndex][strategyIndex][indicatorIndex]);
						if (indicators[indicatorIndex].runScenarios == TRUE) {
							sum = sum + (pow(distance[scenarioIndex][strategyIndex][indicatorIndex], 2) * indicatorWeightings[indicatorIndex]);
						}
					}
					fprintf(llogfp, "Scenario %d, strategy %d. Combined weighted distance before scaling - %e, isEnabled = %d\n", scenarioIndex, strategyIndex,
							sum, CLAMScenarioOptions[scenarioIndex].strategies[strategyIndex].isEnabled);
					sum = sum * CLAMScenarioOptions[scenarioIndex].strategies[strategyIndex].strategyCost;
					euclidDistance[scenarioIndex][strategyIndex] = sum;

					fprintf(llogfp, "Distance after scaling = %e, cost = %e\n", euclidDistance[scenarioIndex][strategyIndex],
							CLAMScenarioOptions[scenarioIndex].strategies[strategyIndex].strategyCost);

					if (sum < minDistance) {
						minDistance = sum;
						/* save the 'best' outcome so we can 'implement' it */
						bestStrategyIndex = strategyIndex;
						bestScenarioIndex = scenarioIndex;
					}
				}
			}
		}
	}

	fprintf(llogfp, "CLAM: bestScenarioIndex %d, bestStrategyIndex %d\n", bestScenarioIndex, bestStrategyIndex);
	if (bestScenarioIndex > 0 && bestStrategyIndex > 0)
		Implement_Strategy(bm, bestScenarioIndex, bestStrategyIndex, llogfp);
	else {
		fprintf(llogfp, "\nCLAM - Cannot find a 'good' scenario to implement so doing nothing.\n\n");
	}
	fprintf(llogfp, "\nCLAM work done\n\n");
}

/**
 * \brief Calculate and return the expect value of this strategy for the given indicator.
 */
static double Get_Expected_Value(MSEBoxModel *bm, int indicatorIndex, int scenarioIndex, int strategyIndex) {
	double sum;
	int index;
	sum = 0.0;

	fprintf(bm->logFile, "Scenario %s, strategy %s \nCLAM outcome array = '", CLAMScenarioOptions[scenarioIndex].scenarioOptionName,
			CLAMScenarioOptions[scenarioIndex].strategies[strategyIndex].strategyName);
	for (index = 0; index < CLAMOutputOptions[indicators[indicatorIndex].CLAMOutputIndex].numStates; index++) {
		sum = sum + (outcomeValueArray[indicatorIndex][scenarioIndex][strategyIndex][index]
				* CLAMOutputOptions[indicators[indicatorIndex].CLAMOutputIndex].OutputStates[index]);

		fprintf(bm->logFile, "%e ", outcomeValueArray[indicatorIndex][scenarioIndex][strategyIndex][index]);

	}

	fprintf(bm->logFile, "'\n");
	return sum;
}

/*
 * \brief
 *
 * Calculate the variance of this strategy for the given indicator.
 *
 */
static double Get_Variance(MSEBoxModel *bm, int indicatorIndex, int scenarioIndex, int strategyIndex, double E) {
	double sum;
	int index;
	sum = 0.0;
	for (index = 0; index < MAXNUMOUTPUTSTATES; index++) {
		sum = sum + (outcomeValueArray[indicatorIndex][scenarioIndex][strategyIndex][index]
				* CLAMOutputOptions[indicators[indicatorIndex].CLAMOutputIndex].OutputStates[index]
				* CLAMOutputOptions[indicators[indicatorIndex].CLAMOutputIndex].OutputStates[index]);
	}

	/* Subtract the expected value ^2*/
	sum = sum - (pow(E, 2));
	return sum;
}

/**
 * \brief Get the standard dev from the assessment module for the given indicator.
 *
 * This code
 */
static double Get_Indicator_SD(MSEBoxModel *bm, int indicatorIndex, int b) {

	/* Work out how to get the indicator value out of the assessment module. */
	if (strcmp(indicators[indicatorIndex].CLAMOutputName, "Lake Salinity") == 0) {
		return Assess_Get_Physical_Sampled_SD(bm, salinity_id, b);
	}
	if (strcmp(indicators[indicatorIndex].CLAMOutputName, "Lake Temperature") == 0) {
		return Assess_Get_Physical_Sampled_SD(bm, temperature_id, b);
	}
	if (strcmp(indicators[indicatorIndex].CLAMOutputName, "Lake Total Suspended Solids") == 0) {
		return Assess_Get_Physical_Sampled_SD(bm, TSS_id, b);
	}
	if (strcmp(indicators[indicatorIndex].CLAMOutputName, "Lake Dissolved Oxygen") == 0) {
		return Assess_Get_Physical_Sampled_SD(bm, o2_wc_id, b);
	}
	if (strcmp(indicators[indicatorIndex].CLAMOutputName, "Lake Total Nitrogen") == 0) {
		return Assess_Get_Physical_Sampled_SD(bm, din_wc_id, b);
	}
	quit("Get_Indicator_SD: indicator %d not recognised\n", indicatorIndex);
	return -1;

}

/**
 * \brief Get the value of the indicator from the assessment module.
 *
 *	The assessment module will be told to sample each site many times a year to give us a nice array of values.
 *
 *
 */
//static double Get_Box_Indicator_Value(MSEBoxModel *bm, int indicatorIndex, int b) {
//
//	/* Work out how to get the indicator value out of the assessment module. */
//
//	if (strcmp(indicators[indicatorIndex].CLAMOutputName, "Lake Total Nitrogen") == 0) {
//		/* Do nothing yet */
//	}
//	if (strcmp(indicators[indicatorIndex].CLAMOutputName, "Lake Total Suspended Solids") == 0) {
//		return Assess_Get_Physical_Sampled_Value(bm, TSS_id, b);
//	}
//	if (strcmp(indicators[indicatorIndex].CLAMOutputName, "Lake Salinity") == 0) {
//		return Assess_Get_Physical_Sampled_Value(bm, salinity_id, b);
//	}
//	if (strcmp(indicators[indicatorIndex].CLAMOutputName, "Lake Pathogens") == 0) {
//		/* Do nothing yet */
//	}
//	if (strcmp(indicators[indicatorIndex].CLAMOutputName, "Lake Temperature") == 0) {
//		return Assess_Get_Physical_Sampled_Value(bm, temperature_id, b);
//	}
//	if (strcmp(indicators[indicatorIndex].CLAMOutputName, "Lake pH") == 0) {
//		/* Do nothing yet */
//	}
//	if (strcmp(indicators[indicatorIndex].CLAMOutputName, "Lake Dissolved Oxygen") == 0) {
//		return Assess_Get_Physical_Sampled_Value(bm, o2_wc_id, b);
//	}
//	quit("Get_Box_Indicator_Value: Monitoring of indicator %s not yet implemented\n", indicators[indicatorIndex].CLAMOutputName);
//
//	return -1;
//}


/**
 * \brief Get the value of the indicator from the assessment module.
 *
 *	The assessment module will be told to sample each site many times a year to give us a nice array of values.
 *
 *
 */
static double Get_Box_Indicator_Value(MSEBoxModel *bm, int indicatorIndex, int b) {

	/* Work out how to get the indicator value out of the assessment module. */

	if (strcmp(indicators[indicatorIndex].CLAMOutputName, "Lake Total Nitrogen") == 0) {
		return Assess_Get_Physical_Sampled_Value(bm, din_wc_id, b);
	}
	if (strcmp(indicators[indicatorIndex].CLAMOutputName, "Lake Algae") == 0) {
		return Assess_Get_Process_Sampled_Value(bm, 2, b);
	}
	if (strcmp(indicators[indicatorIndex].CLAMOutputName, "Lake Total Suspended Solids") == 0) {
		return Assess_Get_Physical_Sampled_Value(bm, TSS_id, b);
	}
	if (strcmp(indicators[indicatorIndex].CLAMOutputName, "Lake Salinity") == 0) {
		return Assess_Get_Physical_Sampled_Value(bm, salinity_id, b);
	}
	//	if (strcmp(indicators[indicatorIndex].CLAMOutputName, "Lake Pathogens") == 0) {
	//		/* Do nothing yet */
	//	}
	if (strcmp(indicators[indicatorIndex].CLAMOutputName, "Lake Temperature") == 0) {
		return Assess_Get_Physical_Sampled_Value(bm, temperature_id, b);
	}
	//	if (strcmp(indicators[indicatorIndex].CLAMOutputName, "Lake pH") == 0) {
	//		/* Do nothing yet */
	//	}
	if (strcmp(indicators[indicatorIndex].CLAMOutputName, "Lake Dissolved Oxygen") == 0) {
		return Assess_Get_Physical_Sampled_Value(bm, o2_wc_id, b);
	}
	if (strcmp(indicators[indicatorIndex].CLAMOutputName, "Change in Seagrass Area") == 0) {
		return Assess_Get_Physical_Sampled_Value(bm, o2_wc_id, b);
	}
	quit("Get_Box_Indicator_Value: Monitoring of indicator %s not yet implemented\n", indicators[indicatorIndex].CLAMOutputName);

	return -1;
}


static void Determine_Indicator_Weightings(MSEBoxModel *bm) {
	int indicatorIndex;
	//double weighting;

	//fprintf(bm->logFile, "CLAM: Assigning weighting to indicators\n");
	/* Assign dummy values for now */
	for (indicatorIndex = 0; indicatorIndex < numIndicators; indicatorIndex++) {

		switch (CLAMManagerType) {
		case RANDOM_MANAGER:
			indicatorWeightings[indicatorIndex] = drandom(0.0, 1.0);

			break;
		case REACTIVE_MANAGER:
			/* Also want to allow for distance to ideal value normalised by the indicator value */
			indicatorWeightings[indicatorIndex] = (indicators[indicatorIndex].idealValue - indicators[indicatorIndex].currentValue)
					/ indicators[indicatorIndex].currentValue;
			break;
		case SET_MANAGER:
			indicatorWeightings[indicatorIndex] = indicators[indicatorIndex].indicatorWeighting;
			break;
		}
		//fprintf(bm->logFile, "CLAM: Indicator %s has weighting %e\n", indicators[indicatorIndex].CLAMOutputName, indicatorWeightings[indicatorIndex]);
	}
}

//
///**
// * \brief Calculate and return the standard deviation of the values in the given array.
// *
// */
//static double Calculate_Array_SD(FILE *llogfp, double *valueArray, int size) {
//
//	int i;
//	double mean, sum;
//	double step1, returnValue;
//
//	sum = 0;
//	for (i = 0; i < size; i++) {
//		sum = sum + valueArray[i];
//		fprintf(llogfp, "valueArray[i] = %e\n", valueArray[i]);
//	}
//	mean = sum / size;
//	//    printf("size = %d\n", size);
//	//    printf("mean = %e\n", mean);
//
//	step1 = 0;
//	for (i = 0; i < size; i++) {
//		step1 += pow((valueArray[i] - mean), 2);
//	}
//	returnValue = step1 / size;
//	return returnValue;
//
//}
