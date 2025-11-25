/**
 * \defgroup atUtil atUtil
 *
 * The utility module of Atlantis
 *
 */

/*
 * atFisheryStruct.h
 *
 *
 * The structure that holds the info about each fishery.
 *
 *
 *  Created on: 12/08/2011
 *      Author: bec
 */

#ifndef ATFISHERYSTRUCT_H_
#define ATFISHERYSTRUCT_H_

/**
 *	The structure that is used to store information about each functional group
 *	in the model.
 */
typedef struct {
    char fisheryCode[20];					/**< The code used to identify this fishery in the input files */
	char name[200];						/**< The verbose name */
	int index;

	int totalEffortTracer;
	int averageCatchSizeTracer;
	int isRec;

}FisheryStruct;
extern FisheryStruct *FisheryArray;

#endif /* ATFISHERYSTRUCT_H_ */
