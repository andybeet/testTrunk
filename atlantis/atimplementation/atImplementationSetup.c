/**
 * 	\ingroup atImplementationLib
 *	\file atImplementationSetup.c
 *
 *  Created on: 07/07/2009
 *      Author: Bec Gorton
 */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <atImplementation.h>
#include <convertXML.h>


double ****TACchange;

/**
 * \brief Initialise the management implementation module.
 *
 */
void Implementation_Init(MSEBoxModel *bm, FILE *llogfp){

	char convertedXMLFileName[50];

	if ( !bm->flag_fisheries_on ){
		return;
	}

	strcpy(convertedXMLFileName, "implementation.xml");

	/* Convert the input file to XML */
	Convert_Implementation_To_XML(bm, bm->fishprmIfname, convertedXMLFileName);

	readModelManageImplParameters(bm, convertedXMLFileName);
}

/**
 * \brief Free the memory allocated by the Implementation library.
 *
 *
 */
void Implementation_Free(MSEBoxModel *bm){

	if ( !bm->flag_fisheries_on ){
		return;
	}
	free4d(TACchange);
}
