/*
 * atBrokerLinkInit.c
 *
 *  Created on: Jun 25, 2012
 *      Author: bec
 */
/*  Open library routines */
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
 * /brief This function reads in the config for the broker linkage.
 *

 *
 */
int Linkage_Read_XML(MSEBoxModel *bm, char *fileName) {

	xmlDocPtr inputDoc;
	xmlXPathObjectPtr polygonList;
	int polygonIndex, numPolygons;
	xmlNodePtr polygonNode;
	char url[STRLEN];

	if (strlen(fileName) == 0) {
		fprintf(stderr,
				"\n\nERROR. Util_Read_Linkage_XML: Trying to read in the config file for the broker linkage information, but no file was provided.\n\n");
		Util_Usage();
		exit(-1);
	}
	if (strstr(fileName, ".xml") == NULL) {
		fprintf(stderr, "\n\nERROR. Util_Read_Linkage_XML: Trying to read in the broker linkage file. This should be read in from an XML input file\n\n");
		Util_Usage();
		exit(-1);
	}

	inputDoc = xmlReadFileDestFolder(bm->destFolder, fileName, NULL, 0);
	if (!inputDoc) {
		fprintf(stderr,
				"\n\nERROR. Util_Read_Linkage_XML: Failed to open the broker linkage config file %s\nThis may be because the file does not exist or it may not be valid XML.\n",
				fileName);
		exit(-1);

	}

	Util_XML_Get_Value_String(fileName, ATLANTIS_ATTRIBUTE, TRUE, inputDoc->children, "URL", url);
	bm->linkageInterface = (LinkageInterface *) malloc(sizeof(LinkageInterface));
	Linkage_Initialise(bm, url);

	/* Read in each of the functional groups */
	polygonList = Util_XML_Get_Node_List(ATLANTIS_ATTRIBUTE_SUB_GROUP, inputDoc->children, "Polygon");

	if (polygonList && polygonList->nodesetval->nodeNr > 0) {

		/* Check the number of groups match the values in the run input file*/
		numPolygons = polygonList->nodesetval->nodeNr;
		printf("numPolygons = %d\n", numPolygons);

		Linkage_SetNumPolygons(bm, numPolygons);

		/* Now read in the groups. */
		for (polygonIndex = 0; polygonIndex < numPolygons; polygonIndex++) {

			polygonNode = polygonList->nodesetval->nodeTab[polygonIndex];
			Util_XML_Get_Value_Integer(fileName, ATLANTIS_ATTRIBUTE, 0, TRUE, polygonNode, integer_check, "id",
					&bm->linkageInterface->polygonList[polygonIndex]);

		}

		xmlXPathFreeObject(polygonList);

	} else {
		quit("ERROR: Error in the broker linkage config XML input file %s\n", fileName);
	}

	xmlFreeDoc(inputDoc);
	/* Shutdown libxml */
	xmlCleanupParser();
	return TRUE;
}

void Linkage_Start(MSEBoxModel *bm) {
	Message *message = NULL;

	while (bm->linkageInterface->shutdown == FALSE) {
		message = Linkage_ReceiveCommand(bm->linkageInterface);
		message = Linkage_DeserialiseCommand(bm, bm->linkageInterface, message);
		if (message != NULL)
			Linkage_SendCommand(bm->linkageInterface, message);
	}
}

void Linkage_End(MSEBoxModel *bm) {
	Linkage_Free(bm);
}
