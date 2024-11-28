/**

 \file
 \brief This file contains functions to convert the harvest prm file to the implemention.XML file. Only the management implementation parameters will be pulled
 into the new XML file.

 Changes:


 */
/*  Open library routines */
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdarg.h>
#include <math.h>
#include <ctype.h>
#include <assert.h>
#include <netcdf.h>
#include <sjwlib.h>
#include <atlantisboxmodel.h>
#include <atUtilLib.h>
#include <atImplementationLib.h>
#include "convertXML.h"

void ManagementImplementationChangeXML(MSEBoxModel *bm, FILE *fp, char *fileName, xmlDocPtr doc, xmlNodePtr rootnode) {

	xmlNodePtr groupingNode;
	groupingNode = Util_XML_Create_Node(ATLANTIS_ATTRIBUTE_SUB_GROUP, rootnode, "Fishing_Catchability_Change", "Changes in catchability through time", "", "");

	Create_Harvest_Fishery_Group_ParamXML(bm, fp, fileName, groupingNode, TAC_num_changes_id, "Number of changes of fishing catchability.", "", XML_TYPE_INTEGERARRAY);

	Create_Fished_Groups_ChangeXML(bm, fp, fileName, doc, groupingNode, "TACchange", "catchability");
}

void Convert_Implementation_To_XML(MSEBoxModel *bm, char *fileName, char *outputFileName) {
	xmlDocPtr doc;
	xmlNodePtr rootnode;
	FILE *fp;

	/* Init libxml */
	xmlInitParser();
	LIBXML_TEST_VERSION

	if ((fp = Open_Input_File(bm->inputFolder,fileName, "rt")) == NULL) {
		quit("Cannot open implementation input file %s%s\n", bm->inputFolder, fileName);
	}

	doc = xmlNewDoc((unsigned char *) "1.0");

	if (doc == NULL) {
		fprintf(stderr, "Document not created successfully. \n");
		return;
	}

	rootnode = xmlNewDocNode(doc, NULL, (const xmlChar*) "AtlantisManagementImplementation", NULL);
	xmlDocSetRootElement(doc, rootnode);

	ManagementImplementationChangeXML(bm, fp, fileName, doc, rootnode);

	/* Save to the output file */
	xmlSaveFormatFileDestFolder(bm->destFolder, outputFileName, doc, 1);
	xmlFreeDoc(doc);
	fclose(fp);

	/* Shutdown libxml */
	xmlCleanupParser();

	return;
}

