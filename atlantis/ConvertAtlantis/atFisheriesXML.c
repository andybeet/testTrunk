/*
 * atFisheriesXML.c
 *
 *  Created on: Dec 9, 2011
 *      Author: bec
 */


/*
 * atGroupsXML.c
 *
 *  Created on: 23/08/2011
 *      Author: bec
 */

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
#include "convertXML.h"

/*
<AttributeSubGroup AttributeGroupName="Fishery" AttributeDescription="The defintion of a fishery">
<Attribute AttributeName="Code" AttributeValue="midwcFP"/>
<Attribute AttributeName="Index" AttributeValue="3"/>
<Attribute AttributeName="Name" AttributeValue="midwater trawl on planktivores"/>
<Attribute AttributeName="IsRec" AttributeValue="0"/>
<Attribute AttributeName="NumSubFleets" AttributeValue="0"/>
</AttributeSubGroup>
*/

void Convert_Fisheries_To_XML(MSEBoxModel *bm, char *fileName, char *outputFileName) {
	xmlDocPtr doc;
	xmlNodePtr rootnode;
	xmlNodePtr groupNode;
	FILE *inputFP;
	int numExpectedTokens = 5;
	int buflen = 2000;
	char ch, buf[2000], seps[] = ",";
	char *varStr;
	int count = 0;

	/* Init libxml */
	xmlInitParser();
	LIBXML_TEST_VERSION

	if ((inputFP = Open_Input_File(bm->inputFolder,fileName, "rt")) == NULL) {
		quit("Cannot open Atlantis fisheries definition input file %s%s\n", bm->inputFolder, fileName);
	}

	doc = xmlNewDoc((unsigned char *) "1.0");

	if (doc == NULL) {
		fprintf(stderr, "Document not created successfully. \n");
		return;
	}

	rootnode = xmlNewDocNode(doc, NULL, (const xmlChar*) AtlantisXMLObjectNAMES[ATLANTIS_ATTRIBUTE_GROUP], NULL);
	Util_XML_Add_Node_Property(rootnode, AtlantisXMLObjectNAMES[ATLANTIS_ATTRIBUTE_GROUP][ATTRIBUTE_NAME], "AtlantisFisheries");

	xmlDocSetRootElement(doc, rootnode);

	/* Now parse each group and create the XML */

	fseek(inputFP, 0L, 0);

	/* Get the first line */
	fgets(buf, buflen, inputFP);

	/* First we need to check that the number of tokens is correct */
	varStr = strtok(buf, seps);
	while (varStr != NULL) {
		varStr = strtok(NULL, seps);
		count++;
	}

	if (count != numExpectedTokens) {
		quit("\n\nERROR: The fisheries input file does not have the correct number of columns %d it has %d columns\n\n", numExpectedTokens, count);

	}

	count = 0;
	while (fgets(buf, buflen, inputFP) != NULL) {

		/* Ignore the first line */
		if (strstr(buf, "Code") != NULL)
			continue;

		ch = buf[0];

		if (ch != '\t' && ch != '\n' && ch != '#') {

			/* Create the node for the group */
			groupNode = Util_XML_Create_Node(ATLANTIS_ATTRIBUTE_SUB_GROUP, rootnode, "Fishery", "The definition of a fishery", "", "");

			/* Group Code value. */
			varStr = strtok(buf, seps);
			Util_XML_Create_Node(ATLANTIS_ATTRIBUTE, groupNode, "Code", "", "", varStr);

			/* Get the index. - ignore this value */
			varStr = strtok(NULL, seps);
			Util_XML_Create_Node(ATLANTIS_ATTRIBUTE, groupNode, "Index", "", "", varStr);

			/* Fishery Name */
			varStr = strtok(NULL, seps);
			Util_XML_Create_Node(ATLANTIS_ATTRIBUTE, groupNode, "Name", "", "", varStr);

			/* IsRec */
			varStr = strtok(NULL, seps);
			Util_XML_Create_Node(ATLANTIS_ATTRIBUTE, groupNode, "IsRec", "", "", varStr);

			/* NumSubFleets */
			varStr = strtok(NULL, seps);
			Util_XML_Create_Node(ATLANTIS_ATTRIBUTE, groupNode, "NumSubFleets", "", "", varStr);

			count++;

		}
	}
	/* Save to the output file */
	xmlSaveFormatFileDestFolder(bm->destFolder, outputFileName, doc, 1);
	xmlFreeDoc(doc);
	fclose(inputFP);

	/* Shutdown libxml */
	xmlCleanupParser();

	return;
}
