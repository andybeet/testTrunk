/**
 * 	\ingroup atImplementationLib
 *	\file atImplementationParamIO.c
 *	\brief Functions to load the management parameters from the manage.xml input document.
 *	Created: 12-08-06-2009 Bec Gorton CSIRO 2009
 *
 *
 *
 */

/*  Open library routines */
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <atImplementation.h>
#include <atManageLib.h>


void readManagementImplementationChangeXML(MSEBoxModel *bm, char *fileName, xmlNodePtr rootnode) {

	xmlNodePtr groupingNode;
	int TAC_max_num_changes;
	groupingNode = Util_XML_Get_Node(ATLANTIS_ATTRIBUTE_SUB_GROUP, rootnode, "Fishing_Catchability_Change");
	if(groupingNode == NULL){
		quit("readManagementImplementationChangeXML node 'Fishing_Catchability_Change' not found.\n");
	}

	Util_XML_Read_Impacted_Group_Param(bm, fileName, groupingNode, TAC_num_changes_id);
	TAC_max_num_changes = Manage_Get_Max_Fishery_Param(bm, TAC_num_changes_id ) + 1;


	/** Allocate the storage to hold the information on the TAC change */
	TACchange = Util_Alloc_Init_4D_Double ( 3, TAC_max_num_changes, bm->K_num_fisheries, bm->K_num_tot_sp, 0.0 );

	/* Now load the data */
	Util_XML_Read_Change_Fished_Groups(bm, fileName, groupingNode, "TACchange", TAC_num_changes_id, TACchange);
}


/**
 *	This routine reads all parameters of the management implementation model
 *
 *  Input:  char *filename: The name of the management implementation xml input file.
 *
 */
int readModelManageImplParameters(MSEBoxModel *bm, char *filename) {
	xmlDocPtr inputDoc;
	inputDoc = xmlReadFileDestFolder(bm->destFolder, filename, NULL, 0);

	readManagementImplementationChangeXML(bm, filename, inputDoc->children);
	xmlFreeDoc(inputDoc);
	/* Shutdown libxml */
	xmlCleanupParser();
	return TRUE;
}
