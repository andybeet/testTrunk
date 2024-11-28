/*
 * LinkageInterface.c
 *
 *  Created on: Jun 7, 2012
 *      Author: bec
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "atBrokerLink.h"

#include <sjwlib.h>
#include <netcdf.h>
#include <atlantisboxmodel.h>
#include <atUtilLib.h>
#include "LinkageInterface.h"
#include "Deserialiser.h"

static void Linkage_Free_Message(Message *message);

int Linkage_Initialise(MSEBoxModel *bm, char *networkURL) {

	bm->linkageInterface->networkInterface = (Network *) malloc(sizeof(Network));
	Network_Initialise(bm->linkageInterface->networkInterface, networkURL);
	bm->linkageInterface->group = -1;
	bm->linkageInterface->polygonList = NULL;
	bm->linkageInterface->shutdown = FALSE;
	bm->linkageInterface->verbose = 1;
	bm->linkageInterface->valueList = NULL;

	/* Linkage arrays */
	bm->linkageInterface->linkageWCDetritusFlux = Util_Alloc_Init_3D_Double(bm->K_num_tot_sp, bm->wcnz, bm->nbox, 0.0);
	bm->linkageInterface->linkageSEDDetritusFlux = Util_Alloc_Init_3D_Double(bm->K_num_tot_sp, bm->sednz, bm->nbox, 0.0);

	bm->linkageInterface->linkageWCMortality = Util_Alloc_Init_4D_Double(bm->K_num_max_cohort * bm->K_num_max_genetypes, bm->K_num_tot_sp, bm->wcnz, bm->nbox, 0.0);
	bm->linkageInterface->linkageSEDMortality = Util_Alloc_Init_4D_Double(bm->K_num_max_cohort * bm->K_num_max_genetypes, bm->K_num_tot_sp, bm->sednz, bm->nbox, 0.0);
	bm->linkageInterface->linkageEPIMortality = Util_Alloc_Init_3D_Double(bm->K_num_max_cohort * bm->K_num_max_genetypes, bm->K_num_tot_sp, bm->nbox, 0.0);



	return 0;
}

void Linkage_Free(MSEBoxModel *bm){
	Network_Close(bm->linkageInterface->networkInterface);
	if (bm->linkageInterface->networkInterface)
		free(bm->linkageInterface->networkInterface);

	if (bm->linkageInterface->polygonList != NULL) {
		free(bm->linkageInterface->polygonList);
	}
	if (bm->linkageInterface->valueList != NULL)
		free(bm->linkageInterface->valueList);

	/* Free the linkage arrays */
	free3d(bm->linkageInterface->linkageWCDetritusFlux);
	free3d(bm->linkageInterface->linkageSEDDetritusFlux);
	free4d(bm->linkageInterface->linkageWCMortality);
	free4d(bm->linkageInterface->linkageSEDMortality);
	free3d(bm->linkageInterface->linkageEPIMortality);

	free(bm->linkageInterface);

}
int Linkage_SetGroup(LinkageInterface *interface, int group) {

	if (interface == NULL) {
		fprintf(stderr, "Linkage_SetGroup. Interface is Null\n");
		return -1;
	}
	interface->group = group;
	return 0;
}

int Linkage_SetNumPolygons(MSEBoxModel *bm, int numPolygons){

	if (bm->linkageInterface == NULL) {
		fprintf(stderr, "Linkage_SetPolygons. Interface is Null\n");
		return -1;
	}
	bm->linkageInterface->numPolygons = numPolygons;
	bm->linkageInterface->polygonList = (int *) malloc(sizeof(int) * numPolygons);

	/* Now allocate the valuesList array as well */
	if (bm->linkageInterface->valueList != NULL)
		free(bm->linkageInterface->valueList);
	bm->linkageInterface->valueList = (double *) malloc(sizeof(double) * numPolygons);

	return 0;
}

//void Linkage_Free(LinkageInterface *interface) {
//
//	Network_Close(interface->networkInterface);
//	if (interface->networkInterface)
//		free(interface->networkInterface);
//
//	if (interface->polygonList != NULL) {
//		free(interface->polygonList);
//	}
//	if (interface->valueList != NULL)
//		free(interface->valueList);
//}

int Linkage_SendCommand(LinkageInterface *interface, Message *message) {
	int sucess;

	sucess = Network_Send(interface->networkInterface, message);

	/* Now free up the message */
	Linkage_Free_Message(message);
	return sucess;
}

void Linkage_Free_Message(Message *message) {

	if (message->bytes != NULL)
		free(message->bytes);
	free(message);
}

Message *Linkage_ReceiveCommand(LinkageInterface *interface) {
	Message *response;

	response = Network_Receive(interface->networkInterface);

	return response;
}

Message *Linkage_DeserialiseCommand(MSEBoxModel *bm, LinkageInterface *interface, Message *message) {
	Message *responseMessage;
	responseMessage = Deserialise_Message(bm,interface, message);

	/* Now free up the message */
	Linkage_Free_Message(message);

	return responseMessage;
}

