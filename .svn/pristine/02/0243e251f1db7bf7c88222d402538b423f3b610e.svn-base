/*
 * Deserialiser.c
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
#include <ImportExportData.h>

/* Generic functions */
static Message *Serialise_Boolean_Response(int responseValue);
static Message *Serialise_Data_Response(double *data, int length);

/* Init requests */
static double Deserialiser_Init(Requests__Request *request);

/* RunNextTimeStep requests */
static void Deserialiser_RunNextTimestep(Requests__Request *request);

/* Shutdown */
static void Deserialiser_ShutDown(Requests__Request *request);
/* GetBiomass*/
static void Deserialiser_GetBiomass(LinkageInterface *interface, Requests__Request *request);
static int Deserialiser_SetBiomass(LinkageInterface *interface, Requests__Request *request);
/* GetMortality*/
static void Deserialiser_GetMortality(LinkageInterface *interface, Requests__Request *request);
static int Deserialiser_SetMortality(LinkageInterface *interface, Requests__Request *baseRequest);
/* GetDetritus */
static void Deserialiser_GetDetritus(LinkageInterface *interface, Requests__Request *request);
static int Deserialiser_SetDetritus(LinkageInterface *interface, Requests__Request *baseRequest);

Message *Deserialise_Message(MSEBoxModel *bm, LinkageInterface *interface, Message *message) {
	Requests__Request *request;
	double time;
	Message *response = NULL;
	double *data;
	unsigned int i;
	int halt;

	// Unpack the message using protobuf-c.
	request = requests__request__unpack(NULL, message->size, message->bytes);
	if (request == NULL) {
		fprintf(stderr, "Deserialiser_Init - error unpacking incoming message\n");
		exit(1);
	}

	switch (request->type) {
	case REQUESTS__REQUEST__TYPE__INIT:

		time = Deserialiser_Init(request);
		printf("time = %e\n", time);
		if (interface->verbose > 0)
			printf("Init request, time - %f\n", time);
		Link_SetTime(bm, time);
		response = Serialise_Boolean_Response(TRUE);
		break;
	case REQUESTS__REQUEST__TYPE__RUNNEXTTIMESTEP:

		if (interface->verbose > 0)
			printf("RunNextTimeStep request\n");
		/* Model run next timestep */
		Deserialiser_RunNextTimestep(request);

		halt = Link_RunNextTimeStep(bm);
		printf("halt = %d\n", halt);

		/* Run model timestep */
		response = Serialise_Boolean_Response(TRUE);
		break;
	case REQUESTS__REQUEST__TYPE__SHUTDOWN:
		if (interface->verbose > 0)
			printf("ShutDown request\n");
		Deserialiser_ShutDown(request);
		interface->shutdown = TRUE;
		/* No response required - but broker sits there listening for one.*/
		response = NULL;
		response = Serialise_Boolean_Response(TRUE);
		break;
	case REQUESTS__REQUEST__TYPE__GETBIOMASS:
		Deserialiser_GetBiomass(interface, request);

		/* Get the response*/
		data = Util_Alloc_Init_1D_Double(interface->numPolygons, 0.0);
		Link_GetBiomass(bm, interface->groupName, interface->numPolygons, interface->polygonList, data);
		//data = Link_GetBiomass(bm, interface->groupName, interface->numPolygons, interface->polygonList);

		for (i = 0; i < interface->numPolygons; i++) {
			printf("data[%d] = %e\n", i, data[i]);
		}

		response = Serialise_Data_Response(data, interface->numPolygons);
		free(data);
		break;
	case REQUESTS__REQUEST__TYPE__GETMORTALITY:
		Deserialiser_GetMortality(interface, request);

		/* Get the response*/
		data = Util_Alloc_Init_1D_Double(interface->numPolygons, 0.0);
		for (i = 0; i < interface->numPolygons; i++) {
			data[i] = (double) i;
		}
		response = Serialise_Data_Response(data, interface->numPolygons);
		free(data);
		break;
	case REQUESTS__REQUEST__TYPE__GETDETRITUS:
		Deserialiser_GetDetritus(interface, request);

		/* Get the response*/

		data = Util_Alloc_Init_1D_Double(interface->numPolygons, 0.0);
		Link_GetBiomass(bm, interface->groupName, interface->numPolygons, interface->polygonList, data);

		for (i = 0; i < interface->numPolygons; i++) {
			printf("data[%d] = %e\n", i, data[i]);
		}

		response = Serialise_Data_Response(data, interface->numPolygons);
		free(data);
		break;
	case REQUESTS__REQUEST__TYPE__SETBIOMASS:
		Deserialiser_SetBiomass(interface, request);

		Link_SetBiomass(bm, interface->groupName, interface->numPolygons, interface->polygonList, interface->valueList);

		/* Run model timestep */
		response = Serialise_Boolean_Response(TRUE);

		break;
	case REQUESTS__REQUEST__TYPE__SETMORTALITY:
		Deserialiser_SetMortality(interface, request);
		/* Run model timestep */
		response = Serialise_Boolean_Response(TRUE);

		break;
	case REQUESTS__REQUEST__TYPE__SETDETRITUS:
		Deserialiser_SetDetritus(interface, request);

		//Link_SetBiomass(bm, interface->groupName, interface->numPolygons, interface->polygonList, interface->valueList);

		/* Run model timestep */
		response = Serialise_Boolean_Response(TRUE);
		break;
	}
	return response;
}

/**
 * Return t.
 */double Deserialiser_Init(Requests__Request *request) {

	Requests__Init *requestInit;
	double returnValue;

	requestInit = request->init;

	// display the message's fields.
	returnValue = atof(requestInit->timestep);

	// Free the unpacked message
	requests__request__free_unpacked(request, NULL);
	return returnValue;

}

Message *Serialise_Boolean_Response(int responseValue) {
	Message *message;
	Responses__Response response = RESPONSES__RESPONSE__INIT; // AMessage
	Responses__BooleanResponse responseInit = RESPONSES__BOOLEAN_RESPONSE__INIT;

	message = (Message *) malloc(sizeof(Message));

	response.bool_ = &responseInit;
	response.type = RESPONSES__RESPONSE__TYPE__BOOL;
	responseInit.response = responseValue;

	message->size = responses__response__get_packed_size(&response);
	message->bytes = malloc(sizeof(void) * (message->size));

	responses__response__pack(&response, message->bytes);

	return message;
}

Message *Serialise_Data_Response(double *data, int length) {
	Message *message;
	Responses__Response response = RESPONSES__RESPONSE__INIT; // AMessage
	Responses__GetResponse responseGet = RESPONSES__GET_RESPONSE__INIT;

	message = (Message *) malloc(sizeof(Message));

	response.get = &responseGet;
	response.type = RESPONSES__RESPONSE__TYPE__GET;

	/* Allocate memory and copy data across */
	responseGet.n_values = length;
	responseGet.values = data;

	message->size = responses__response__get_packed_size(&response);
	message->bytes = malloc(sizeof(void) * (message->size));

	responses__response__pack(&response, message->bytes);

	return message;
}

/**
 * Simply unpacking incase there is an error.
 */
void Deserialiser_RunNextTimestep(Requests__Request *request) {

	Requests__RunNextTimeStep *requestNextTimeStep;

	requestNextTimeStep = request->runnext;
	if (requestNextTimeStep == NULL) {
		fprintf(stderr, "Deserialiser_RunNextTimestep: requestNextTimeStep message is null.\n");
	}

	// Free the unpacked message
	requests__request__free_unpacked(request, NULL);
}

/**
 * Simply unpacking incase there is an error.
 */
void Deserialiser_ShutDown(Requests__Request *request) {

	Requests__Shutdown *shutDown;

	shutDown = request->shutdown;
	if (shutDown == NULL) {
		fprintf(stderr, "Deserialiser_ShutDown: shutDown message is null.\n");
	}

	// Free the unpacked message
	requests__request__free_unpacked(request, NULL);
}

/**
 * Unpack the getBiomass request..
 */
void Deserialiser_GetBiomass(LinkageInterface *interface, Requests__Request *request) {

	Requests__GetBiomass *biomassRequest;
	char *group;

	biomassRequest = request->getbiomass;
	if (biomassRequest == NULL) {
		fprintf(stderr, "Deserialiser_GetBiomass: getbiomass message is null.\n");
	}
	group = biomassRequest->groupid;
	strcpy(interface->groupName, group);
	if (interface->verbose) {
		printf("Deserialiser_GetBiomass. Group %s\n", interface->groupName);
	}

	/* Free the unpacked message */
	requests__request__free_unpacked(request, NULL);
}

/**
 * Unpack the setbiomass request..
 */
int Deserialiser_SetBiomass(LinkageInterface *interface, Requests__Request *baseRequest) {

	Requests__SetBiomass *request;

	size_t length;
	unsigned int i;

	request = baseRequest->setbiomass;
	if (request == NULL) {
		fprintf(stderr, "Deserialiser_SetBiomass: setbiomass message is null.\n");
	}

	/* Check the length specified is what we expect */
	length = request->n_values;
	if (length != interface->numPolygons) {
		fprintf(stderr, "Deserialiser_SetBiomass. Number of values sent from Broker %ld does not match number of expected polygon values %d\n", length,
				interface->numPolygons);
		exit(-1);
	}
	/* Copy the groupName and valueList data across */
	strcpy(interface->groupName, request->groupid);
	memcpy(interface->valueList, request->values, (sizeof(double) * length));

	if (interface->verbose) {
		printf("SetBiomass - %s, [", interface->groupName);
		for (i = 0; i < length; i++) {
			printf("%s%e", i > 0 ? ", " : "", interface->valueList[i]);
		}
		printf("]\n");
	}

	/* Free the unpacked message */
	requests__request__free_unpacked(baseRequest, NULL);
	return TRUE;
}

/**
 * Unpack the getMortality request..
 */
void Deserialiser_GetMortality(LinkageInterface *interface, Requests__Request *baseRequest) {

	Requests__GetMortality *request;

	request = baseRequest->getmortality;
	if (request == NULL) {
		fprintf(stderr, "Deserialiser_GetBiomass: getMortality message is null.\n");
	}

	strcpy(interface->groupName, request->groupid);
	if (interface->verbose) {
		printf("Deserialiser_GetMortality. Group %s\n", interface->groupName);
	}

	/* Free the unpacked message */
	requests__request__free_unpacked(baseRequest, NULL);
}

/**
 * Unpack the setmortality request..
 */
int Deserialiser_SetMortality(LinkageInterface *interface, Requests__Request *baseRequest) {

	Requests__SetMortality *request;
	size_t length;
	unsigned int i;

	request = baseRequest->setmortality;
	if (request == NULL) {
		fprintf(stderr, "Deserialiser_SetMortality: setmortality message is null.\n");
	}

	/* Check the length specified is what we expect */
	length = request->n_values;
	if (length != interface->numPolygons) {
		fprintf(stderr, "Deserialiser_SetMortality. Number of values sent from Broker does not match number of expected polygon values\n");
		exit(-1);
	}
	/* Copy the groupName and valueList data across */
	strcpy(interface->groupName, request->groupid);
	memcpy(interface->valueList, request->values, (sizeof(double) * length));

	if (interface->verbose) {
		printf("SetMortality - %s, [", interface->groupName);
		for (i = 0; i < length; i++) {
			printf("%s%e", i > 0 ? ", " : "", interface->valueList[i]);
		}
		printf("]\n");
	}
	/* Free the unpacked message */
	requests__request__free_unpacked(baseRequest, NULL);
	return TRUE;
}

/**
 * Unpack the getMortality request..
 */
void Deserialiser_GetDetritus(LinkageInterface *interface, Requests__Request *baseRequest) {

	Requests__GetDetritus *request;

	request = baseRequest->getdetritus;
	if (request == NULL) {
		fprintf(stderr, "Deserialiser_GetDetritus: getdetritus message is null.\n");
	}

	strcpy(interface->groupName, request->groupid);
	if (interface->verbose) {
		printf("Deserialiser_GetDetritus. Group %s\n", interface->groupName);
	}

	/* Free the unpacked message */
	requests__request__free_unpacked(baseRequest, NULL);
}

/**
 * Unpack the setdetritus request..
 */
int Deserialiser_SetDetritus(LinkageInterface *interface, Requests__Request *baseRequest) {

	Requests__SetDetritus *request;
	size_t length;
	unsigned int i;

	request = baseRequest->setdetritus;
	if (request == NULL) {
		fprintf(stderr, "Deserialiser_SetDetritus: setdetritus message is null.\n");
	}

	/* Check the length specified is what we expect */
	length = request->n_values;
	if (length != interface->numPolygons) {
		fprintf(stderr, "Deserialiser_SetDetritus. Number of values sent from Broker does not match number of expected polygon values\n");
		exit(-1);
	}
	/* Copy the groupName and valueList data across */
	strcpy(interface->groupName, request->groupid);
	memcpy(interface->valueList, request->values, (sizeof(double) * length));

	if (interface->verbose) {
		printf("SetDetritus - %s, [", interface->groupName);
		for (i = 0; i < length; i++) {
			printf("%s%e", i > 0 ? ", " : "", interface->valueList[i]);
		}
		printf("]\n");
	}
	/* Free the unpacked message */
	requests__request__free_unpacked(baseRequest, NULL);
	return TRUE;
}

