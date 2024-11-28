/*
 * atBrokerLink.h
 *
 *  Created on: Jun 25, 2012
 *      Author: bec
 */

#ifndef ATBROKERLINK_H_
#define ATBROKERLINK_H_

#include "requests.pb-c.h"
#include "responses.pb-c.h"
#include <zmq.h>
#include "Message.h"
#include "network.h"

typedef struct{
	Network *networkInterface;
	unsigned int numPolygons;
	int *polygonList;
	unsigned int group;
	char groupName[1024];
	int shutdown;
	int verbose;
	double *valueList;

	/*
	 * The additional mortality values that are used when atlantis links to other models such as EwE
	 */
	double ****linkageWCMortality;
	double ****linkageSEDMortality;
	double ***linkageEPIMortality;

	/* The additional detritus flux values that are used to change the atlantis detritus biomass values
	 * when Atlantis links with other models such as EwE
	 */

	double ***linkageWCDetritusFlux;
	double ***linkageSEDDetritusFlux;

}LinkageInterface;



#endif /* ATBROKERLINK_H_ */
