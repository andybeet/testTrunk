/*
 * network.c
 *
 *  Created on: Jun 6, 2012
 *      Author: bec
 */

#include <zmq.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include "Message.h"
#include "network.h"
#include "NetworkError.h"

void Network_Initialise(Network *network, char *modelurl) {
	network->modelURL = strdup(modelurl);
	network->context = zmq_init(1);
	network->socket = zmq_socket(network->context, ZMQ_REP);
	zmq_bind(network->socket, network->modelURL);
}

Message *Network_Receive(Network *network) {

	int result;
	zmq_msg_t request;
	Message *messageReceived;

	zmq_msg_init(&request);
	result = zmq_recv(network->socket, &request, 0);
	if (result == -1) {
		printErrorMessage("Network_Receive()");
	}

	messageReceived = (Message *) malloc(sizeof(Message));
	messageReceived->size = zmq_msg_size(&request);
	messageReceived->bytes = malloc(sizeof(void) * (messageReceived->size));
	memcpy(messageReceived->bytes, zmq_msg_data(&request),
			messageReceived->size);

	zmq_msg_close(&request);

	return messageReceived;
}
int Network_Send(Network *network, Message *message) {

	zmq_msg_t response;
	int result;

	/* Send the data */
	zmq_msg_init_size(&response, message->size);
	memcpy(zmq_msg_data(&response), message->bytes, message->size);
	result = zmq_send(network->socket, &response, 0);

	/* Check to see if we sent successfully */
	if (result == -1) {
		printf("error message\n");
		printErrorMessage("Network_Send()");
		return -1;
	}
	zmq_msg_close(&response);

	return 0;
}

void Network_Close(Network *network) {
	if(network->modelURL != NULL)
		free(network->modelURL);
	zmq_close(network->socket);
	zmq_term(network->context);
}
