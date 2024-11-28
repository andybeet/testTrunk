/*
 * network.h
 *
 *  Created on: Jun 6, 2012
 *      Author: bec
 */

#ifndef NETWORK_H_
#define NETWORK_H_

typedef struct
{
	/**
	 * the URL of the model server to connect to; eg tcp://localhost:555
	 */
	char *modelURL;
	void *context;
	void *socket;
}Network;


void Close_Message(Network *network);

/**
 * This method is called once during broker initialisation, and may be used
 * for any socket initialisation, etc. It is an error to call this more than
 * once.
 */
void Network_Initialise(Network *network, char *modelurl);

/**
 * Used to send binary data over this connection; results are also returned
 * as a byte array. This is a blocking method. It is an error to call this
 * before calling {@link #init()} or after {@link #close()}.
 *
 * @param data
 *            The data to transmit. Should not be null.
 * @return The response data; never null, but may be zero-length.
 */
int Network_Send(Network *network, Message *message);

Message *Network_Receive(Network *network);
/**
 * This method is called once at the end of the broker's life-cycle. It may
 * be used to tear down connections, release resources, etc. It is an error
 * to call this more than once, and undefined what happens if
 * {@link #init()} is called after calling {@link #close()}.
 */
void Network_Close(Network *network);


#endif /* NETWORK_H_ */
