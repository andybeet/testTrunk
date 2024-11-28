/*
 * NetworkError.c
 *
 *  Created on: Jun 6, 2012
 *      Author: bec
 */


#include <zmq.h>
#include <stdio.h>
#include <unistd.h>
#include <string.h>

void printErrorMessage(char *functionName){

	int errorNumber =  zmq_errno();
	fprintf (stderr, "Error occurred during %s: %s\n", functionName, zmq_strerror (errorNumber));
	//abort ();
}
