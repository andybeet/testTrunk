/*
 * LinkageInterface.h
 *
 *  Created on: Jun 7, 2012
 *      Author: bec
 */

#ifndef LINKAGEINTERFACE_H_
#define LINKAGEINTERFACE_H_


int Linkage_Read_XML(MSEBoxModel *bm, char *fileName);
void Linkage_Start(MSEBoxModel *bm);
void Linkage_Free(MSEBoxModel *bm);
int Linkage_Initialise(MSEBoxModel *bm, char *networkURL);
int Linkage_SetGroup(LinkageInterface *interface, int group);
int Linkage_SetPolygons(MSEBoxModel *bm, int numPolygons, const int *polygons);
int Linkage_SetNumPolygons(MSEBoxModel *bm, int numPolygons);
//void Linkage_Free(LinkageInterface *interface);
void *Linkage_SerialiseCommand(LinkageInterface *interface, const Requests__Request__Type cmd);
int Linkage_SendCommand(LinkageInterface *interface, Message *message);
Message *Linkage_ReceiveCommand(LinkageInterface *interface);
Message *Linkage_DeserialiseCommand(MSEBoxModel *bm, LinkageInterface *interface, Message *message) ;




#endif /* LINKAGEINTERFACE_H_ */
