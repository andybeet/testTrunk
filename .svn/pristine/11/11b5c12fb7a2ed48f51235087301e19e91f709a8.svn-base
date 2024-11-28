/*
 * atSS3LinkLib.h
 *
 *  Created on: May 2, 2013
 *      Author: bec
 */

#ifndef ATSS3LINKLIB_H_
#define ATSS3LINKLIB_H_

int SS3Link_Init(MSEBoxModel *bm, int species);

void SS3Link_Free(MSEBoxModel *bm);
//void Start_SS3_Link(MSEBoxModel *bm, char *inputFileName);

void WriteSSFiles(MSEBoxModel *bm, int species, int year, FILE *llogfp);
void WriteHistory(MSEBoxModel *bm, int groupIndex);

void FindWhichData(MSEBoxModel *bm, int groupIndex, int iyr);
void GenData(MSEBoxModel *bm, int groupIndex, int yearIndex);
void GetLengthWeightData(MSEBoxModel *bm, int groupIndex, int iyr);

#endif /* ATSS3LINKLIB_H_ */
