/*
 * atSS3LinkPrivate.h
 *
 *  Created on: May 2, 2013
 *      Author: bec
 */

#ifndef ATSS3LINKPRIVATE_H_
#define ATSS3LINKPRIVATE_H_


typedef struct{
	char SS3Path[STRLEN];
	char starterFileName[STRLEN];
	char dataFileName[STRLEN];
	char controlFileName[STRLEN];
	int doRun;

} SS3LinkStructure;


double Read_SS3_Report_File(MSEBoxModel *bm, int species, char *folderName);
double Read_SS3_Par_File(MSEBoxModel *bm, char *folderName);

void Create_Starter_File(MSEBoxModel *bm, char *dirName, int groupIndex, int versionID);
void Write_Forecast_File(MSEBoxModel *bm, char *dirName, int maxyr, int groupIndex, int versionID);
void Write_SS_Data_File(MSEBoxModel *bm, char *dirName, char *fileName, int numYears, int groupIndex, int versionID);
void Write_SS_Control_File(MSEBoxModel *bm, char *dirName, char *fileName, int numYears, int groupIndex, int versionID);
void SS3_Data_Gen(MSEBoxModel *bm);
void SS3_Find_Hist_Data(MSEBoxModel *bm, int groupIndex);

int SS3_Test_Init(MSEBoxModel *bm);

#endif /* ATSS3LINKPRIVATE_H_ */
