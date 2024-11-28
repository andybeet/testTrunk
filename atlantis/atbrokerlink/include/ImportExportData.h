/*
 * ImportExportData.h
 *
 *  Created on: Jun 27, 2012
 *      Author: bec
 */

#ifndef IMPORTEXPORTDATA_H_
#define IMPORTEXPORTDATA_H_

int Link_SetTime(MSEBoxModel *bm, double t);

int Link_RunNextTimeStep(MSEBoxModel *bm);

int Link_SetBiomass(MSEBoxModel *bm, char *groupCode, int numBoxes, int *boxIndexValues, double *BoxValues);
int Link_GetDetritus(MSEBoxModel *bm, char *groupCode, int numBoxes, int *boxIndexValues, double *boxValues);
int Link_GetBiomass(MSEBoxModel *bm, char *groupCode, int numBoxes, int *boxIndexValues, double *boxValues);
int Link_SetMortality(MSEBoxModel *bm, char *groupCode, int numBoxes, int *boxIndexValues, double *boxValues);

#endif /* IMPORTEXPORTDATA_H_ */
