/**
 *	\file The C file containing the functions that allocation and free the atAssessStructure.
 *
 *
 *
 *	@author Bec Gorton 28/Nov/2008
 */

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <math.h>
#include <sjwlib.h>
#include <atlantisboxmodel.h>
#include <atUtilLib.h>

void allocatePhysicsModule(MSEBoxModel *bm){

	bm->atPhysicsModule = (atPhysicsStructure *)malloc(sizeof(atPhysicsStructure));
	bm->atPhysicsModule->inpfp = NULL;
	bm->atPhysicsModule->masstosed = NULL;
	bm->atPhysicsModule->expfp = NULL;
	bm->atPhysicsModule->totinp = NULL;
}

void closePhysicsFile(FILE *fp)
{
	if(fp != NULL)
		fclose(fp);

}
void freePhysicsStruct(MSEBoxModel *bm){

	closePhysicsFile(bm->atPhysicsModule->inpfp);
	closePhysicsFile(bm->atPhysicsModule->expfp);

	if(bm->atPhysicsModule->masstosed != NULL)
		free2d(bm->atPhysicsModule->masstosed);

	if(bm->atPhysicsModule->totinp != NULL)
		free1d(bm->atPhysicsModule->totinp);

	free(bm->atPhysicsModule);
}
