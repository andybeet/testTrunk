/*******************************************************************//**
\file
\brief Header file for the physics module.

    File:           atPhysics.h

    Created:		19/12/2008

    Author:         Bec Gorton
                    CSIRO

    Purpose:       	Header file for the physics module.

    Arguments:      See below

    Returns:        void

    Revisions:      24-07-2009 Bec Gorton
					Added the atlantisUtil header file to this file so that we can use the array functions defined in the atlantisUtil lib.

*********************************************************************/

/* Functions to deal with the allocation and freeing of the atPhysicsStructure */
void allocatePhysicsModule(MSEBoxModel *bm);
void freePhysicsStruct(MSEBoxModel *bm);

void freeTempSalt(MSEBoxModel *bm);

int Load_Ice_Timeseries(MSEBoxModel *bm);
int CDF_iceBM(MSEBoxModel *bm);

void init_PhyPropertyData(MSEBoxModel *bm, FILE *fp, PhyPropertyData *propInput, char *variableName, char *shortName, char *longName, double minValue, double maxValue, int is_valid_z);
void get_property(MSEBoxModel *bm, PhyPropertyData *propInput);
