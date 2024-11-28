/**
\file
\brief Header file for the Atlantis-EwE communication functions.

    File:		atComms.c
	Created:	4/9/2008
	Author:		Bec Gorton,   CSIRO Marine Laboratories,  2008
	Purpose:	Header file for the Atlantis-EwE communication functions.
    Note:

    Revisions:

*************/

#ifndef ATCOMMS_H_
#define ATCOMMS_H_

/* The set functions */
int setBiomass(char *species, double **dataArray, int numRows, int numColumns);
int setProduction(char *species, double **dataArray, int numRows, int numColumns);
int setMortality(char *species, double **dataArray, int numRows, int numColumns);
int setNutrients(char *species, double **dataArray, int numRows, int numColumns);
int setProperties(char *species, double **dataArray, int numRows, int numColumns);

/* The get functions */
int getBiomass(char *species, double **values, int numRows, int numColumns);
int getProduction(char *species, double **values, int numRows, int numColumns);
int getMortality(char *species, double **values, int numRows, int numColumns);
int getNutrients(char *species, double **values, int numRows, int numColumns);
int getProperties(char *species, double **values, int numRows, int numColumns);

/* The function to actually start the EwE model */
int runModel(void);
int setupEwECommunication(void);
int TestLocal(void);
void sayHello(int argc, char **argv);
void sayEwEHello(int argc, char **argv);

#endif /*ATCOMMS_H_*/
