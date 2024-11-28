/*******************************************************************//**
\file
\brief  Include file for assessment code arrays for Box Model
    File:           atsamplesetup.h

    Created:        18/5/2006

    Author:         Beth Fulton

    Purpose:        Include file for assessment code arrays for Box Model

    Revisions:


    08-04-2010 Bec Gorton
    Removed the CON array as we can now just use the isPredator flag in the functional group definition file.
	Resized the availfood and altavailfood arrays to store stage and habitat information .

	16-04-2010 Bec Gorton
	Added the cvphys array so we can store the coefficient variance values for the physical properties that are
	sampled.

	26-05-2010 Bec Gorton
	Removed the spNetName, SPmappreyid, SPAguildid, SPpreyid and SPfishedid arrays.

*********************************************************************/

/* Parameter arrays */
int    *InOffshorei = 0;

double **detrit_import = 0;
double ****dietpref = 0;
double **DiscardFate = 0;
double *TrigPts = 0;

/* Data storage arrays */
double *samplearea = 0;
double *totalarea = 0;
double *zonearea = 0;
int    *zoneboxes; /* The number of boxes in a zone */

double ***biolbiom = 0;
double **bioleat = 0;
double **biolprod = 0;
double ****biolVERTinfo = 0;
double **phys = 0;
double **physprocess = 0;

/* Sample storage arrays */
double *****individVERTinfo = 0;
double ****invstockinfo = 0;
double ***pop_fraction = 0;
double ***samplebiom = 0;
double ****sampledetbiom;
double ***sampleeat = 0;
double ***samplephys = 0;
double ***sampleprocess = 0;
double ***sampleprod = 0;
double ****stockinfo = 0;

/* Fisheries data arrays */
double *****agebins = 0;
double *****agelengthkey = 0;
double *****fishery = 0;
double **samplenums = 0;
double ******sizebins = 0;

/* Diet info arrays */
double ****inshorediet = 0;
double ****intruediet = 0;
double ****offshorediet = 0;
double ****offtruediet = 0;
double ****totareadiet = 0;
double ****tottruediet = 0;
double ***TL = 0;
double ***trueTL = 0;

/* Data processing arrays */
double ****altavailfood = 0;
double ****availfood = 0;

double **new_altavailfood = 0;
double **new_availfood = 0;

double **contribvert = 0;
double ***biom = 0;
double ***globalnums = 0;
double ***oldbaby = 0;
double **nums = 0;
double *totn = 0;
double ***totnums = 0;
double **totsamplebiom = 0;

/* Network calculation arrays */
double ***biomnetwk = 0;
double ***eatnetwk = 0;
double ***exportnetwk = 0;
double ***importnetwk = 0;
double ***mortnetwk = 0;
double ***prodnetwk = 0;
double ***respnetwk = 0;

/* Indice calculation arrays */
int    **divsp = 0;
int    **divfn = 0;
int    ****nsq = 0;
int    *numsort = 0;
int    *permpivot = 0;
int	   **stomachs = 0;

double *upper = 0;
double *lower = 0;

int *MultiNom_binIndex = 0;
int *MultiNom_cohortIndex = 0;

double ***cvsample = 0;
double ***cvphys = 0; /* Array to store coefficient variance values for the physical parameters that are sampled */
double **globalfledge = 0;
double **globalpd = 0;
double **mineat = 0;
double ***oi = 0;
double **PBRglobal = 0;
double *zasum = 0;

/* Arrays to store the standard deviation of sampled values */
double **physicalSigma;

/* Indices arrays */
double **avgtl = 0;
double *cvt = 0;
double ***endnums = 0;
double ***endnumsbig = 0;
double *disrate = 0;
double **fledge = 0;
double **habindx = 0;
double **max_lngth = 0;
double **monbet = 0;
double **netwkindx = 0;
double *nppb = 0;
double **PBRcat = 0;
double **pd = 0;
double **pdcat = 0;
double ***pelbin = 0;
double **sedbin = 0;
double *soi = 0;
double ***trophspect = 0;

int    **divindx = 0;

/* Asessment arrays */
double *BootResu = 0;
double **CData = 0;
double *CTData = 0;
double **CPUEtrend = 0;
double **F = 0;
double **IData = 0;
double *ITData = 0;
double **IDatahat = 0;
double **NEst = 0;
double ***NResult = 0;
double **num_nyr = 0;
double **P = 0;
double **PP = 0;
double **Resu = 0;
double *ResultToSort = 0;
double *ResultSorted = 0;
double *X = 0;
double *XX = 0;
double *Y = 0;
double *YY = 0;
double ****zoneVERTpopratio = 0;

int *nbsSort = 0;

/* Extra calculation arrays */
int *checkedz = 0;


ErrorStructure *spErrorStructure;

int readModelAssessmentParameters(MSEBoxModel *bm, char *filename);
