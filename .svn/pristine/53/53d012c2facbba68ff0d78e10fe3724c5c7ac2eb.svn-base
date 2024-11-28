/**
 * \defgroup atUtilLib atUtilLib
 *  \brief The utility library in the Atlantis code.
 *
 * This file contains the headers for the public functions in the atlantis utility library.
 *
 *
 * 	Changes:
 *
 *
 * 	05-07-2012 Bec Gorton
 * 	Added the do_BrokerLinkage variable. This will be set to true if the user provides a linkage config file and linkage with the broker will occur.
 */

#include <atXMLUtil.h>
#include <atTracer.h>
#include <atlantisutil.h>
#include "atFisheryXMLUtil.h"
#include <atlantisMem.h>
#include <regex.h>

/*Global variables */
extern int verbose;
extern int fishtest;
extern int do_biology;
extern int do_physics;
extern int do_assess;
extern int do_manage;
extern int do_economics;
extern int do_BrokerLinkage;	/** This will be set to true if the user provides a linkage config file and linkage with the broker will occur. */
extern int do_CLAMLinkage;
extern int not_tracking_flux;

extern FILE *logfp;
extern FILE *haltfp;

extern int *groupEatOrder;
extern int *groupMovementOrder;
extern int *TLorder;
extern int *TLorder2;
extern int *TLorder3;
extern int *wcOrder;
extern int *sedOrder;
extern int *epiOrder;

/* Storage for new values at each time step */
extern double ***newwctr; /* new water column tracer */
extern double ***newsedtr; /* new sediment tracer */
extern double ***newicetr; /* new ice tracer */
extern double **newlandtr; /* new land tracer */

extern char **paramStrings;
extern char **cohortParamStrings;
extern char **spawnParamStrings;
extern char **RBCParamStrings;

/** Function declarations */
void Util_Close_Output_File(FILE *fp);

void Util_Setup_Species_Param_Strings(MSEBoxModel *bm);
char *Util_Get_Fishery_Name(MSEBoxModel *bm, int fisheryIndex);
int Util_Get_Fishery_Index(MSEBoxModel *bm, char *str);
int Util_Get_FG_Index_From_Token(MSEBoxModel *bm, char *speciesName, char *fileName, char *token, int quitOnFail);
int Util_Get_FG_Index(MSEBoxModel *bm, char *speciesName);
int Util_Get_Detritus_Index(MSEBoxModel *bm, int detIndex);
int Util_Get_Current_Stock_Index(MSEBoxModel *bm, int sp, int boxNumber, int relk);
void Util_Calculate_StockID(MSEBoxModel *bm);

/* Tracer lookup functions */
int Util_Get_Tracer_ID(MSEBoxModel *bm, char *tracerName);

/* Parameter change functions - these are used in functionality like the temperature and salinity change */
double Util_Get_Change_Scale(MSEBoxModel *bm, int num_changes, double **changeArray);
double Util_Get_Accumulative_Change_Scale(MSEBoxModel *bm, int num_changes, double **changeArray);
void Util_Copy_Change_Values(MSEBoxModel *bm, double **originalArray, double **newArray, int size, double start, double period, double mult);

void Util_Usage(void);
void Util_Free(MSEBoxModel *bm);
FILE *Util_fopen(MSEBoxModel *bm, const char *name, const char *mode);
int Util_ncopen(const char *destFolder, const char* name ,int mode);

/* General useage subroutine prototypes */
double Util_Lognorm_Distrib(double mu, double sigma, double x_b);
double Util_Logx_Result(double mu, double sigma);
double Util_Mich_Ment(double x, double m);

int Util_Read_Functional_Group_XML(MSEBoxModel *bm, char *fileName, FILE *llogfp);

/** Array functions **/
void Util_Init_5D_Double(double *****array, int dim1, int dim2, int dim3, int dim4, int dim5, double value);
void Util_Init_4D_Double(double ****array, int dim1, int dim2, int dim3, int dim4, double value);
void Util_Init_3D_Double(double ***array, int dim1, int dim2, int dim3, double value);
void Util_Init_2D_Double(double **array, int dim1, int dim2, double value);
void Util_Init_1D_Double(double *array, int dim1, double value);

void Util_Init_4D_Long_Double(long double ****array, int dim1, int dim2, int dim3, int dim4, long double value);
void Util_Init_3D_Long_Double(long double ***array, int dim1, int dim2, int dim3, long double value);
void Util_Init_2D_Long_Double(long double **array, int dim1, int dim2, long double value);

long double *Util_Alloc_Init_1D_Long_Double(int dim1, long double value);
long double **Util_Alloc_Init_2D_Long_Double(int dim1, int dim2, long double value);
long double ***Util_Alloc_Init_3D_Long_Double(int dim1, int dim2, int dim3, long double value);
long double ****Util_Alloc_Init_4D_Long_Double(int dim1, int dim2, int dim3, int dim4, long double value);

double *Util_Alloc_Init_1D_Double(int dim1, double value);
double **Util_Alloc_Init_2D_Double(int dim1, int dim2, double value);
double ***Util_Alloc_Init_3D_Double(int dim1, int dim2, int dim3, double value);
double ****Util_Alloc_Init_4D_Double(int dim1, int dim2, int dim3, int dim4, double value);
double *****Util_Alloc_Init_5D_Double(int dim1, int dim2, int dim3, int dim4, int dim5, double value);

void Util_Init_5D_Int(int ****array, int dim1, int dim2, int dim3, int dim4, int dim5, int value);
void Util_Init_4D_Int(int ****array, int dim1, int dim2, int dim3, int dim4, int value);
void Util_Init_3D_Int(int ***array, int dim1, int dim2, int dim3, int value);
void Util_Init_2D_Int(int **array, int dim1, int dim2, int value);
void Util_Init_1D_Int(int *array, int dim1, int value);

int *Util_Alloc_Init_1D_Int(int dim1, int value);
int **Util_Alloc_Init_2D_Int(int dim1, int dim2, int value);
int ***Util_Alloc_Init_3D_Int(int dim1, int dim2, int dim3, int value);
int ****Util_Alloc_Init_4D_Int(int dim1, int dim2, int dim3, int dim4, int value);
int *****Util_Alloc_Init_5D_Int(int dim1, int dim2, int dim3, int dim4, int dim5, int value);

/* Read in the fisheries definition input file */
int Util_Read_Fisheries_XML(MSEBoxModel *bm, char *fileName, FILE *llogfp);

/* Regular Expressions */
#define elementsof(x)	(sizeof(x)/sizeof(x[0]))
int at_compileRegExpression(regex_t *regBuffer, char *str);

void Util_XML_Read_Species_Param(MSEBoxModel *bm, char *fileName, xmlNodePtr parent, int paramID);
int Util_XML_Get_Param_Index(SpeciesParamStruct array[], int size, int paramID);

double Util_XML_Get_Species_Node_Value(xmlNodePtr parent, int paramID, int guild, SpeciesParamStruct *paramStruct);
extern SpeciesParamStruct speciesParamStructArray[tot_prms];
extern SpeciesParamStruct cohortSpeciesParamStructArray[cohortDepParams];
extern SpeciesParamStruct spawnSpeciesParamStructArray[spawnDepParams];
extern SpeciesParamStruct RBCSpeciesParamStructArray[num_rbc_species_params_id];

/* Reallocations */
double *Util_ReAlloc_1D_Double(int newdim, int olddim, double *oldarray, double value);
double ***Util_ReAlloc_3D_Double(int newdim1, int olddim1, int newdim2, int olddim2, int newdim3, int olddim3, double ***oldarray, double value);

int *Util_ReAlloc_1D_Int(int newdim, int olddim, int *oldarray, int value);
int **Util_ReAlloc_2D_Int(int newdim1, int olddim1, int newdim2, int olddim2, int **oldarray, int value);


/* Utility functions */
void Util_GenMnorm(double *vec, double *means, int *iseed, int np, double **tt, double *sg);
double Util_xnorm(double mean, double sigg, int *iiseed);
int Util_Check_NetCDF_Size(MSEBoxModel *bm, int fid, int *dump, char *fileName, int *index, int type);
