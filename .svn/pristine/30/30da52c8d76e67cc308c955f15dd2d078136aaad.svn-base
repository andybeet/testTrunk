/**
\file
\brief Header file for conversion code to take Atlantis variables to pass to EwE
\groupdef atLink

    File:		atlink.h
	Created:	5/9/2008
	Author:		Beth Fulton,   CSIRO Marine Laboratories,  2008
    Revisions:

*************/

typedef enum
{
	type_id,
	value_id
}LINK_DATA_ENTRIES;

typedef enum
{
	WEIGHTED_AVERAGE,
	FIXED_POINT
}LINK_CALC_TYPES;

/* Declarations so enums can be used in code */
LINK_CALC_TYPES link_calc_types;
DATA_TYPES link_data_types;
LINK_DATA_ENTRIES link_data_entries;

/* Array definitions */
char **ExportName;

/* Routine prototypes */
int GetArrays(MSEBoxModel *bm);
int SendArrays(MSEBoxModel *bm);

void CreateConvertedArray(MSEBoxModel *bm);
void LinkToModel(MSEBoxModel *bm);
void LoadCoversionParameters(MSEBoxModel *bm);
void UnpackArrays(MSEBoxModel *bm);
void WeightedAvg(MSEBoxModel *bm);
