
typedef enum
{
	TYPE_BIOMASS,
	TYPE_PRODUCTION,
	TYPE_MORTALITY,
	TYPE_NUTRIENTS,
	TYPE_PHYSCHEM
}DATA_TYPES;

typedef struct 
{
	/* Link parameters */
    int refresh;
	double numdt;
	int ingroup;
	int outgroup;
	int io_nrow;
	int io_ncol;
	double cumday;
	double EwEdt;
	int timehandling;
	int spacehandling;
	double ulx;
	double uly;
	double ullat;
	double ullong;
	double ulref;
	double io_area;	
}convertPRMSStruct;

typedef struct {
	char 		*parameterName;
	double 		**values;
	double 		**area;
	DATA_TYPES 	dataType;	

}LinkTransportStruct;

typedef struct {
	
	/* The link parameters structure */
    convertPRMSStruct convertPRMS;
    LinkTransportStruct *ImportArray;
	LinkTransportStruct *ExportArray;		
    
}EwELinkStruct;

EwELinkStruct *eweLinkModule;
