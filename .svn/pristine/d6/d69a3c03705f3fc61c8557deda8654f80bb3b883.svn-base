
#define MAXNUMOUTPUTSTATES 10   /* The max number of output states - we are assuming the max is 10 */
#define numClamIndicators 7     /* The number of output indicators */
#define maxStrategyOptions 7    /* The max number of stategies in any scenario */
#define MAXSCENARIOS 13

typedef enum{
PSS_CHANGE,
Q_CHANGE,
HABITAT_CHANGE,
LINEAR_MORT_CHANGE
} CLAMLevers;


/**
 * The type of manager - this determines the weight allocated to each indicator.
 */
typedef enum{
	RANDOM_MANAGER,
	REACTIVE_MANAGER,
	SET_MANAGER
} CLAMMANAGERTYPE;

typedef struct{
	char leverName[200];
	CLAMLevers leverType;
	double start;
	double period;
	double mult;
	int groupCode;
	int cohort;
	int location;
	char targetName[200]; /* This could be the group code or the tracer name in the case of a PSS change */

} CLAMLeverStructure;

/* Information about each strategy */
typedef struct {
	char strategyName[2000];
	int CLAMStrategyID;
	int isEnabled;
	double lockoutPeriod;
	double timeLockOutStarted;
	int numLevers;
	int isLocked;
	double strategyCost;	/* The cost of this strategy */
	CLAMLeverStructure strategyLevers[20]; /* Assuming the max number of levers for each strategy */
} CLAMStrategyOptions;
/**
 *	These are in effect the radio buttons in the scenario tab of the CLAM GUI.
 *	These will change depending on the CLAM implementation. If we are going to use another CLAM model
 *	then we would want these to be read in from an input file.
 */
typedef struct{
	char scenarioOptionName[200];
	int maxOptionValue;
	int optionSelected;
	int doRun;
	CLAMStrategyOptions strategies[maxStrategyOptions];
	int clamScenarioIndex;
} CLAMScenarioOptionStructure;

extern CLAMScenarioOptionStructure *CLAMScenarioOptions;


typedef enum{
    STANDARD_DEV,
    PERCENT,
} CLAMOutputType;
/**
 *	These are the output indicators from the CLAM.
 *	These will change depending on the CLAM implementation. If we are going to use another CLAM model
 *	then we would want these to be read in from an input file.
 */
typedef struct{
	char   outputOptionName[200];
	int    optionSelected;
	double   OutputStates[MAXNUMOUTPUTSTATES];
	CLAMOutputType  outputType;    /* Either standard dev or percent */
	int    numStates;
} CLAMOutputOptionStructure;

/* Each of the scenarios and strategies - these are used to check the values in the input file are valid */
typedef struct{
	char scenarioName[STRLEN];
	char strategyNames[maxStrategyOptions][STRLEN];

} CLAMScenarioNames;

extern CLAMOutputOptionStructure CLAMOutputOptions[];
extern int numScenarioOptions;

extern CLAMScenarioNames CLAMScenarios[];


/*
 *  The information about each indicator we are watching within Atlantis.
 *  There should be a corresponding value in the CLAMOutputOptions array for each of these indicators
 */
typedef struct{
    double minValue;
    double maxValue;
    double idealValue;
    int indicatorLocation;     /* This is a bitwise indicator. 1 is wc, 2 is SED and 4 is Epi */
    int CLAMOutputIndex;    /* This links the indicator to the values in the CLAMOutputOptions array */
    char CLAMOutputName[500];
    int numScenarioLinks;
    int scenarioLinks[MAXSCENARIOS];
    int runScenarios;
    double targetChange;
    double currentValue;
    double indicatorWeighting;
} atCLAMIndicator;

extern char outputFolderPath[200];
extern char outputFileName[200];
extern char inputFileName[200];
extern char iniFileName[200];
extern char startedFileName[200];
extern char errorFileName[200];
extern char finishedFileName[200];
extern char commandLinesStr[1000];
extern char pathToICMSBuilder[200];
extern char ICMSProject[200];
extern char ICMSDLL[200];
extern double CLAMStartDay;
extern double CLAMFreq;
extern CLAMMANAGERTYPE CLAMManagerType;

extern double processWaitTime; /* the number of ms we will wait for the CLAM process to finish before we terminate it */
extern double startedWaitTime; /* the number of ms we will wait for the CLAM process to create the started file */

extern double *indicatorWeightings;  /* Weight for each of the indicators against each management objective */
extern double **impact; /* The impact value assigned to each outcome of the CLAM strategy */

extern int numIndicators;
extern atCLAMIndicator *indicators;
extern int indicatorCalculationFlag; /* If 1 then we are using an integration through the water column, else we are using a single layer value*/
extern int numIndicatorScenarioLinks;
extern double  ****outcomeValueArray;    /* Store the outcomes for each scenario */
extern double *scenarioOutcome;
extern int maxNumOutputStates;

/* Calculate results on the returned probabilities from the CLAM */
extern double ***expectedValue;
extern double ***variance;
extern double ***distance;

/* The euclidian distance values */
extern double **euclidDistance;

/* The 'best' outcome */
extern int bestScenarioIndex;
extern int bestStrategyIndex;
extern int considerAllIndicators;

/* Function definitions */
void Allocate_CLAM_Arrays(MSEBoxModel *bm);
void Free_CLAM_Arrays(MSEBoxModel *bm);
int Get_CLAM_IndicatorID(MSEBoxModel *bm, char *name);
void Start_CLAM(MSEBoxModel *bm, char *inputFileName);
int Read_CLAM_Link_Params(MSEBoxModel *bm, char *fileName);
void Create_CLAM_Ini_File(MSEBoxModel *bm);
void Parse_CLAM_OutputFile(MSEBoxModel *bm, int scenarioIndex, int strategyIndex);

/* Implementation Functions */
void Implement_Strategy(MSEBoxModel *bm, int scenarioIndex, int strategyIndex, FILE *llogfp);
