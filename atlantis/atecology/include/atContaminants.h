void Allocate_Contaiminants(MSEBoxModel *bm);
void Free_Contaminants(MSEBoxModel *bm);
void Init_Contaminants(MSEBoxModel *bm);

int Species_Contaminant_Uptake(MSEBoxModel *bm, BoxLayerValues *boxLayerInfo, HABITAT_TYPES habitat, double dtsz, int cIndex);
int Calculate_Species_Contaminant_Decay(MSEBoxModel *bm, BoxLayerValues *boxLayerInfo, HABITAT_TYPES habitat, double dtsz, int cIndex, double cGroupLevel);
int Init_Contaminant_Transfer_Values(MSEBoxModel *bm);
int Reconcile_Global_Contaminant_Values(MSEBoxModel *bm,HABITAT_TYPES habitatType);
int Group_Transfer_Contaminant(MSEBoxModel *bm, BoxLayerValues *boxLayerInfo, HABITAT_TYPES globalHabitat, HABITAT_TYPES habitat, int toGuild, int toCohort, int fromGuild, int fromCohort, double amountEaten, double ***spSPinfo, double initialBiomass, double dtsz, int need_prop, int caseGTC);

int Degrade_Contaminants(MSEBoxModel *bm, BoxLayerValues *boxLayerInfo, HABITAT_TYPES habitat, double dtsz, int cIndex, double cLevel);
//int Calculate_Contaminants_Flux(MSEBoxModel *bm,  double *tracerArray, double *fluxArray, HABITAT_TYPES habitat);
int Calculate_Contaminants_Flux(MSEBoxModel *bm, BoxLayerValues *boxLayerInfo, HABITAT_TYPES habitatType);
int Calculate_Species_Contaminant_Effects(MSEBoxModel *bm, int box, int clayer, double dtsz, HABITAT_TYPES habitatType);

void Contaminant_Record_Death(MSEBoxModel *bm, int sp, int cohort, double amount);
void Contaminant_Write_Contact_Record(MSEBoxModel *bm);
void Contaminant_Close_Contact_Record(MSEBoxModel *bm);
void Change_Contaminant_Levels(MSEBoxModel *bm, BoxLayerValues *boxLayerInfo, HABITAT_TYPES habitat, double dtsz, int offset);

void Contaminant_Update_ContactMort_Record(MSEBoxModel *bm, int sp, int cohort);
void Calculate_Contaminant_Q10_Corrections(MSEBoxModel *bm, BoxLayerValues *boxLayerInfo, HABITAT_TYPES habitat);

int Loose_Contaminant(MSEBoxModel *bm, BoxLayerValues *boxLayerInfo, HABITAT_TYPES globalHabitat, HABITAT_TYPES habitat, int species, int cohort, double amountLost, double ***spSPinfo, double initialBiomass);
int Gain_Contaminants(MSEBoxModel *bm, BoxLayerValues *boxLayerInfo, HABITAT_TYPES globalHabitat, HABITAT_TYPES habitat, int species, int cohort, double amountLost, double ***spSPinfo, double initialBiomass);

double Avoid_Contaminants(MSEBoxModel *bm, int groupIndex, int cohort, int box, int layer);

void Get_ContamMoveEffects(MSEBoxModel *bm, int species, int cohort, int box, int layer);
void Move_Vert_Contaminated(MSEBoxModel *bm, int sp, int cohort, double ****currentden);

typedef enum{
	linear_contaminant_uptake_id = 1,
	sigmoidal_uptake_id,
	piecewise_linear_id

} CONTAMINANT_UPTAKE_OPTION;

typedef enum {
	ongoingC_id,
	finalC_id
} CONTAMINANT_MORT_OPTION;
