/**
 *	\file
 *	\brief The public interface for the Assessment model.
 *
 *	This file contains all of the public routines in the Assessment module. All other functions
 *	within this module are private and shouldn't be called by any other module.
 *
 **/

/* Error types  - might be used by other models */
#define uniform_err 0
#define normal_err 1
#define lognorm_err 2
#define negbiom_type 3   // currently only used for a cdf


/* Assessment output */
#define est_med_recruit_id 0
#define est_top_recruit_id 1
#define est_bot_recruit_id 2
#define est_med_prm2_id 3
#define est_top_prm2_id 4
#define est_bot_prm2_id 5
#define est_med_stock_id 6
#define est_top_stock_id 7
#define est_bot_stock_id 8
#define est_SS_id 9
#define est_FrefA_id 10
#define est_FrefC_id 11
#define est_Fcurr_id 12						/* Fishing mortality estimate */
#define est_M_id 13							/* Mortality estimates */
#define est_Ccurr_id 14						/* Current catch estimate */
#define est_Cslope_id 15
#define est_numyr_id 16
#define est_CPUEnow_id 17
#define est_CPUEstart_id 18
#define est_balpha_id 19
#define est_bbeta_id 20
#define est_trend_id 21

#define K_num_assess_prms 22

/* Assessment selectivity curve types */
#define const_q 0
#define logistic_q 1
#define normal_q 2
#define lognorm_q 3
#define gamma_q 4

/* For phys[id][box][record_in_time]
 property                 id
 water column depth        0
 sediment depth            1
 salinity                  2
 light in water column     3
 light at sediment surface 4
 oxygen in sediment        5
 don in water column       6
 din in water column       7
 volume                    8
 */

#define wc_depth_id 0
#define sed_depth_id 1
#define salinity_id 2
#define light_wc_id 3
#define light_sed_id 4
#define o2_sed_id 5
#define o2_wc_id 6
#define TSS_id 7
#define don_wc_id 8
#define din_wc_id 9
#define vol_id 10
#define temperature_id 11
#define num_sampled_phy_id 12

//#define temperature_id 11
//#define chl_a_id 12
//#define num_sampled_phy_id 13

/** Initialise the assessment module */
void Assess_Init(MSEBoxModel *bm, FILE *llogfp);
void Assess_Free(MSEBoxModel *bm);

int Tier_Assessment_Free(MSEBoxModel *bm);

/* Sampling error routines */
//double Assess_errfun (double xlim);
double Assess_Get_Physical_Sampled_Value(MSEBoxModel *bm, int sample, int zone);
double Assess_Get_Physical_Sampled_SD(MSEBoxModel *bm, int sample, int zone);
double Assess_Get_Process_Sampled_Value(MSEBoxModel *bm, int sample, int b);
double Assess_Get_Biomass_Sampled_Value(MSEBoxModel *bm, int sample, int b);

/* Assessment model routines */
void Assess_Resources(MSEBoxModel *bm, FILE *llogfp);
void Assess_Pseudo_Estimate_Prm(MSEBoxModel *bm, int sp, FILE *ofp);

/** Annual routines */
void Assess_Annual_Schedule(MSEBoxModel *bm, FILE *llogfp);


void Assess_Set_Name_ID(MSEBoxModel *bm);


