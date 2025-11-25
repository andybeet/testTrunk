/*******************************************************************//**
 \ingroup atUtil

 File:           atlantisutil.c

 Created:        21/9/2004

 Author:         Beth Fulton
 CSIRO Division of Marine Resaech

 Purpose:        Utility routines for atlantis mse box model


 Revisions:      2/1/2008 Moved all code out of atlantismain except the main function.


 15-Feb_2008 Bec Gorton

 Fixed a bug in the initialisation of the bm->regionalData array.
 It was indexing the array as

 bm->RegionalData[i][b][sp] = 1.0;

 and it should be

 bm->RegionalData[sp][b][i] = 1.0;

 26-05-2008 Bec Gorton
 Changed all references from bm->VERTind to the tracer arrays associated with
 each functional group.

 14-06-2008 Bec Gorton
 Updated code to dynamically build the bm->dinfo array instead of
 reading in a blank file on start up.

 23-12-2008 Resized the following arrays:

 bm->mS = Util_Alloc_Init_3D_Double(2, 4, bm->K_num_vert_sp, 0.0); to

 bm->mS = Util_Alloc_Init_2D_Double(4, bm->K_num_vert_sp, 0.0); As the implicit
 mortality due to seabirds and fishing is now being combined into a single parameter.

 The age diet availability of invertebrates for consumption by vertebrates
 bm->pSPageeat = Util_Alloc_Init_3D_Double(Not_age_specific_id, bm->K_num_max_cohort, bm->K_num_vert_sp, 0.0); to:
 bm->pSPageeat = Util_Alloc_Init_3D_Double(bm->K_num_tot_sp, bm->K_num_max_cohort, bm->K_num_tot_sp, 0.0);
 This means this array can be referenced by the groupCode instead of the AGE_PREY_id values.

 02-02-2009 Bec Gorton
 Moved the freeing of the economic arrays in the bm structure into the Economic_Free
 function in this file instead of them being freed in the shutdownmodel function. This
 means that we can check to see if the economic module is on - if not these
 arrays should not be freed.
 Also added checks before calling freeassess and freeeconomic so these are
 only called if these modules are on in the model.

 02-04-2009 Bec Gorton
 Moved the function call that frees up the box geometry to after the call
 that frees up the eoclogy structure. The ecology structure now needs to
 free up the home range memory allocated in setupHomeRanges. A new function
 has been added to atbiolsetup.c that will free this up. This needs to be
 called before the code to free up the boxes.


 22-04-2009 Bec Gorton
 Removed the bm->fisheryParamNAME and SP_FISHERYprmsName arrays.
 Moved all the IO related functions to atIOUtil.c

 10-06-2009 Bec Gorton
 Added support for rolling log files. The log file will be closed when it gets larger
 than MAX_LOG_FILE_SIZE and a new one will be opened. The log file names are simply logD where D
 starts at 1 and is incremented when a new log file is opened.


 22-May-2009 Bec Gorton
 Changed the ncopen calls to use sjw_ncopen. This function checks that the
 netcdf file exits before calling ncopen which will crash if the file is not
 found.


 07-07-2009 Bec Gorton
 Moved the speciesParamStrings into this file - these are now setup in the atlantisUtil
 lib so they can be used in the fishery and management libs.
 Also started adding support for the fishery definition file to be read in. So there is now
 an additional command line argument - the fishery definition file. Useage now includes:
 -p fisheries.csv

 28-10-2009 Bec Gorton
 Changed array sizes and references to get rid of the K_num_fished_sp, K_num_impacted_sp
 and K_num_fishedTac references.

 30-10-2009 Bec Gorton
 Added the Util_Close_Output_File function.

 03-11-2009 Bec Gorton
 Removed the SPtoCATid array as these values are no longer required

 04-11-2009 Bec Gorton
 Merged in Beths bycatch incentive code - revision 961.

 05-11-2009 Bec Gorton
 Trunk merge 1064. Added support for the Q10 parameter to be specified for each functional group.

 14-12-2009 Bec Gorton
 Trunk merge 1424. Added the calcMLinearMort, calcMQuadMort and calcMPredMort arrays to store different mortality values.


 15-12-2009 Bec Gorton
 Trunk merge 1439 - Added the code to read in the clamlinkage input file.

 19-01-2010 Bec Gorton
 Removed the bm->PREYid array as its no longer used!
 Updated the headers that are included.
 Removed the code that read in the fstatistic_blank.nc netcdf file. This is replaced
 with a call to the harvest library to build the tracers as the biology tracers are now created.
 Moved the Distance_to_Port call into the Manage_init function so it can be private to the
 management library.
 Deleted the bm->SPnid2guild, bm->SPmove2guild and bm->PREYid arrays.
 The FFCDR array has been moved out of the bm and into the harvest library,
 ( as have the bm->tscatchid, bm->tsdiscardid, bm->harvestindx, bm->harvestindxNAME and bm->TotCumCatch variables).

 29-01-2010 Bec Gorton
 Removed the bm->fisheryName array. How using the FisheryArray[nf].fisheryCode.

 16-02-2010 Bec Gorton
 Added the Util_Setup_Species_Param_Strings function to this library - it was in the biology lib.

 25-02-2010 Bec Gorton
 Moved the model setup code back into the atlantis main module. We no longer need to use the test interface that
 meant they had to separate.

 26-02-2010 Bec Gorton
 Moved the functions that were defined in the atlog.h header file:
 Util_Logx_Result
 Util_Lognorm_Distrib
 Util_Mich_Ment
 out of the eoclogy lib into the util library and deleted the atlog.h file from the repository/.


 04-05-2010 Bec Gorton
 Added new functions to get calculate the scaling factor from the given change parameter arrays.

 05-05-2010 Bec Gorton
 Added code to Util_Setup_Species_Param_Strings() that will set up the paramString values
 for the new Q10 algorithm.

 There is now an additional method of implementing temperature correction based on work by Gary G.
 The new parameters to implement this are now read in from the input file. These are not compulsory unless the
 q10_method_id method is set to 1 in which case ecol_readSpeciesParam will display an error message and quit.

 *********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <signal.h>
#include <math.h>
#include <ctype.h>
#include <time.h>
#include <string.h>
#include <sjwlib.h>
#include <netcdf.h>
#include <sys/stat.h>
#include <atlantisboxmodel.h>
#include <atUtilLib.h>
#include <atEcologyLib.h>
#include <atHarvestLib.h>
#include <atManageLib.h>
#include <atImplementationLib.h>


char **paramStrings;
char **cohortParamStrings;
char **spawnParamStrings;
char **RBCParamStrings;

char *cohortStrings[] =
	{ "juv", "adult" };

char *sexStrings[] =
    { "F", "M" };

void Util_Free(MSEBoxModel *bm){
	Free_Fishery_Def_Memory(bm);
	Free_Functional_Group_Memory(bm);
}
/**
 *	\brief These functions are to get around the issue that the economic input files
 *	use 'recfish' for the recfish_id parameter but the manage input files use 'REC'.
 *
 */

char *Util_Get_Fishery_Name(MSEBoxModel *bm, int fisheryIndex) {
//	if (fisheryIndex == recfish_id)
//		return "REC";
	return FisheryArray[fisheryIndex].fisheryCode;
}

/**
 * \brief Close the given output file if the given pointer is not null.
 */
void Util_Close_Output_File(FILE *fp) {
	if (fp != NULL){
		fflush(fp);
		fclose(fp);
	}
}

/**
 * \brief Return the index of the given fishery code.
 *
 * This function will quit if the given string is not recognised.
 *
 *
 */
int Util_Get_Fishery_Index(MSEBoxModel *bm, char *strPtr) {
	int returnValue = -1;
	int nf;
//
//	if (strcmp(strPtr, "REC") == 0)
//		return recfish_id;
//	else {
		for (nf = 0; nf < bm->K_num_fisheries; nf++) {
			if (strcmp(strPtr, FisheryArray[nf].fisheryCode) == 0)
				return nf;
		}
	//}

	if (returnValue == -1)
		quit("Util_Get_Fishery_Name - fishery code '%s' not recognised", strPtr);

	return returnValue;
}

/**
 *	\brief Setup the array of species parameter strings.
 *
 *	These are used to read in values from the input file.
 *
 */
void Util_Setup_Species_Param_Strings(MSEBoxModel *bm) {
	// Need to move this into a seperate function.
	paramStrings = (char **) c_alloc2d(25, tot_prms);
	cohortParamStrings = (char **) c_alloc2d(25, cohortDepParams);
	spawnParamStrings = (char **) c_alloc2d(25, spawnDepParams);

	sprintf(paramStrings[flag_id], "%s", "flag");
	sprintf(paramStrings[flagdem_id], "%s", "flagdem");
	sprintf(paramStrings[flagplankfish_id], "%s", "flagplankfish");
	sprintf(paramStrings[flagrecruit_id], "%s", "flagrecruit");
	sprintf(paramStrings[flagrecpeak_id], "%s", "flagrecpeak");
	//sprintf(paramStrings[flaglocalrecruit_id], "%s", "flaglocalrecruit");
	sprintf(paramStrings[flagbearlive_id], "%s", "flagbearlive");
	sprintf(paramStrings[flagmother_id], "%s", "flagmother");
	sprintf(paramStrings[feed_while_spawn_id], "%s", "feed_while_spawn");
	sprintf(paramStrings[flagstocking_id], "%s", "flagstocking");
    sprintf(paramStrings[flagkeep_plusgroup_id], "%s", "flagkeep_plusgroup");
	
    sprintf(paramStrings[flagq10receff_id], "%s", "flagq10receff");
	sprintf(paramStrings[flagq10eff_id], "%s", "flagq10eff");
	sprintf(paramStrings[q10_id], "%s", "q10_");
	sprintf(paramStrings[q10_method_id], "%s", "q10_method_");
	sprintf(paramStrings[q10_optimal_temp_id], "%s", "q10_optimal_temp_");
	sprintf(paramStrings[q10_correction_id], "%s", "q10_correction_");
	sprintf(paramStrings[q10_method_id], "%s", "q10_method_");
	sprintf(paramStrings[optimal_pH_id], "%s", "optimal_pH_");
	sprintf(paramStrings[salt_correction_id], "%s", "salt_correction_");
	sprintf(paramStrings[pH_correction_id], "%s", "pH_correction_");
	sprintf(paramStrings[flagtempsensitive_id], "%s", "flagtempsensitive");
	sprintf(paramStrings[temp_coefftA_id], "%s", "temp_coefftA_");
	sprintf(paramStrings[flagpHsensitive_id], "%s", "flagpHsensitive_");
	sprintf(paramStrings[flagfecundsensitive_id], "%s", "flagfecundsensitive_");
	sprintf(paramStrings[flagSaltSensitive_id], "%s", "flagSaltSensitive_");
	sprintf(paramStrings[flagnutvaleffect_id], "%s", "flagnutvaleffect_");
	sprintf(paramStrings[flagpredavaileffect_id], "%s", "flagpredavaileffect_");
	sprintf(paramStrings[flagcontract_tol_id], "%s", "flagcontract_tol_");
	sprintf(paramStrings[contract_tol_id], "%s", "contract_tol_");
	sprintf(paramStrings[KN_pH_id], "%s", "KN_pH_");
	sprintf(paramStrings[pH_constA_id], "%s", "pH_constA_");
	sprintf(paramStrings[pH_constB_id], "%s", "pH_constB_");
    sprintf(paramStrings[pH_constC_id], "%s", "pH_constC_");
	sprintf(paramStrings[min_pH_id], "%s", "min_pH_");
	sprintf(paramStrings[max_pH_id], "%s", "max_pH_");
	sprintf(paramStrings[pHsensitive_model_id], "%s", "pHsensitive_model_");
	sprintf(paramStrings[flagchannel_id], "%s", "flagchannel");
	sprintf(paramStrings[pHmortstart_id], "%s", "pHmortstart_");
	sprintf(paramStrings[pHmortA_id], "%s", "pHmortA_");
	sprintf(paramStrings[pHmortB_id], "%s", "pHmortB_");
	sprintf(paramStrings[pHmortmid_id], "%s", "pHmortmid_");
    
    sprintf(paramStrings[light_coefft_id], "%s", "light_coefft_");
    sprintf(paramStrings[noise_coefft_id], "%s", "noise_coefft_");

	sprintf(paramStrings[catcheater_id], "%s", "catcheater");
    sprintf(paramStrings[flagactive_id], "%s", "day");
	sprintf(paramStrings[vla_T15_id], "%s", "vla_");

	sprintf(paramStrings[Kcov_juv_id], "%s", "Kcov_juv_");
	sprintf(paramStrings[Bcov_juv_id], "%s", "Bcov_juv_");
	sprintf(paramStrings[Acov_juv_id], "%s", "Acov_juv_");
	sprintf(paramStrings[Kcov_ad_id], "%s", "Kcov_ad_");
	sprintf(paramStrings[Bcov_ad_id], "%s", "Bcov_ad_");
	sprintf(paramStrings[Acov_ad_id], "%s", "Acov_ad_");

	sprintf(paramStrings[KL_id], "%s", "KL_");
	sprintf(paramStrings[KU_id], "%s", "KU_");
	sprintf(paramStrings[KUP_id], "%s", "KUP_");
	sprintf(paramStrings[KLP_id], "%s", "KLP_");
	sprintf(paramStrings[Kmax_coefft_id], "%s", "Kmax_coefft_");
	sprintf(paramStrings[extra_feed_id], "%s", "_extra_feed");
	
	sprintf(paramStrings[KDEP_id], "%s", "KDEP_");
	sprintf(paramStrings[KWSR_id], "%s", "KWSR_");
	sprintf(paramStrings[KWRR_id], "%s", "KWRR_");
	sprintf(paramStrings[recover_start_id], "%s", "recover_start_");
	sprintf(paramStrings[recover_mult_id], "%s", "recover_mult_");
	sprintf(paramStrings[BHbeta_id], "%s", "BHbeta_");
	sprintf(paramStrings[BHalpha_id], "%s", "BHalpha_");
	sprintf(paramStrings[Rbeta_id], "%s", "Rbeta_");
	sprintf(paramStrings[Ralpha_id], "%s", "Ralpha_");
	sprintf(paramStrings[PP_id], "%s", "PP_");
	sprintf(paramStrings[hta_id], "%s", "hta_");
	sprintf(paramStrings[htb_id], "%s", "htb_");
	sprintf(paramStrings[pR_id], "%s", "pR_");
	sprintf(paramStrings[prop_spawn_lost_id], "%s", "prop_spawn_lost_");
	sprintf(paramStrings[jack_a_id], "%s", "jack_a_");
	sprintf(paramStrings[jack_b_id], "%s", "jack_b_");
	sprintf(paramStrings[intersp_depend_recruit_id], "%s", "intersp_depend_recruit_");
	sprintf(paramStrings[intersp_depend_sp_id], "%s", "intersp_depend_sp_");
	sprintf(paramStrings[intersp_depend_scale_id], "%s", "intersp_depend_scale_");
	sprintf(paramStrings[aquacult_fry_id], "%s", "_aquacult_fry");

	sprintf(paramStrings[li_a_id], "%s", "li_a_");
	sprintf(paramStrings[li_b_id], "%s", "li_b_");
    sprintf(paramStrings[linf_id], "%s", "linf_");
    sprintf(paramStrings[Kbert_id], "%s", "Kbert_");
    sprintf(paramStrings[tzero_id], "%s", "tzero_");
    
    sprintf(paramStrings[min_li_mat_id], "%s", "min_li_mat_");
	sprintf(paramStrings[KA_id], "%s", "KA_");
	sprintf(paramStrings[KB_id], "%s", "KB_");

    sprintf(paramStrings[RSmax_id], "%s", "RSmax_");
    sprintf(paramStrings[RSmid_id], "%s", "RSmid_");
    sprintf(paramStrings[RSslope_id], "%s", "RSslope_");
    sprintf(paramStrings[RSprop_id], "%s", "RSprop_");
    sprintf(paramStrings[SNcost_id], "%s", "SNcost_");
    sprintf(paramStrings[RNcost_id], "%s", "RNcost_");
    sprintf(paramStrings[RSstarve_id], "%s", "RNstarve_");
    
	sprintf(paramStrings[overwinterStartTofY_id], "%s", "overwinterStartTofY_");
	sprintf(paramStrings[overwinterEndTofY_id], "%s", "overwinterEndTofY_");
	sprintf(paramStrings[overwinterStartTemp_id], "%s", "overwinterStartTemp_");
	sprintf(paramStrings[overwinterEndTemp_id], "%s", "overwinterEndTemp_");
	sprintf(paramStrings[crit_mum_id], "%s", "crit_mum_");
	sprintf(paramStrings[crit_nut_id], "%s", "crit_nut_");
	sprintf(paramStrings[crit_temp_id], "%s", "crit_temp_");
	sprintf(paramStrings[encyst_rate_id], "%s", "encyst_rate_");
	sprintf(paramStrings[hatch_rate_id], "%s", "hatch_rate_");
	sprintf(paramStrings[encyst_period_id], "%s", "encyst_period_");
	sprintf(paramStrings[flagencyst_id], "%s", "flagencyst_");

	sprintf(paramStrings[max_prop_shift_id], "%s", "max_prop_shift_");
	sprintf(paramStrings[inheritance_id], "%s", "inheritance_");
	sprintf(paramStrings[trait_variance_id], "%s", "trait_variance_");
	sprintf(paramStrings[min_trait_variance_id], "%s", "min_trait_variance_");

	sprintf(paramStrings[E1_id], "%s", "E_");
	sprintf(paramStrings[E2_id], "%s", "EPlant_");
	sprintf(paramStrings[E3_id], "%s", "EDL_");
	sprintf(paramStrings[E4_id], "%s", "EDR_");

	sprintf(paramStrings[KSPA_id], "%s", "KSPA_");
	sprintf(paramStrings[FSP_id], "%s", "FSP_");
	sprintf(paramStrings[Recruit_Time_id], "%s", "_Recruit_Time");
    
	sprintf(paramStrings[Recruit_Period_id], "%s", "Recruit_Period_");
    sprintf(paramStrings[cohort_recruit_entry_id], "%s", "_cohort_recruit_entry");    
	sprintf(paramStrings[spawn_period_id], "%s", "_spawn_period");
	sprintf(paramStrings[age_mat_id], "%s", "_age_mat");
	sprintf(paramStrings[FDMort_id], "%s", "FDM_");
	sprintf(paramStrings[FDG_id], "%s", "FDG_");

	sprintf(paramStrings[log_mult_id], "%s", "_log_mult");
    sprintf(paramStrings[norm_sigma_id], "%s", "_norm_sigma");
    sprintf(paramStrings[flag_recruit_stochastic_id], "%s", "_flag_recruit_stochastic");
	sprintf(paramStrings[min_spawn_temp_id], "%s", "_min_spawn_temp");
	sprintf(paramStrings[max_spawn_temp_id], "%s", "_max_spawn_temp");
	sprintf(paramStrings[min_spawn_salt_id], "%s", "_min_spawn_salt");
	sprintf(paramStrings[max_spawn_salt_id], "%s", "_max_spawn_salt");
	sprintf(paramStrings[min_move_temp_id], "%s", "_min_move_temp");
	sprintf(paramStrings[max_move_temp_id], "%s", "_max_move_temp");
	sprintf(paramStrings[min_move_salt_id], "%s", "_min_move_salt");
	sprintf(paramStrings[max_move_salt_id], "%s", "_max_move_salt");
	sprintf(paramStrings[min_O2_id], "%s", "_min_O2");
    
    sprintf(paramStrings[K_temp_const_id], "%s", "_K_temp_const");
    sprintf(paramStrings[K_salt_const_id], "%s", "_K_salt_const");
    sprintf(paramStrings[K_o2_const_id], "%s", "_K_o2_const");
    
	sprintf(paramStrings[predcase_id], "%s", "predcase_");

	sprintf(paramStrings[ddepend_move_id], "%s", "_ddepend_move");

	sprintf(paramStrings[vlb_id], "%s", "vlb_");

    sprintf(paramStrings[turbid_refuge_id], "%s", "turbid_refuge_");
    
    sprintf(paramStrings[FDGDL_id], "%s", "FDGDL_");
	sprintf(paramStrings[FDGDR_id], "%s", "FDGDR_");
	sprintf(paramStrings[age_structured_prey_id], "%s", "age_structured_prey_");
	sprintf(paramStrings[p_split_id], "%s", "p_split_");
	sprintf(paramStrings[KTUR_id], "%s", "KTUR_");
	sprintf(paramStrings[KIRR_id], "%s", "KIRR_");
	sprintf(paramStrings[vl_id], "%s", "vl_");
	sprintf(paramStrings[ht_id], "%s", "ht_");
    sprintf(paramStrings[hvm_id], "%s", "hvm_");

	sprintf(paramStrings[KN_id], "%s", "KN_");
	sprintf(paramStrings[KS_id], "%s", "KS_");
	sprintf(paramStrings[KF_id], "%s", "KF_");
	sprintf(paramStrings[KP_id], "%s", "KP_");
	sprintf(paramStrings[KI_T15_id], "%s", "KI_*_T15");
	sprintf(paramStrings[Beta_D_id], "%s", "Beta_D_");
	sprintf(paramStrings[PBmax_D_id], "%s", "PBmax_D_");
	sprintf(paramStrings[ICE_KDEP_id], "%s", "ICE_KDEP_");

	sprintf(paramStrings[KI_L_T15_id], "%s", "L_KI_*_T15");
	sprintf(paramStrings[Kext_id], "%s", "Kext_");
	sprintf(paramStrings[Ksub_id], "%s", "Ksub_");
	sprintf(paramStrings[KNepi_id], "%s", "KN_epi_");
	sprintf(paramStrings[KsubEpi_id], "%s", "KsubEpi_");
	sprintf(paramStrings[Ktrans_id], "%s", "Ktrans_");

		/* Coral and rugosity related parameters */
	sprintf(paramStrings[bleach_periodA_id], "%s", "_bleach_periodA");
	sprintf(paramStrings[bleach_periodB_id], "%s", "_bleach_periodB");
	sprintf(paramStrings[mBleach_id], "%s", "_mBleach");
	sprintf(paramStrings[bleaching_rate_id], "%s", "_bleaching_rate");
	sprintf(paramStrings[bleach_recovery_rate_id], "%s", "_bleach_recovery_rate");
	sprintf(paramStrings[bleach_tempshift_id], "%s", "_bleach_tempshift");
	sprintf(paramStrings[bleach_growshift_id], "%s", "_bleach_growshift");
	sprintf(paramStrings[bleach_temp_id], "%s", "_bleach_temp");
    sprintf(paramStrings[min_bleach_temp_id], "%s", "_min_bleach_temp");
    sprintf(paramStrings[prop_zooxanth_id], "%s", "_prop_zooxanth");
	sprintf(paramStrings[DHW_thresh_id], "%s", "_DHW_thresh");
    sprintf(paramStrings[threshdepth_id], "%s", "_threshdepth_id");
    sprintf(paramStrings[depmum_scalar_id], "%s", "_depmum_scalar_id");

    sprintf(paramStrings[min_bleach_salt_id], "%s", "_min_bleach_salt");
    sprintf(paramStrings[max_bleach_salt_id], "%s", "_max_bleach_salt");
    
    sprintf(paramStrings[rugFeedScalar_id], "%s", "_rugFeedScalar");
	sprintf(paramStrings[HostRemin_id], "%s", "_HostRemin");
	sprintf(paramStrings[calcifRefBaseline_id], "%s", "_calcifRefBaseline");
	sprintf(paramStrings[calcifTconst_id], "%s", "_calcifTconst");
	sprintf(paramStrings[calcifTcoefft_id], "%s", "_calcifTcoefft");
	sprintf(paramStrings[calcifTopt_id], "%s", "_calcifTopt");
	sprintf(paramStrings[calcifLambda_id], "%s", "_calcifLambda");
	sprintf(paramStrings[FeedLightThresh_id], "%s", "_FeedLightThresh");
	sprintf(paramStrings[PropLightFeed_id], "%s", "_PropLightFeed");
	sprintf(paramStrings[rug_erode_id], "%s", "_rug_erode");
	sprintf(paramStrings[rug_bleacherode_id], "%s", "_rug_bleacherode");
	sprintf(paramStrings[rug_factor_id], "%s", "_rug_factor");
	sprintf(paramStrings[colony_ha_id], "%s", "_colony_ha");
	sprintf(paramStrings[coral_overgrow_id], "%s", "_coral_overgrow");
	sprintf(paramStrings[coral_compete_id], "%s", "_coral_compete");
	sprintf(paramStrings[coral_max_accel_trans_id], "%s", "_coral_max_accel_trans");
	sprintf(paramStrings[coral_max_accelA_id], "%s", "_coral_max_accelA");
	sprintf(paramStrings[coral_max_accelB_id], "%s", "_coral_max_accelB");
	sprintf(paramStrings[CrecruitA_id], "%s", "_CrecruitA");
	sprintf(paramStrings[CrecruitB_id], "%s", "_CrecruitB");
	sprintf(paramStrings[CrecruitC_id], "%s", "_CrecruitC");
	sprintf(paramStrings[rec_HabDepend_id], "%s", "rec_HabDepend");
	sprintf(paramStrings[RugCover_scalar_id], "%s", "_RugCover_scalar");

    sprintf(paramStrings[sponge_overgrow_id], "%s", "_sponge_overgrow");
    sprintf(paramStrings[sponge_compete_id], "%s", "_sponge_compete");
    sprintf(paramStrings[Ksmother_A_id], "%s", "_Ksmother_A");
    sprintf(paramStrings[Ksmother_B_id], "%s", "_Ksmother_B");
    sprintf(paramStrings[Vmax_deltaSi_id], "%s", "_Vmax_deltaSi");
    sprintf(paramStrings[Km_deltaSi_id], "%s", "_Km_deltaSi");
    sprintf(paramStrings[rug_erode_sponge_id], "%s", "_rug_erode_sponge");
    
    sprintf(cohortParamStrings[rugosity_inc_id], "%s", "_rugosity_inc");
	sprintf(cohortParamStrings[rugosity_dec_id], "%s", "_rugosity_dec");
	sprintf(cohortParamStrings[colony_diam_id], "%s", "_colony_diam");

	/* Primary producer P variables */
	sprintf(paramStrings[P_max_uptake_id], "%s", "P_uptake_");
	sprintf(paramStrings[P_uptake_scale_id], "%s", "P_scale_uptake_");
	sprintf(paramStrings[P_concp_id], "%s", "P_concp_");
	sprintf(paramStrings[C_max_uptake_id], "%s", "C_uptake_");
	sprintf(paramStrings[C_uptake_scale_id], "%s", "C_scale_uptake_");
	sprintf(paramStrings[C_concp_id], "%s", "C_concp_");

	sprintf(paramStrings[P_min_internal_id], "%s", "P_min_internal_");
	sprintf(paramStrings[P_max_internal_id], "%s", "P_max_internal_");

	sprintf(paramStrings[Psa_min_id], "%s", "PSA_min_");
	sprintf(paramStrings[C_min_id], "%s", "C_min_");
	sprintf(paramStrings[phyto_resp_rate_id], "%s", "Phyto_Resp_Rate_");

	sprintf(paramStrings[KLYS_id], "%s", "KLYS_");
	sprintf(paramStrings[mD_id], "%s", "mD_");
	sprintf(paramStrings[mStarve_id], "%s", "mStarve_");
	sprintf(paramStrings[mT_id], "%s", "mT_");
	sprintf(paramStrings[mO_id], "%s", "mO_");
	sprintf(paramStrings[mS_id], "%s", "mS_");
	sprintf(paramStrings[mS_T15_id], "%s", "mS_*_T15");
	sprintf(paramStrings[KO2_id], "%s", "KO2_");
	sprintf(paramStrings[KO2LIM_id], "%s", "KO2LIM_");
	sprintf(paramStrings[sp_remin_contrib_id], "%s", "_remin_contrib");

	sprintf(paramStrings[max_id], "%s", "max");
	sprintf(paramStrings[low_id], "%s", "_low");
	sprintf(paramStrings[thresh_id], "%s", "thresh");
	sprintf(paramStrings[sat_id], "%s", "_sat");
	sprintf(paramStrings[FSBDR_id], "%s", "FSBDR_");

	sprintf(paramStrings[flux_thresh_id], "%s", "thresh");
	sprintf(paramStrings[flux_damp_id], "%s", "damp");
	sprintf(paramStrings[flag_lim_id], "%s", "flaglim");

	sprintf(cohortParamStrings[mL_T15_id], "%s", "_mL");
	sprintf(cohortParamStrings[mQ_T15_id], "%s", "_mQ");
	sprintf(spawnParamStrings[Time_Spawn_id], "%s", "_Time_Spawn");
	sprintf(spawnParamStrings[Time_Age_id], "%s", "_Time_Age");

    sprintf(cohortParamStrings[L_turbid_id], "%s", "_L_turbid");
    sprintf(cohortParamStrings[a_turbid_id], "%s", "_a_turbid");
    sprintf(cohortParamStrings[b_turbid_id], "%s", "_b_turbid");
    
    sprintf(paramStrings[mindepth_id], "%s", "_mindepth");
	sprintf(paramStrings[maxdepth_id], "%s", "_maxdepth");
    sprintf(paramStrings[maxtotdepth_id], "%s", "_maxtotdepth");
	sprintf(paramStrings[homerangerad_id], "%s", "_homerangerad");
	sprintf(paramStrings[rangeoverlap_id], "%s", "_overlap");
	sprintf(paramStrings[Speed_id], "%s", "Speed_");

	/* The fisheries params */
	sprintf(paramStrings[flagfish_id], "%s", "flagfish");
	sprintf(paramStrings[access_thru_wc_id], "%s", "flag_access_thru_wc_");
	sprintf(paramStrings[tier_id], "%s", "tier");
	sprintf(paramStrings[regionalSP_id], "%s", "regionalSP");
	sprintf(paramStrings[basketSP_id], "%s", "basketSP");
	sprintf(paramStrings[basket_size_id], "%s", "basket_size");
	sprintf(paramStrings[coType_id], "%s", "coType_");
	sprintf(paramStrings[tac_resetperiod_id], "%s", "tac_resetperiod");
	sprintf(paramStrings[aquacult_age_harvest_id], "%s", "_age_harvest");
    sprintf(paramStrings[cpue_cdf_poor_r_id], "%s", "cpue_cdf_poor_r_");
    sprintf(paramStrings[cpue_cdf_poor_p_id], "%s", "cpue_cdf_poor_p_");
    sprintf(paramStrings[cpue_cdf_top_r_id], "%s", "cpue_cdf_top_r_");
    sprintf(paramStrings[cpue_cdf_top_p_id], "%s", "cpue_cdf_top_p_");
    
    sprintf(paramStrings[whichRAssess_id], "%s", "whichRAssess_");

	bm->SP_FISHERYprmsName = (char **) c_alloc2d(50, tot_sp_specif_fishing_prms);
	sprintf(bm->SP_FISHERYprmsName[sel_id], "%s", "sel_");
	sprintf(bm->SP_FISHERYprmsName[flagQchange_id], "%s", "flagQchange_");
	sprintf(bm->SP_FISHERYprmsName[Q_num_changes_id], "%s", "_q_changes");
	sprintf(bm->SP_FISHERYprmsName[q_id], "%s", "q_");
	sprintf(bm->SP_FISHERYprmsName[mFC_start_age_id], "%s", "_mFC_startage");
    sprintf(bm->SP_FISHERYprmsName[mFC_end_age_id], "%s", "_mFC_endage");
	sprintf(bm->SP_FISHERYprmsName[mFC_num_changes_id], "%s", "_mFC_changes");
	sprintf(bm->SP_FISHERYprmsName[mFC_id], "%s", "mFC_");
    sprintf(bm->SP_FISHERYprmsName[flagFchange_id], "%s", "flagFchange_");
	sprintf(bm->SP_FISHERYprmsName[flaghabitat_id], "%s", "flaghabitat_");
	sprintf(bm->SP_FISHERYprmsName[flagescapement_id], "%s", "flagescapement_");
    sprintf(bm->SP_FISHERYprmsName[spawn_closure_id], "%s", "spawn_closure_");
	sprintf(bm->SP_FISHERYprmsName[flagF_id], "%s", "flagF_");
	sprintf(bm->SP_FISHERYprmsName[Ka_escape_id], "%s", "Ka_escape_");
	sprintf(bm->SP_FISHERYprmsName[Kb_escape_id], "%s", "Kb_escape_");
	sprintf(bm->SP_FISHERYprmsName[flagdiscard_id], "%s", "flagdiscard_");
	sprintf(bm->SP_FISHERYprmsName[flagchangeDISCRD_id], "%s", "flagchangeDISCRD_");
	sprintf(bm->SP_FISHERYprmsName[DISCRD_num_changes_id], "%s", "_discard_changes");
	sprintf(bm->SP_FISHERYprmsName[FFCDR_id], "%s", "FFCDR_");
	sprintf(bm->SP_FISHERYprmsName[FC_thresh_id], "%s", "FC_thresh_");
	sprintf(bm->SP_FISHERYprmsName[FC_high_thresh_id], "%s", "FC_high_thresh_");
	sprintf(bm->SP_FISHERYprmsName[FCthreshli_id], "%s", "FCthreshli_");
	sprintf(bm->SP_FISHERYprmsName[k_retain_id], "%s", "k_retain_");
	sprintf(bm->SP_FISHERYprmsName[incidmort_id], "%s", "incidmort_");
	sprintf(bm->SP_FISHERYprmsName[k_waste_id], "%s", "k_waste_");
	sprintf(bm->SP_FISHERYprmsName[p_escape_id], "%s", "p_escape_");
	sprintf(bm->SP_FISHERYprmsName[TAC_num_changes_id], "%s", "_TAC_changes");
	sprintf(bm->SP_FISHERYprmsName[flagimposecatch_id], "%s", "flagimposecatch_");
	sprintf(bm->SP_FISHERYprmsName[imposecatchstart_id], "%s", "imposecatchstart_");
	sprintf(bm->SP_FISHERYprmsName[imposecatchend_id], "%s", "imposecatchend_");
	sprintf(bm->SP_FISHERYprmsName[FC_reportscale_id], "%s", "reportscale_");
	sprintf(bm->SP_FISHERYprmsName[trip_lim_id], "%s", "TripLimit_");
	sprintf(bm->SP_FISHERYprmsName[co_sp_catch_id], "%s", "co_sp_catch_");
	sprintf(bm->SP_FISHERYprmsName[co_sp_catch2_id], "%s", "co_sp_catch2_");
	sprintf(bm->SP_FISHERYprmsName[prop_spawn_close_id], "%s", "prop_spawn_close_");
	sprintf(bm->SP_FISHERYprmsName[TAC_id], "%s", "TAC_");
	sprintf(bm->SP_FISHERYprmsName[FC_case_id], "%s", "FC_case");
	sprintf(bm->SP_FISHERYprmsName[avail_id], "%s", "avail_");
	sprintf(bm->SP_FISHERYprmsName[saleprice_id], "%s", "saleprice");
	sprintf(bm->SP_FISHERYprmsName[tax_id], "%s", "tax");
	sprintf(bm->SP_FISHERYprmsName[deemedvalue_id], "%s", "deemed");
    sprintf(bm->SP_FISHERYprmsName[assess_nf_id], "%s", "assess_nf_");
    sprintf(bm->SP_FISHERYprmsName[flagPerShotCPUE_id], "%s", "flagPerShotCPUE_");
    sprintf(bm->SP_FISHERYprmsName[flagRecordCPUE_id], "%s", "flagRecordCPUE_");
    
	//#define catch_allowed 23				/* Catch currently allowed to take (quota - cumulative catch) */
	//#define flagquota_id 24
	//#define marketwgt_id 25			    /* Market weighting for amount of fish feed to each market */
	//#define desired_chrt_id 31
	//#define origprice_id 32
	//#define TACvsMPA_id 38				/* Whether to reduce q to reflect effects of seasonal spawning season closures */
	//#define phase_out_id 41				/* Day of run when TAC should be phased out completely */
	//#define phase_start_id 42				/* Day of run when TAC phase out begins */


	RBCParamStrings = (char **) c_alloc2d(25, num_rbc_species_params_id);
	sprintf(RBCParamStrings[DiscType_id], "%s", "DiscType_");
	sprintf(RBCParamStrings[MaxH_id], "%s", "MaxH_");
	sprintf(RBCParamStrings[Growthage_L1_id], "%s", "Growthage_L1_");
	sprintf(RBCParamStrings[Growthage_L2_id], "%s", "Growthage_L2_");
	sprintf(RBCParamStrings[MinCatch_id], "%s", "MinCatch_");
	sprintf(RBCParamStrings[Tier4_Cmaxmult_id], "%s", "Tier4_Cmaxmult_");
    sprintf(RBCParamStrings[AssessStart_id], "%s", "AssessStart_");
    sprintf(RBCParamStrings[NumRegions_id], "%s", "NumRegions_");
    sprintf(RBCParamStrings[Nsexes_id], "%s", "Nsexes_");

	sprintf(RBCParamStrings[Maturity_Inflect_id], "%s", "Maturity_Inflect_");
	sprintf(RBCParamStrings[Maturity_Slope_id], "%s", "Maturity_Slope_");
	sprintf(RBCParamStrings[Tier3_Fcalc_id], "%s", "Tier3_Fcalc_");
	sprintf(RBCParamStrings[Tier3_time_id], "%s", "Tier3_time_");
	sprintf(RBCParamStrings[Tier3_maxage_id], "%s", "Tier3_maxage_");
	sprintf(RBCParamStrings[T1_steep_phase_id], "%s", "T1_steep_phase_");
	sprintf(RBCParamStrings[tiertype_id], "%s", "tiertype_");
	sprintf(RBCParamStrings[Tier3_M_id], "%s", "Tier3_M_");
	sprintf(RBCParamStrings[Tier3_S25_id], "%s", "Tier3_S25_");
	sprintf(RBCParamStrings[Tier3_S50_id], "%s", "Tier3_S50_");
	sprintf(RBCParamStrings[Tier3_F_id], "%s", "Tier3_F_");
	sprintf(RBCParamStrings[Tier3_h_id], "%s", "Tier3_h_");
	sprintf(RBCParamStrings[Tier3_matlen_id], "%s", "Tier3_matlen_");
	sprintf(RBCParamStrings[Tier3_maxF_id], "%s", "Tier3_maxF_");
	sprintf(RBCParamStrings[Tier4_avtime_id], "%s", "Tier4_avtime_");
	sprintf(RBCParamStrings[Tier4_CPUEyrmin_id], "%s", "Tier4_CPUEyrmin_");
	sprintf(RBCParamStrings[Tier4_CPUEyrmax_id], "%s", "Tier4_CPUEyrmax_");
    sprintf(RBCParamStrings[Tier4_Bo_correct_id], "%s", "Tier4_Bo_correct_");    
	sprintf(RBCParamStrings[Tier4_m_id], "%s", "Tier4_m_");
	sprintf(RBCParamStrings[Tier4_alpha_id], "%s", "Tier4_alpha_");
	sprintf(RBCParamStrings[Tier5_length_id], "%s", "Tier5_length_");
	sprintf(RBCParamStrings[Tier5_S50_id], "%s", "Tier5_S50_");
	sprintf(RBCParamStrings[Tier5_cv_id], "%s", "Tier5_cv_");
	sprintf(RBCParamStrings[Tier5_flt_id], "%s", "Tier5_flt_");
	sprintf(RBCParamStrings[Tier5_reg_id], "%s", "Tier5_reg_");
	sprintf(RBCParamStrings[Tier5_p_id], "%s", "Tier5_p_");
	sprintf(RBCParamStrings[Tier5sel_id], "%s", "Tier5sel_");
	sprintf(RBCParamStrings[Tier5q_id], "%s", "Tier5q_");
	sprintf(RBCParamStrings[PostRule_id], "%s", "PostRule_");
    
	sprintf(RBCParamStrings[Tier1Sig_id], "%s", "Tier1Sig_");
	sprintf(RBCParamStrings[Tier2Sig_id], "%s", "Tier2Sig_");
	sprintf(RBCParamStrings[Tier3Sig_id], "%s", "Tier3Sig_");
    
    sprintf(RBCParamStrings[UseRBCAveraging_id], "%s", "UseRBCAveraging_");

    sprintf(RBCParamStrings[isTriggerSpecies_id], "%s", "isTriggerSpecies_");
	sprintf(RBCParamStrings[TriggerResponseScen_id], "%s", "TriggerResponseScen_");

	sprintf(RBCParamStrings[CPUEmult_id], "%s", "CPUEmult_");
	sprintf(RBCParamStrings[MaxChange_id], "%s", "MaxChange_");

	sprintf(RBCParamStrings[Hsteep_id], "%s", "Hsteep_");
	sprintf(RBCParamStrings[Agesel_Pattern_id], "%s", "Agesel_Pattern_");
	sprintf(RBCParamStrings[AssessFreq_id], "%s", "AssessFreq_");
	sprintf(RBCParamStrings[CCsel_years_id], "%s", "CCsel_years_");

	sprintf(RBCParamStrings[MG_offset_id], "%s", "MG_offset_");
	sprintf(RBCParamStrings[RecDevBack_id], "%s", "RecDevBack_");
	sprintf(RBCParamStrings[Regime_shift_assess_id], "%s", "Regime_shift_assess_");
    sprintf(RBCParamStrings[Regime_year_assess_id], "%s", "Regime_year_assess_");
    sprintf(RBCParamStrings[BallParkF_id], "%s", "BallParkF_");
	sprintf(RBCParamStrings[BallParkYr_id], "%s", "BallParkYr_");
	sprintf(RBCParamStrings[NumChangeLambda_id], "%s", "NumChangeLambda_");
    sprintf(RBCParamStrings[num_enviro_obs_id], "%s", "num_enviro_obs_");
    sprintf(RBCParamStrings[num_growth_morphs_id], "%s", "num_growth_morphs_");
    
	sprintf(RBCParamStrings[Nsex_samp_id], "%s", "Nsex_samp_");
	sprintf(RBCParamStrings[MaxAge_id], "%s", "MaxAge_");
	sprintf(RBCParamStrings[Nyfuture_id], "%s", "Nyfuture_");
	sprintf(RBCParamStrings[Nlen_id], "%s", "Nlen_");
	sprintf(RBCParamStrings[Lbin_id], "%s", "Lbin_");
    sprintf(RBCParamStrings[NumFisheries_id], "%s", "NumFisheries_");    
    sprintf(RBCParamStrings[thresh_mat_id], "%s", "thresh_mat_");
    sprintf(RBCParamStrings[femsexratio_id], "%s", "femsexratio_");

	sprintf(RBCParamStrings[flagLAdirect_id], "%s", "flagLAdirect_");
	sprintf(RBCParamStrings[flagSLAdirect_id], "%s", "flagSLAdirect_");
	sprintf(RBCParamStrings[flagWAdirect_id], "%s", "flagWAdirect_");

	sprintf(RBCParamStrings[SigmaR1_id], "%s", "SigmaR1_");
	sprintf(RBCParamStrings[SigmaR2_id], "%s", "SigmaR2_");
	sprintf(RBCParamStrings[SigmaR_future_id], "%s", "SigmaR_future_");
	sprintf(RBCParamStrings[PSigmaR1_id], "%s", "PSigmaR1_");
	sprintf(RBCParamStrings[Regime_year_id], "%s", "Regime_year_");
	sprintf(RBCParamStrings[RecDevMinYr_id], "%s", "RecDevMinYr_");
	sprintf(RBCParamStrings[RecDevMaxYr_id], "%s", "RecDevMaxYr_");
	sprintf(RBCParamStrings[RecDevFlag_id], "%s", "RecDevFlag_");
	sprintf(RBCParamStrings[AutoCorRecDev_id], "%s", "AutoCorRecDev_");
    
    sprintf(RBCParamStrings[LFSSlim_id], "%s", "LFSSlim_");
    sprintf(RBCParamStrings[AFSSlim_id], "%s", "AFSSlim_");
    sprintf(RBCParamStrings[NumSurvey_id], "%s", "NumSurvey_");
    sprintf(RBCParamStrings[NblockPattern_id], "%s", "NblockPattern_");
    sprintf(RBCParamStrings[SRBlock_id], "%s", "SRBlock_");
    sprintf(RBCParamStrings[assRecDevMinYear_id], "%s", "assRecDevMinYear_");
    
    sprintf(RBCParamStrings[MultispAssessType_id], "%s", "MultispAssessType_");
    sprintf(RBCParamStrings[mgt_indicator_id], "%s", "mgt_indicator_");
    sprintf(RBCParamStrings[init_mgt_category_id], "%s", "init_mgt_category_");
    sprintf(RBCParamStrings[init_mgt_sp_id], "%s", "init_mgt_sp_");
    
    return;
}

/*
 * \brief Get the index of the given tracer
 *
 * return -1 if not found.
 */
int Util_Get_Tracer_ID(MSEBoxModel *bm, char *tracerName) {
	int returnValue = -1;
	int i;

	for (i = 0; i < bm->ntracer; i++) {
		if (strcmp(bm->tinfo[i].name, tracerName) == 0)
			return i;
	}

	return returnValue;
}

/**
 * \brief Return the scaling factor that should be used at this timestep.
 *
 *
 */
double Util_Get_Change_Scale(MSEBoxModel *bm, int num_changes, double **changeArray) {

	int now_change = 0, past_change, i;
	double scale = 0.0, end_date, multB, multA, start, period;

	if ((bm->dayt >= changeArray[0][start_id]) && (changeArray[0][start_id] != 0)) {
		now_change = 0;
		for (i = 0; i < num_changes; i++) {
			if (bm->dayt >= changeArray[i][start_id])
				now_change = i;
		}

		end_date = (changeArray[now_change][start_id] + changeArray[now_change][period_id]);
		past_change = now_change - 1;
		multA = changeArray[now_change][mult_id];

		if (now_change > 0)
			multB = changeArray[past_change][mult_id];
		else
			multB = 1.0;
		start = changeArray[now_change][start_id];
		period = changeArray[now_change][period_id] + small_num;

		if (end_date < bm->dayt)
			scale = multA;
		else {
			if (multA > multB)
				scale = multB + (multA - multB) * (bm->dayt - start) / period;
			else
				scale = multB - (multB - multA) * (bm->dayt - start) / period;
		}
	} else {
		scale = 1.0;
	}
	return scale;
}

/**
 * \brief Return the scaling factor that should be used at this timestep.
 *
 *
 */
double Util_Get_Accumulative_Change_Scale(MSEBoxModel *bm, int num_changes, double **changeArray) {

	int changeIndex;
	double scale = 1.0, end_date, step1  = 0.0;

	for (changeIndex = 0; changeIndex < num_changes; changeIndex++) {

		if ((bm->dayt >= changeArray[changeIndex][start_id]) && (changeArray[0][start_id] != 0)) {

			end_date = changeArray[changeIndex][start_id] + changeArray[changeIndex][period_id];

			if (bm->dayt >= end_date) {
                /* If the increase has finished then simply apply */
				scale = scale * changeArray[changeIndex][mult_id];
			} else {
                /* Else apply partial outcome */
                if (changeArray[changeIndex][mult_id] < 1.0) {
                    step1 = (1.0 - changeArray[changeIndex][mult_id]) * ((bm->dayt - changeArray[changeIndex][start_id]) / changeArray[changeIndex][period_id]);
                    scale = scale * (1.0 - step1);
                } else {
                    step1 = (changeArray[changeIndex][mult_id] - 1.0) * ((bm->dayt - changeArray[changeIndex][start_id]) / changeArray[changeIndex][period_id]);
                    scale = scale * (step1 + 1.0);
                }
			}
		}
	}
	return scale;
}

void Util_Copy_Change_Values(MSEBoxModel *bm, double **originalArray, double **newArray, int size, double start, double period, double mult) {

	int arrayIndex, i;
	int beenInserted = FALSE;

	/* Now copy the values across inserting the new value in the correct place.*/
	arrayIndex = 0;
	if (size > 1) {
		for (i = 0; i < size - 1; i++) {
			if (beenInserted == FALSE) {
				if (originalArray[i][start_id] > start) {
					/* insert the new values into the array here */
					newArray[arrayIndex][start_id] = start + bm->dayt;
					newArray[arrayIndex][period_id] = period;
					newArray[arrayIndex][mult_id] = mult;
					arrayIndex++;
					beenInserted = TRUE;
				}
			}
			newArray[arrayIndex][start_id] = originalArray[i][start_id];
			newArray[arrayIndex][period_id] = originalArray[i][period_id];
			newArray[arrayIndex][mult_id] = originalArray[i][mult_id];

			arrayIndex++;
		}
	}
	/* Make sure the value has been inserted else insert it at the end */
	if (beenInserted == FALSE) {
		newArray[arrayIndex][start_id] = start + bm->dayt;
		newArray[arrayIndex][period_id] = period;
		newArray[arrayIndex][mult_id] = mult;
	}
}

/**
 * \brief Get the current stock id for this species in the given box and depth.
 *
 * Only age structured groups have stock structure typically so the stocj array
 * for other groups will be zeroed on initialisation. So to check against that
 * use max_id to make sure there is a stock map, if not default to stock_id 0
 */
int Util_Get_Current_Stock_Index(MSEBoxModel *bm, int sp, int boxNumber, int nzk) {
	int sp_stock_type;
	int stock_id = 0;
	int relk, diffdz;


	diffdz = bm->wcnz - bm->boxes[boxNumber].nz;
	relk = nzk + diffdz;
	sp_stock_type = (int) (FunctGroupArray[sp].speciesParams[stockstruct_type_id]);
	switch (sp_stock_type) {
	case horiz_only:
		stock_id = bm->stock_struct[boxNumber][sp] - 1;
		break;
	case vert_only:
		stock_id = bm->v_stock_struct[relk][sp] - 1;
		break;
	case mixed_stock:
		stock_id = (bm->stock_struct[boxNumber][sp] - 1) + (bm->v_stock_struct[relk][sp] - 1) * ((int) (FunctGroupArray[sp].speciesParams[hstockstruct_id]));
		break;
	default:
		quit("No such stock structure case (%d) defined. Recode\n");
		break;
	}



	if(stock_id < 0) {
		// Default to a value of zero - will be for biomass pools (age structured or otherwise)
		stock_id = 0;
		//quit("Stock_id %d is not a valid value for %s (stock type: %d)", stock_id, FunctGroupArray[sp].groupCode, sp_stock_type);
	}

	if(stock_id != bm->group_stock[sp][boxNumber][nzk]){
			printf("bm->group_stock[%d][%d][%d] = %d,  stock_id = %d\n", sp, boxNumber, nzk, bm->group_stock[sp][boxNumber][nzk], stock_id);
			abort();
		}
	return stock_id;
}

/**
 * 
 * Calculate the stock id for each group in each box/layer combination.
 */
void Util_Calculate_StockID(MSEBoxModel *bm){
	int sp, box, layer;
	int stock_type, diff, relk, stock = 0;

	for(sp = 0; sp < bm->K_num_tot_sp; sp++){
		stock_type = (int) (FunctGroupArray[sp].speciesParams[stockstruct_type_id]);
		for(box = 0; box < bm->nbox; box++){
			diff = bm->wcnz - bm->boxes[box].nz;

            for(layer = 0; layer < bm->boxes[box].nz; layer++){
				relk = layer + diff;
				switch (stock_type) {
				case horiz_only:
					stock = bm->stock_struct[box][sp] - 1;
					break;
				case vert_only:
					stock = bm->v_stock_struct[relk][sp] - 1;
					break;
				case mixed_stock:
					stock = (bm->stock_struct[box][sp] - 1) + (bm->v_stock_struct[relk][sp] - 1)
							* ((int) (FunctGroupArray[sp].speciesParams[hstockstruct_id]));
					break;
				default:
					quit("No such stock structure case (%d) defined. Recode\n", stock_type);
					break;
				}
				if(stock < 0) {
					// Default to a value of zero - will be for biomass pools (age structured or otherwise)
					stock = 0;
					//quit("Stock_id %d is not a valid value for %s (stock type: %d)", stock_id, FunctGroupArray[sp].groupCode, sp_stock_type);
				}

				bm->group_stock[sp][box][layer] = stock;
			}
		}
	}

}

/******** HELPER ROUTINES ************************************************************

 These general routines are used by all the modules, mainly to do with feeding
 or are other helper routines.
 **/

/**
 *	\brief Michaelis-Menten Kinetics
 *
 */
double Util_Mich_Ment(double x, double m) {
	double result_MM = 0;

	if (x > 0.0)
		if (m > 0.0)
			result_MM = x / (m + x);
		else
			result_MM = 1.0;
	else
		result_MM = 0.0;

	return (result_MM);
}

/**
 *    \brief random number from normal distribution
 *
 */
double Util_Normx_Result(double mu, double sigma) {
    double result, step_a, step_b, step_c;

    step_a = drandom(0.0, 1.0);
    step_b = drandom(0.0, 1.0);
    step_c = sqrt(-2.0 * log(step_a)) * cos(2.0 * 3.1415926 * step_b);
    result = step_c * sigma + mu;

    return (result);

}

/**
 *	\brief random number from lognormal distribution
 *
 */
double Util_Logx_Result(double mu, double sigma) {
	double x_b, result, step_a, step_b, step_c;

	x_b = drandom(0.0, 1.0);
	step_a = 1.0 / (x_b * sigma * sqrt(2.0 * 3.141592654));
	step_b = (log(x_b) - mu) * (log(x_b) - mu);
	step_c = 2.0 * sigma * sigma;
	result = step_a * exp(-step_b / step_c);

	return (result);

}

/**
 *	\brief lognormal distribution
 *
 */
double Util_Lognorm_Distrib(double mu, double sigma, double x_b) {
	double result, step_a, step_b, step_c;

	step_a = 1.0 / (sigma * sqrt(2.0 * 3.141592654));
	step_b = (log(x_b) - mu) * (log(x_b) - mu);
	step_c = 2.0 * sigma * sigma;
	result = step_a * exp(-step_b / step_c);

	return (result);

}

int at_compileRegExpression(regex_t *regBuffer, char *str) {
	int returnValue;
	regBuffer->re_nsub = 1;
	returnValue = regcomp(regBuffer, str, REG_EXTENDED);

	if (returnValue != 0) {
		printf("error compiling regular expression");
		return -1;
	}
	return 0;
}


/**
 *	Open an Atlantis output file in the given destination folder.
 *
 */
FILE	*Util_fopen(MSEBoxModel *bm, const char *name, const char *mode){
	char fileName[BMSLEN];


	sprintf(fileName, "%s%s", bm->destFolder, name);
	return fopen(fileName,mode);
}


/**
 *	\brief - Wrapper around the ncopen function to check that a file with the given filename
 *	exists before ncopen is call.
 *
 *	@returns The FILE * pointer to the opened file.
 *
 *	@param name The name of the file to open
 *	@param mode The netcdf mode to use when opening the file. See the ncopen documentation
 *	for more information.
 *
 *
 */
int Util_ncopen(const char *destFolder, const char* name ,int mode){
	int fp;
	FILE *fileFp;
	char fileName[BMSLEN];

	if(strlen(destFolder) > 0){
		sprintf(fileName, "%s%s", destFolder, name);
	}else{
		strcpy(fileName, name);
	}
	trim(fileName);
    
    printf("Opening %s\n", fileName);

	/**
	 *	Open the file to check that it exists if the mode is NC_NOWRITE.
	 *	The netcdf library should return an error if the file
	 *	does not exits but its spitting the dummy.
	 *
	 */
	if(mode == NC_NOWRITE){
		if((fileFp = fopen(fileName, "r")) == NULL)
			quit("fopen: Can't open netcdf model input data fileno readable - reason 1 \n '%s' \n",fileName);
		fclose(fileFp);
		/* Open the file */
		if( (fp=ncopen(fileName,mode)) < 0 )
		   quit("Util_ncopen: Can't open netcdf model input data file - ncopen fail '%s'\n",fileName);
	}
	else{
		if( (fp=ncopen(fileName,mode)) < 0 ){
			return -1;
		}
			//quit("Util_ncopen: Can't create netcdf model input data file %s\n",fileName);
	}
    
    if(verbose) {
        printf("Finished Util_ncopen\n");
    }
    
	return fp;
}

//  mean is mean of normal dsn and sigg is std devn
//
//  Calls :   ran1  (which is in cplib library)
//  Called by:
//  Created: translated from Gavin's Fortran by Sally 10/10/2007
//
// ******************************************************************************

double Util_xnorm(double mean, double sigg, int *iiseed)
{
	double z1,z2,xn;

// ignore zero calls
	if (sigg==0.0) return 0.0;

	z1 = ran3(iiseed);
	z2 = ran3(iiseed);

	xn = sin(6.283185*z1) * sqrt(-2.0*log(z2)) * sigg + mean;

	return xn;

}

/******************************************************************************
 *
 *  Name: GenMnorm
 *  Description: generate from a multivariate normal
 *  Parameters :
 *      vec is the answer
 *      means is the vector of means
 *      np is number of parameters (array dimensions)
 *      tt is np x np correlation matrix
 *      sg is 1 x np vector - standard deviation of random walk
 *
 *  Calls :    xnorm
 *  Created:  translated from Gavin's Fortran by Sally 10/10/07
 *
 ******************************************************************************/
void Util_GenMnorm(double *vec, double *means, int *iseed, int np, double **tt, double *sg){
	double *aa;
	double x,sigs;
	int i,j;

	// initial errors to zero
	aa = Util_Alloc_Init_1D_Double(np, 0.0);

	// generate the multivariate normal - sg is the same for all lengths
	for (i=0; i < np; i++){
		sigs = sg[i];
		x = Util_xnorm(0.0, sigs, iseed);
		for (j=0; j < np; j++)
			aa[j] += x * tt[j][i];    // shouldn't this be sqrt(tt(j,i))?? NO correlation, not var-covar
	}

	// complete the generation
	for (i=0; i < np; i++)
		vec[i] = means[i] + aa[i];

	// clean-up
	free1d(aa);

	return;
}

int Util_Check_NetCDF_Size(MSEBoxModel *bm, int fid, int *dump, char *fileName, int *index, int type){
	struct stat buffer;
	char fname[STRLEN*2];
	char endname[] = ".nc";
	char tempStr[STRLEN];
	char *pdest;
	int newFid;

	if(fid < 0)
		return -1;

	if((*index) > 0){
		strcpy(tempStr, fileName);
		pdest = strstr(tempStr, endname);
		*pdest = '\0';

		sprintf(fname, "%s%s%d.nc", bm->destFolder, tempStr, *index);
	}else{
		sprintf(fname, "%s%s", bm->destFolder, fileName);
	}


	/* Might need to close the file and then open it after we have checked its size */
	/* Check the size of the diff file. */
	if(stat(fname, &buffer) < 0){
		quit("Unable to get information about file %s\n", fileName);
	}

	if (buffer.st_size > MAX_NETCDF_FILE_SIZE) {
		strcpy(tempStr, fileName);
		pdest = strstr(tempStr, endname);
		*pdest = '\0';

		(*index)++;
		sprintf(fname, "%s%d.nc", tempStr, *index);
		ncclose(fid);

		/* Create file anew */
		if(type == -1){
			newFid = createBMSummaryDataFile(bm->destFolder, fname, bm);
		} else if (type == 4) {
            newFid = createBMAnnAgeBioDataFile(bm->destFolder, fname, bm);
        } else if (type == 5) {
            newFid = createBMAnnAgeCatDataFile(bm->destFolder, fname, bm);
        } else if (type == 6) {
            newFid = createBMDietDataFile(bm->destFolder, fname, bm);
        } else {
			newFid = createBMDataFile(bm->destFolder, fname, bm, type);
		}
		*dump = 0;


		return newFid;

	}
	return fid;
}
