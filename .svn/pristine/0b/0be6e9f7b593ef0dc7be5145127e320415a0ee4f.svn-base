/**
 *  \ingroup atManageLib
 *	\file atManagePrivate.h
 *   Definitions of the variables that are private to the atManageLib.
 *
 *   Author: Beth Fulton
 *   CSIRO Marine Laboratories
 *   2004
 *
 *   Revisions: 8/8/2004 Created from the sebiolsetup.h file
 *
 *				16/8/2004 Cut back to include only those integers and doubles actually
 *                 used in msemanagesetup.c
 *
 *              21/9/2004 Renamed FSL, FSO, FSF, FST and FSG to generic names so apply to
 *				   other models/generic case.
 *				   FSL -> FDC  (Ling -> Demersal_DC_Fish)
 *				   FSO -> FDO  (Orange_Roughy -> Demersal_O_Fish)
 *				   FSF -> FDB  (Flathead -> Demersal_B_Fish)
 *				   FST -> FVO  (Bluefin_Tuna -> Pisciv_O_Fish)
 *				   FSG -> SHB  (Gummy_Shark -> Shark_B)
 *
 *              6/10/2004 Added parameter initalisation for new fisheries (ptrawlPWN,
 *				   dtrawlFBP, midwcZL, trapFDE, dlineFDE, netFDE) and biological
 *				   groups (BML, BFF, ZL, PWN, FBP, FDE, FDF, FVB, SSK, SHR, SP, REP, WDG)
 *
 *				23/10/2004 Added midwcPWN and additional fisheries links
 *
 *				18/1/2005 Modfified age specific selectivity names so now just one
 *					addults and one juvenile entry for each group (rather than a fishery
 *					entry for each group)
 *
 *				8/6/2005 Added mowMA and the 3 groups FDM, FDP, WHS.
 *
 *				18/03/2008 Bec Gorton
 *				Removed all the now un-used variables.
 *
 *				22-04-2009 Bec Gorton
 *				Removed the ports_popnumchangei and ports_popchangei arrays.
 *
 *				28-10-2009 Bec Gorton
 *				Removed reference to array that are no longer global but are now in the atManageModule.
 *				Also removed the temp arrays that are used in param loading.
 */


/******************************************************************
Modelling variables for control of processes within the model
*/

int *flagdropeffort = 0, *MPAKeyMap = 0, *checkedbox = 0, need_discard;

double ***p_fishi = 0, **k_cover = 0, ***Effort_vdistrib = 0, ****p_fish_origi = 0,
	***EFFORTchange = 0, ***qSTOCK = 0,
	**FC_hdistrib = 0, **MPAendangered = 0, **SEASONAL = 0, **effort_scale = 0,
	*oldFishEndDay = 0, *scale_effort = 0, *prev_mult = 0,
    **FC_case, *DistPeak = 0, *FrefAi = 0, *FrefHi = 0, *FreStarti = 0, *LeverUsei = 0,
    *estErrori = 0, *estCVi = 0, *estBiasi = 0;


double **mEff;
double **gear_conflict;
double **origEffort_vdistrib;
double *whichrefi;
double ***CAPchange;
double ***EFFORTchange;

double ***MPAendangeredlist;        /**< List of potential zoning (spatial management) based on
								           PET species needs */
double **MPAoverfishedlist;         /* List of potential zoning (spatial management) based on
								           overfished species (for US-like system) */
double ***LeverQueue;				/* MSY-based management queue */

double ***cpue_cdf;                  /* Cumulative distributions for the probabilty distribution associate
                                     with catch in CPUE generation - typically a negative binomial
                                     filled during set-up */
double *tot_alloced;



