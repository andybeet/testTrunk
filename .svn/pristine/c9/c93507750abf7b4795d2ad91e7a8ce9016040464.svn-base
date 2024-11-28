/**
 \file
 \brief C file for economic models for Atlantis
 \ingroup ateconomic

 File:	ateconomic.c
 Author:	Beth Fulton
 Created:	27/10/2004
 Purpose:	C file for economic models for Atlantis
 Revisions: 27/10/2004 Created the file

 29/9/2005 Added Simple_Econ_Indicators() routine

 2/10/2005 Added structure of Economics() routine

 14-04-2008 Beth Fulton
 Moved bm->SUBFLEET_ECONprms weekly reinitialisation into the ns loop

 28-10-2009 Bec Gorton
 Removed references to the K_num_fished_sp and K_num_impacted_sp. Arrays that were this
 size have been resized to K_num_tot_sp.

 04-11-2009 Bec Gorton
 Merged in Beths new bycatch incentive code - revision 961.
 **/

/*  Open library routines */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sjwlib.h>
#include <atEconomic.h>

void Economic_Output_Temporal_Catches(MSEBoxModel *bm);
void Economic_Output_ExpectedCatchPerMonth(MSEBoxModel *bm);

/**
 *
 *	\brief Boxmodel annual economic routines - this is called each year.
 *
 */
void Economic_Annual(MSEBoxModel *bm, FILE *llogfp) {
	int nf, ns, sp, regid;

	if (verbose > 1)
		printf("Doing annual economics\n");

	/** Get current fish prices - this way initialised before needed at start of year and run **/
	Economic_Response(bm, llogfp);
	Get_Fish_Prices(bm, llogfp);

	/** Update Economic Indices and Output them **/
	Update_Econ_Indx(bm, llogfp);
	Update_Port_Indx(bm, llogfp);
	Economic_Output_Indices(bm, llogfp);

	fprintf(llogfp, "Day %e, output econ measure\n", bm->dayt);

	Economic_Output_Meas(bm);
	Economic_Output_Port_Meas(bm);

	/** Quota allocations **/
	/* Need to know what quota you have to make your effort decisions */
	Annual_Quota_Allocation(bm, llogfp);

	/** Annual scale economic decisions - also write out econ indicators at this scale **/
	Annual_Effort_Schedule(bm, llogfp);

	/* Reinitialise indices for the new year */
	for (nf = 0; nf < bm->K_num_fisheries; nf++) {
		bm->econindx[nf][fleet_switch_id] = 0;
		bm->econindx[nf][overallrent_id] = 0;
		bm->econindx[nf][SubfleetCash_id] = 0;
		bm->econindx[nf][LeaseVal_id] = 0;
		bm->econindx[nf][InvestReturn_id] = 0;
		bm->econindx[nf][GrossVal_id] = 0;
		bm->econindx[nf][GrossTax_id] = 0;
		bm->econindx[nf][GrossDV_id] = 0;

		/* Initialise regional cumulative catch */
		for (ns = 0; ns < bm->FISHERYprms[nf][nsubfleets_id]; ns++) {

			if (!bm->MultiPlanEffort) {
				/* Only do this for Dan Holland code here, all else do it in Update_Vessel_Numbers() */
				bm->SUBFLEET_ECONprms[nf][ns][lasttot_cash_id] = bm->SUBFLEET_ECONprms[nf][ns][tot_cash_id];
				bm->SUBFLEET_ECONprms[nf][ns][tot_cash_id] = 0;
			}
			bm->SUBFLEET_ECONprms[nf][ns][gvp_ind_id] = 0;
			bm->SUBFLEET_ECONprms[nf][ns][tax_ind_id] = 0;
			bm->SUBFLEET_ECONprms[nf][ns][dv_ind_id] = 0;
			for (sp = 0; sp < bm->K_num_tot_sp; sp++) {
				if (FunctGroupArray[sp].isFished == TRUE) {
					bm->QuotaAlloc[nf][ns][sp][DVtriggered_id] = 0;
					bm->QuotaAlloc[nf][ns][sp][TAXtriggered_id] = 0;
					bm->QuotaAlloc[nf][ns][sp][deemed_value_id] = 0;
					for (regid = 0; regid < bm->K_num_active_reg; regid++) {
						RegCatch[nf][ns][regid][sp] = 0;
					}
				}
			}

			tot_marg_rent[nf][ns] = 0;

		}
	}

	return;
}

/**
 *
 *	\brief Boxmodel economic routines - this is called each timestep.
 *
 *	Harvest code doesn't really have the subfleets. So each time info is passed between
 *	harvest and economics they need to be aggregates.
 *
 */
void Economics(MSEBoxModel *bm, FILE *llogfp) {
	int nf, ns, sp;
	int trip_on = 1;

	if (verbose > 1)
		printf("Doing economics\n");

	/* Divide current catch amongst subfleets - not needed in first timestep
	 as done in initialisation in initEconomic() */
	if (bm->dayt > 0)
		Allocate_Catch(bm, llogfp);
	else
		fprintf(bm->logFile, "catch not allocated yet \n");

	/** Coarse scale economic decisions
	 Roughly monthly level calculations, some part of quota calulations, but
	 also seasonal-level effort allocation decisions
	 **/

	/* Call indicator calculations */
	if ((!bm->TofY || bm->newmonth) && bm->flagday) {
		Economic_Indicators(bm, llogfp);

		/* Update other economic and social concerns - not needed in first month as already done
		 in Economic_Annual()
		 */
		if (bm->TofY)
			Economic_Response(bm, llogfp);

		Update_Econ_Indx(bm, llogfp);

		bm->reset_month_val = 1;
	}
	/* Allocate effort accordingly  - where are we at against our plan for the year.*/
	if (bm->newmonth && bm->flagday)
		Monthly_Effort_Schedule(bm, llogfp);

	/** Fine scale economic decisions
	 Quota trading happens at this scale, time scale is typically about one week
	 **/
	if (!bm->econweekly) {
		/* If effort allocation has non-weekly update (e.g. as per Dan Holland)
		 then have non-week updating capability. Check done for updating within
		 the actual routine now as may vary between fleets.
		 */
		trip_on = 1;
	} else {
		/* For the base economic model the update is assumed to happen once per week,
		 which is caught by the bm->newweek check and so set trip_on = 0 here so only
		 tripped once per week in that case.
		 */
		trip_on = 0;
	}

	//	fprintf(llogfp,"TIme: %e, trip_on: %d, bm->newweek: %d, bm->flagday: %d\n",
	//		bm->dayt, trip_on, bm->newweek, bm->flagday);

	if ((trip_on || bm->newweek) && bm->flagday) {
		/* If using ITE then trade effort */
		if (bm->TradeEffort)
			Effort_Quota_trade(bm, llogfp);

		/* Trade catch quota */
		if(bm->quota_trading)
			Quota_trade(bm, llogfp);
        else
            Consolidate_Allowed_Catches(bm, llogfp);  // Or trip limits never reset

		/** Calculate final realised effort allocation - distribute the spatial effort.
		 * Aggregate the subfleet data to pass to the harvest module.
		 **/
		Allocate_Final_Effort(bm, llogfp);
	}

	/* Initialise weekly catch records if necessary - do this after done all other calculations so
	 always working with a week of fisheries data in the economics (if did it at the start
	 would only have one days worth, which throws all the cost calculations out as they're based around
	 trips, which are typically a week or more
	 */
	if (bm->newweek) {
		for (nf = 0; nf < bm->K_num_fisheries; nf++) {
			for (ns = 0; ns < bm->FISHERYprms[nf][nsubfleets_id]; ns++) {
				for (sp = 0; sp < bm->K_num_tot_sp; sp++) {
					if (FunctGroupArray[sp].isFished == TRUE) {
						bm->QuotaAlloc[nf][ns][sp][boatcatch_id] = 0;
					}
				}
				bm->SUBFLEET_ECONprms[nf][ns][discards_ind_id] = 0;
			}
		}
	}
	if (bm->newmonth && bm->reset_month_val) {

		/* Write out monthly catches */
		Economic_Output_Temporal_Catches(bm);
		if(bm->dayt > 0) /* will actually write out data from the previous month */
			Economic_Output_ExpectedCatchPerMonth( bm );

		/* Reinit the arrays */
		for (nf = 0; nf < bm->K_num_fisheries; nf++) {
			for (ns = 0; ns < bm->FISHERYprms[nf][nsubfleets_id]; ns++) {
				for (sp = 0; sp < bm->K_num_tot_sp; sp++) {
					if (FunctGroupArray[sp].isFished == TRUE) {
						bm->QuotaAlloc[nf][ns][sp][cummonthcatch_id] = 0;
					}
				}
				for (sp = 0; sp < bm->K_num_tot_sp; sp++) {
					if (FunctGroupArray[sp].isImpacted == TRUE) {
						bm->QuotaAlloc[nf][ns][sp][cummonthbycatch_id] = 0;
					}
				}
				bm->SUBFLEET_ECONprms[nf][ns][CurrentCatch_id] = 0;
				bm->SUBFLEET_ECONprms[nf][ns][CurrentEffort_id] = 0;
			}
		}
		bm->reset_month_val = 0;

		if (((bm->debug == debug_econeffort) || (bm->debug == debug_econ)) && (bm->dayt > bm->checkstart))
			fprintf(llogfp, "Time: %e Reset cummonthcatch\n", bm->dayt);
	}

	return;
}
