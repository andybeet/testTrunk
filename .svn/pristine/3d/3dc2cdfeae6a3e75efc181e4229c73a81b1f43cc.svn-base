/**
 \file
 \brief C file for trade of effort quota in economic models for Atlantis
 \ingroup ateconomic

 File:	atquota.c
 Author:	Beth Fulton
 Created:	7/10/2005
 Purpose:	C file for economic effort quota trade models for Atlantis
 Revisions: 3/8/2009 Created the file

 **/

/*  Open library routines */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sjwlib.h>
#include <atEconomic.h>

/**
 *
 *	\brief Get overall effort quota prices
 *
 */
void Total_Effort_Quota_Price(MSEBoxModel *bm, FILE *llogfp) {
	int nf, ns;
	if (verbose > 1)
		printf("Calculating effort quota prices\n");

	/* Using Rich Little's lease price proxy */
	for (nf = 0; nf < bm->K_num_fisheries; nf++) {
		for (ns = 0; ns < bm->FISHERYprms[nf][nsubfleets_id]; ns++) {
			/* If no boats in the subfleet currently skip ahead */
			if (!bm->SUBFLEET_ECONprms[nf][ns][nboat_id]) {
				continue;
			}

			bm->SUBFLEET_ECONprms[nf][ns][effortquotaprice_id]
					= max(0.0,bm->SUBFLEET_ECONprms[nf][ns][margprofit_id] * (1.0 + 1.0 / (bm->interest_rate + small_num)));

			/* Cumulative lease cost through time */
			bm->SUBFLEET_ECONprms[nf][ns][eleased_Q_val_ind_id] += bm->SUBFLEET_ECONprms[nf][ns][effortquotaprice_id];
		}
	}

	return;
}

/**
 *
 *	\brief Effort quota trading model - based on Rich Little's quota model
 *
 */
void Effort_Quota_trade(MSEBoxModel *bm, FILE *llogfp) {
	//int do_debug_base = 0;
	//int tnf = bm->K_num_fisheries;
	//int maxmonth = 12;
	//int nm = maxmonth - bm->MofY;
	//int ascendflag = 0;
	//int totpi = bm->K_num_subfleet * bm->K_num_fisheries * maxmonth;
	//int checkmonth = 0;
	//int nsort = 0;

//	if (bm->debug == debug_quota)
//		do_debug_base = 1;
//	else
//		do_debug_base = 0;

	if (verbose > 1)
		printf("Trading effort quota\n");

	//	/** Trading **/
	//
	//	/* Following work is based on Rich Little's quota trading model *
	//
	//	for(nf=0; nf<bm->K_num_fisheries; nf++){
	//		spp_id = bm->SPid[sp];
	//		quota_sp = 0;
	//		aggDemand = 0;
	//		aggSupply = 0;
	//		within_min = 0;
	//		any_spare = 0;
	//
	//		/* Initialise linearPI *
	//		for(ns=0; ns<totpi; ns++){
	//			linearPI[nf][linPI_id][ns] = 0;
	//			linearPI[nf][idd_id][ns] = 0;
	//			ax[ns] = 0;
	//			bx[ns] = 0;
	//		}
	//
	//		idp = 0;
	//		if(bm->FISHERYprms[nf][flageffortmodel_id] > rec_econ_model){ /* Only commerical fisheries trade quota */
	//			/* Only bother continuing if an ITE present *
	//			if(bm->FISHERYprms[nf][flagcap_id]){
	//				for(ns=0; ns<bm->FISHERYprms[nf][nsubfleets_id]; ns++){
	//					/* If no boats in the subfleet currently skip ahead *
	//					if(!bm->SUBFLEET_ECONprms[nf][ns][nboat_id]){
	//						continue;
	//					}
	//					TripLength = bm->SUBFLEET_ECONprms[nf][ns][max_trip_length_id];
	//					TripCost = bm->SUBFLEET_ECONprms[nf][ns][cost_ind_id];
	//
	//
	//					// FROM HERE GET CPUE PER MONTH AND PUT MOST EFFORT IN MONTH WITH HIGHEST CPUE
	//
	//
	//
	//						/* Identify if anyone actually needs to trade this species rightly now *
	//						if(bm->QuotaAlloc[nf][ns][sp][within_id])
	//							within_min = 1;
	//
	//						/* Identify if anyone actually has some to spare *
	//						if(bm->QuotaAlloc[nf][ns][sp][spare2sell_id])
	//							any_spare = 1;
	//
	//						/* Calculate profit (pi) *
	//						cumCatch = bm->QuotaAlloc[nf][ns][sp][cumboatcatch_id];
	//						prop_olease = bm->QuotaAlloc[nf][ns][sp][permolease_id];
	//						prop_tlease = bm->QuotaAlloc[nf][ns][sp][templease_id];
	//						ownQuota = (1.0 - prop_olease - prop_tlease) * bm->QuotaAlloc[nf][ns][sp][owned_id];
	//						leaseQuota = bm->QuotaAlloc[nf][ns][sp][leased_id];
	//						priceQuota = bm->QuotaAlloc[nf][ns][sp][spmarg_profit_id] * (1.0 + 1.0 / (bm->interest_rate + small_num));
	//						for(month=bm->MofY; month<12; month++){
	//							expectedCatch = bm->MonthAlloc[nf][ns][sp][month];
	//							scheduledEffort = bm->EffortSchedule[nf][ns][month][expect_id];
	//							fishprice = bm->SP_FISHERYprms[spp_id][nf][saleprice_id];
	//							expectedCPUE = bm->BlackBook[nf][ns][sp][bm->MofY][expect_id] / (bm->EffortSchedule[nf][ns][bm->MofY][hist_id] + small_num);
	//							pi[nf][ns][month][sp] = scheduledEffort	* ((fishprice * expectedCPUE) - (TripCost / (TripLength + small_num)));
	//							pi[nf][ns][month][sp] -= priceQuota * (cumCatch + expectedCatch - ownQuota - leaseQuota);
	//
	//							idd = ns*tnf*nm + nm*nf + month - bm->MofY;
	//							linearPI[sp][linPI_id][idd] = bm->QuotaAlloc[nf][ns][sp][finalutility_id];
	//							linearPI[sp][idd_id][idd] = idd;
	//
	//							if(bm->QuotaAlloc[nf][ns][sp][within_id] || bm->QuotaAlloc[nf][ns][sp][spare2sell_id]){
	//								/* To keep sort times to a minimum only include those
	//									fishery components that are actually buying or selling
	//									at this time.
	//								*
	//								ax[idp] = linearPI[sp][linPI_id][idd];
	//								bx[idp] = linearPI[sp][idd_id][idd];
	//								cx[idp] = 1.0;		// Not needed in this case
	//								idp++;
	//							}
	//						}
	//
	//						/** Calculate aggregate demand and supply **/
	//						/* Aggregate demand - is sum of the differences between the cumulative catch
	//						+ expected future catch and the quota held already *
	//						aggDemand += cumCatch + expectedCatch - ownQuota - leaseQuota;
	//
	//						/* Aggregate supply - here purchase quota is sold quota so need is negative *
	//						aggSupply += ownQuota + leaseQuota - cumCatch;
	//
	//					}
	//				}
	//			}
	//		}
	//
	//		/* Only continue if quota species that someone has nearly filled their quota on
	//			and someone else has some spare to sell. Also don't bother continuing if not
	//			using the sp_by_sp code. More complex case is handled below.
	//		*
	//		if(!quota_sp || !within_min || !any_spare || !bm->sp_by_sp)
	//			continue;
	//
	//		/* Sort participating subfleet and fisheries based on marginal rent only - set nsort first
	//			(so don't waste time sorting non-participating members *
	//		nsort = idp;
	//		Quicksort_Dir(ax, bx, cx, dx, ex, nsort, ascendflag);
	//
	//		/* Copy ax and bx across to linearPIP for use below *
	//		for(ns=0; ns<totpi; ns++){
	//			linearPI[sp][linPI_id][ns] = ax[ns];
	//			linearPI[sp][idd_id][ns] = bx[ns];
	//		}
	//
	//		b = 0;
	//		s = nsort;
	//		while((s >= 0) && (b < nsort && b >= 0)){
	//			b = 0;
	//			demand = 1;
	//			while((demand > 0) && (b < nsort && b >= 0)){
	//				demand = 0;
	//				/* Deconstructing ids from single sorted array - assuming ids in
	//				the sorted array are of the form ns*tnf*(12 - month) + (12 - month)*nf + month - bm->MofY *
	//				step1 = (int)floor(ROUNDGUARD + (linearPI[sp][idd_id][b] / (nm + small_num)));	/* First step in finding fishery_id and subfleet_id *
	//				buyerns = (int)floor(ROUNDGUARD + (step1 / (tnf + small_num)));				/* Finding subfleet_id *
	//				buyernf = step1 % tnf;								/* Finding the fishery_id *
	//				buyermonth = (int)floor(ROUNDGUARD + (linearPI[sp][idd_id][b])) % nm;		/* Finding the month_id *
	//				down_time = bm->SUBFLEET_ECONprms[buyernf][buyerns][down_time_id];
	//				max_month_effort = bm->month_scalar * (1 - down_time) * bm->SUBFLEET_ECONprms[buyernf][buyerns][nboat_id];
	//
	//				/* The expected increase in catch by filling the current month *
	//				if((max_month_effort - bm->EffortSchedule[buyernf][buyerns][buyermonth][expect_id] >= 1) && (bm->MargRent[buyernf][buyerns][sp][buyermonth] >= 0)){
	//					scheduledEffort = bm->EffortSchedule[buyernf][buyerns][buyermonth][expect_id];
	//					ExpCatch = bm->BlackBook[buyernf][buyerns][sp][buyermonth][expect_id];
	//					ExpEffort = bm->EffortSchedule[buyernf][buyerns][buyermonth][hist_id];
	//					expectedIncreasedCatch = (max_month_effort - scheduledEffort) * ExpCatch / (ExpEffort + small_num);
	//					expectedCatch = 0;
	//					prop_olease = bm->QuotaAlloc[buyernf][buyerns][sp][permolease_id];
	//					prop_tlease = bm->QuotaAlloc[buyernf][buyerns][sp][templease_id];
	//					ownQuota = (1.0 - prop_olease - prop_tlease) * bm->QuotaAlloc[buyernf][buyerns][sp][owned_id];
	//					leaseQuota = bm->QuotaAlloc[buyernf][buyerns][sp][leased_id];
	//					cumCatch = bm->QuotaAlloc[buyernf][buyerns][sp][cumboatcatch_id];
	//					for(m=bm->MofY; m<12; m++){
	//						scheduledEffort = bm->EffortSchedule[buyernf][buyerns][m][expect_id];
	//						ExpCatch = bm->BlackBook[buyernf][buyerns][sp][m][expect_id];
	//						ExpEffort = bm->EffortSchedule[buyernf][buyerns][m][hist_id];
	//						expectedCatch += scheduledEffort * ExpCatch / (ExpEffort + small_num);
	//					}
	//					demand = cumCatch + expectedCatch + expectedIncreasedCatch - ownQuota - leaseQuota;
	//				}
	//				if(demand <= 0)
	//					b++;
	//			}
	//
	//			s = nsort;
	//			supply = 0;
	//			while(supply <= 0 && s >= 0){
	//				/* Deconstructing ids from single sorted array - assuming ids in
	//				the sorted array are of the form ns*tnf*(12 - month) + (12 - month)*nf + month - bm->MofY *
	//				step1 = (int)floor(ROUNDGUARD + (linearPI[sp][idd_id][s] / (nm + small_num)));	/* First step in finding fishery_id and subfleet_id *
	//				sellerns = (int)floor(ROUNDGUARD + (step1 / (tnf + small_num)));			/* Finding subfleet_id *
	//				sellernf = step1 % tnf;								/* Finding the fishery_id *
	//				sellermonth = (int)floor(ROUNDGUARD + (linearPI[sp][idd_id][s])) % nm;	/* Finding the month_id *
	//
	//				/* Expected annual quota leftover if no fishing in this month *
	//				if((bm->EffortSchedule[sellernf][sellerns][sellermonth][expect_id] > 0) && (bm->MargRent[sellernf][sellerns][sp][sellermonth] < 0)){
	//					scheduledEffort = bm->EffortSchedule[sellernf][sellerns][sellermonth][expect_id];
	//					ExpCatch = bm->BlackBook[sellernf][sellerns][sp][sellermonth][expect_id];
	//					ExpEffort = bm->EffortSchedule[sellernf][sellerns][sellermonth][hist_id];
	//					expectedDecreaseCatch = scheduledEffort * ExpCatch / (ExpEffort + small_num);
	//					prop_olease = bm->QuotaAlloc[sellernf][sellerns][sp][permolease_id];
	//					prop_tlease = bm->QuotaAlloc[sellernf][sellerns][sp][templease_id];
	//					ownQuota = (1.0 - prop_olease - prop_tlease) * bm->QuotaAlloc[sellernf][sellerns][sp][owned_id];
	//					leaseQuota = bm->QuotaAlloc[sellernf][sellerns][sp][leased_id];
	//					cumCatch = bm->QuotaAlloc[sellernf][sellerns][sp][cumboatcatch_id];
	//					/* Need supply to cover the rest of the year *
	//					expectedCatch = 0;
	//					for(m = bm->MofY; m<12; m++){
	//						scheduledEffort = bm->EffortSchedule[sellernf][sellerns][m][expect_id];
	//						ExpCatch = bm->BlackBook[sellernf][sellerns][sp][m][expect_id];
	//						ExpEffort = bm->EffortSchedule[sellernf][sellerns][m][hist_id];
	//						expectedCatch += scheduledEffort * ExpCatch / (ExpEffort + small_num);
	//					}
	//					remainderQ = ownQuota + leaseQuota - cumCatch - (expectedCatch - expectedDecreaseCatch);
	//					supply = min(expectedDecreaseCatch, remainderQ);
	//				}
	//				if(supply <= 0)
	//					s--;
	//			}
	//
	//			if((demand <= bm->recon_buffer) || (supply <= bm->recon_buffer))
	//				break;
	//
	//			down_time = bm->SUBFLEET_ECONprms[buyernf][buyerns][down_time_id];
	//			max_month_effort = bm->month_scalar * (1 - down_time) * bm->SUBFLEET_ECONprms[buyernf][buyerns][nboat_id];
	//			trade(bm, sp, bm->MofY, buyernf, buyerns, sellernf, sellerns, demand, supply, buyermonth, sellermonth, max_month_effort, &remain_demand, &remain_supply, llogfp);
	//
	//			/* Need to recalculate expected catch and planned effort as just changed schedule in
	//			quota trading routine trade() *
	//			bm->SUBFLEET_ECONprms[sellernf][sellerns][totPlanEffort_id] = 0;
	//			bm->SUBFLEET_ECONprms[sellernf][sellerns][totPlanCatch_id] = 0;
	//			for(m=bm->MofY; m<12; m++){
	//				scheduledEffort = bm->EffortSchedule[sellernf][sellerns][m][expect_id];
	//				ExpCatch = bm->BlackBook[sellernf][sellerns][sp][m][expect_id];
	//				ExpEffort = bm->EffortSchedule[sellernf][sellerns][m][hist_id];
	//				bm->SUBFLEET_ECONprms[sellernf][sellerns][totPlanEffort_id] += scheduledEffort;
	//				bm->SUBFLEET_ECONprms[sellernf][sellerns][totPlanCatch_id] += scheduledEffort * ExpCatch / ExpEffort;
	//			}
	//
	//			bm->SUBFLEET_ECONprms[buyernf][buyerns][totPlanEffort_id] = 0;
	//			bm->SUBFLEET_ECONprms[buyernf][buyerns][totPlanCatch_id] = 0;
	//			for(m=bm->MofY; m<12; m++){
	//				scheduledEffort = bm->EffortSchedule[buyernf][buyerns][m][expect_id];
	//				ExpCatch = bm->BlackBook[buyernf][buyerns][sp][m][expect_id];
	//				ExpEffort = bm->EffortSchedule[buyernf][buyerns][m][hist_id];
	//				bm->SUBFLEET_ECONprms[buyernf][buyerns][totPlanEffort_id] += scheduledEffort;
	//				bm->SUBFLEET_ECONprms[buyernf][buyerns][totPlanCatch_id] += scheduledEffort * ExpCatch / ExpEffort;
	//			}
	//
	//			/* Get new aggregate supply and demand *
	//			aggDemand = 0;
	//			aggSupply = 0;
	//			for(nf=0; nf<bm->K_num_fisheries; nf++){
	//				for(ns=0; ns<bm->FISHERYprms[nf][nsubfleets_id]; ns++){
	//					/* If no boats in the subfleet currently skip ahead *
	//					if(!bm->SUBFLEET_ECONprms[nf][ns][nboat_id]){
	//						continue;
	//					}
	//					cumCatch = bm->QuotaAlloc[nf][ns][sp][cumboatcatch_id];
	//					expectedCatch = bm->MonthAlloc[nf][ns][sp][bm->MofY];
	//					prop_olease = bm->QuotaAlloc[nf][ns][sp][permolease_id];
	//					prop_tlease = bm->QuotaAlloc[nf][ns][sp][templease_id];
	//					ownQuota = (1.0 - prop_olease - prop_tlease) * bm->QuotaAlloc[nf][ns][sp][owned_id];
	//					leaseQuota = bm->QuotaAlloc[nf][ns][sp][leased_id];
	//
	//					/** Calculate aggregate demand and supply **/
	//					/* Aggregate demand - is sum of the differences between the cumulative catch
	//					+ expected future catch and the quota held already *
	//					aggDemand += cumCatch + expectedCatch - ownQuota - leaseQuota;
	//
	//					/* Aggregate supply - here purchase quota is sold quota so need is negative *
	//					aggSupply += ownQuota + leaseQuota - cumCatch;
	//				}
	//			}
	//		}
	//
	//		/* Step 2 see if there is anyone with full effort scedule but is expected to not fill their quota *
	//
	//		s = nsort;
	//		b = 0;
	//
	//		while(s >= 0 && (b < nsort && b >= 0)){
	//			s = nsort;
	//			supply = 0;
	//			while((supply <= 0) && (s >= 0)){
	//				/* Deconstructing ids from single sorted array - assuming ids in
	//				the sorted array are of the form ns*tnf*(12 - month) + (12 - month)*nf + month - bm->MofY *
	//				step1 = (int)floor(ROUNDGUARD + (linearPI[sp][idd_id][s] / (nm + small_num)));	/* First step in finding fishery_id and subfleet_id *
	//				sellerns = (int)floor(ROUNDGUARD + (step1 / (tnf + small_num)));			/* Finding subfleet_id *
	//				sellernf = step1 % tnf;								/* Finding the fishery_id *
	//
	//				/* Use totPlanEffort and totPlanCatch calculated elsewhere *
	//				supply = 0;
	//				totPlannedEffort = bm->SUBFLEET_ECONprms[sellernf][sellerns][totPlanEffort_id];
	//				totPlannedCatch = bm->SUBFLEET_ECONprms[sellernf][sellerns][totPlanCatch_id];
	//				cumCatch = bm->QuotaAlloc[sellernf][sellerns][sp][cumboatcatch_id];
	//				expectedCatch = bm->MonthAlloc[sellernf][sellerns][sp][bm->MofY];
	//				prop_olease = bm->QuotaAlloc[sellernf][sellerns][sp][permolease_id];
	//				prop_tlease = bm->QuotaAlloc[sellernf][sellerns][sp][templease_id];
	//				ownQuota = (1.0 - prop_olease - prop_tlease) * bm->QuotaAlloc[sellernf][sellerns][sp][owned_id];
	//				leaseQuota = bm->QuotaAlloc[sellernf][sellerns][sp][leased_id];
	//				effThresh = 30.0; /* Effort threshold as max can actually do
	//									- if effort = days fishing then 30 is thresh,
	//									if seconds then 86400*30 etc *
	//				if((totPlannedEffort >= (12.0 - bm->MofY) * effThresh) && ((totPlannedCatch + cumCatch) < (ownQuota + leaseQuota))){
	//					supply = (ownQuota + leaseQuota) - (totPlannedCatch + cumCatch);
	//				}
	//
	//				if(supply <= 0)
	//					s--;
	//			}
	//
	//			b = 0;
	//			demand = 0;
	//			while((demand <= 0) && (b < nsort && b >= 0)){
	//				demand = 0;
	//				/* Deconstructing ids from single sorted array - assuming ids in
	//				the sorted array are of the form ns*tnf*(12 - month) + (12 - month)*nf + month - bm->MofY *
	//				step1 = (int)floor(ROUNDGUARD + (linearPI[sp][idd_id][b] / (nm + small_num)));	/* First step in finding fishery_id and subfleet_id *
	//				buyerns = (int)floor(ROUNDGUARD + (step1 / (tnf + small_num)));				/* Finding subfleet_id *
	//				buyernf = step1 % tnf;								/* Finding the fishery_id *
	//				buyermonth = (int)floor(ROUNDGUARD + (linearPI[sp][idd_id][b])) % nm;	/* Finding the month_id *
	//
	//				/* The expected increase in catch by filling the current month *
	//				if((max_month_effort - bm->EffortSchedule[buyernf][buyerns][buyermonth][expect_id] >= 1) && (bm->MargRent[buyernf][buyerns][sp][buyermonth] >= 0)){
	//					scheduledEffort = bm->EffortSchedule[buyernf][buyerns][buyermonth][expect_id];
	//					ExpCatch = bm->BlackBook[buyernf][buyerns][sp][buyermonth][expect_id];
	//					ExpEffort = bm->EffortSchedule[buyernf][buyerns][buyermonth][hist_id];
	//					prop_olease = bm->QuotaAlloc[buyernf][buyerns][sp][permolease_id];
	//					prop_tlease = bm->QuotaAlloc[buyernf][buyerns][sp][templease_id];
	//					ownQuota = (1.0 - prop_olease - prop_tlease) * bm->QuotaAlloc[buyernf][buyerns][sp][owned_id];
	//					leaseQuota = bm->QuotaAlloc[buyernf][buyerns][sp][leased_id];
	//					cumCatch = bm->QuotaAlloc[buyernf][buyerns][sp][cumboatcatch_id];
	//					expectedIncreasedCatch = (max_month_effort - scheduledEffort) * ExpCatch / (ExpEffort + small_num);
	//					expectedCatch = 0;
	//					for(m=bm->MofY; m<12; m++){
	//						scheduledEffort = bm->EffortSchedule[buyernf][buyerns][m][expect_id];
	//						ExpCatch = bm->BlackBook[buyernf][buyerns][sp][m][expect_id];
	//						ExpEffort = bm->EffortSchedule[buyernf][buyerns][m][hist_id];
	//			            expectedCatch += scheduledEffort * ExpCatch / (ExpEffort + small_num);
	//					}
	//					demand = cumCatch + expectedCatch + expectedIncreasedCatch - ownQuota - leaseQuota;
	//
	//				}
	//				if(demand <= 0)
	//					b++;
	//			}
	//
	//			if((demand <= bm->recon_buffer) || (supply <= bm->recon_buffer))
	//				break;
	//
	//			down_time = bm->SUBFLEET_ECONprms[buyernf][buyerns][down_time_id];
	//			max_month_effort = bm->month_scalar * (1 - down_time) * bm->SUBFLEET_ECONprms[buyernf][buyerns][nboat_id];
	//			trade(bm, sp, bm->MofY, buyernf, buyerns, sellernf, sellerns, demand, supply, buyermonth, sellermonth, max_month_effort, &remain_demand, &remain_supply, llogfp);
	//		}
	//	}
	//
	//	*/
	return;
}
