/*
 * atAdditionalTracer.h
 *
 *  Created on: Aug 30, 2012
 *      Author: bec
 */

#ifndef ATADDITIONALTRACER_H_
#define ATADDITIONALTRACER_H_

double Calculate_Phosphorus_Uptake(MSEBoxModel *bm, int guild, double grow, double PRatio);
double Calculate_Carbon_Uptake(MSEBoxModel *bm, int guild, double grow, double CRatio);
int Calculate_Phosphorus(MSEBoxModel *bm, BoxLayerValues *boxLayerInfo, HABITAT_TYPES habitatType);
int Calculate_Carbon(MSEBoxModel *bm, BoxLayerValues *boxLayerInfo, HABITAT_TYPES habitatType);

void Integrate_Ratio_Variables(MSEBoxModel *bm, double *tracerArray, double *fluxArray, double tsz, double dtsz, HABITAT_TYPES habitat);
void Calculate_Element_Release(MSEBoxModel *bm, BoxLayerValues *boxLayerInfo, int guild, double brokenDown, HABITAT_TYPES habitat, HABITAT_TYPES globalHabitat);
void Calculate_Element_Flux(MSEBoxModel *bm, double *tracerArray, double *fluxArray, HABITAT_TYPES habitat, HABITAT_TYPES callHabitat);
void PP_uptake(MSEBoxModel *bm, int sp, double uptakeN, double uptakeP, double uptakeC);

void Transfer_Element(MSEBoxModel *bm, BoxLayerValues *boxLayerInfo, HABITAT_TYPES habitatType, int toGuild, int toCohort, int fromGuild, int fromCohort, double amount, HABITAT_TYPES globalHabitatType, int isGlobal);
void Loose_Element(MSEBoxModel *bm, BoxLayerValues *boxLayerInfo, HABITAT_TYPES habitat, int guild, int cohort, long double amount, HABITAT_TYPES globalHabitat, int isGlobal);
void Loose_Element_From_Prey(MSEBoxModel *bm, BoxLayerValues *boxLayerInfo, HABITAT_TYPES habitat, int guild, int cohort, long double amount, HABITAT_TYPES globalHabitat, int isGlobal);
void Gain_Element(MSEBoxModel *bm, BoxLayerValues *boxLayerInfo, HABITAT_TYPES habitat,  int toGuild, int toCohort, int fromGuild, int fromCohort, double amount, HABITAT_TYPES globalHabitat, int isGlobal);
void Transfer_To_Pred(MSEBoxModel *bm, BoxLayerValues *boxLayerInfo, double ***spGRAZEinfo, int currentHabitat, int guild, int cohort, double scalar, double E_SP, double EDR_SP, double EDL_SP, HABITAT_TYPES globalHabitatType);
void Vertebrates_Transfer_To_Pred(MSEBoxModel *bm, BoxLayerValues *boxLayerInfo, double ***spGRAZEinfo, int guild, int cohort, double growth_scalar, double E_SP, double E_Plant, double SPgrazeDR, double E3_sp, double SPgrazeDL, double E4_sp, HABITAT_TYPES globalHabitatType, double density);

void Gain_Fish_Element(MSEBoxModel *bm, BoxLayerValues *boxLayerInfo, HABITAT_TYPES habitat, int fromGuild, int fromCohort, int toGuild, int toCohort, double amount, HABITAT_TYPES globalHabitat, int isGlobal);

#endif /* ATADDITIONALTRACER_H_ */
