/**
 * \ingroup atManageLib
 * \file atCloseKin.c
 * \brief Implementation of CK simulator
 *     \author    Beth Fulton 1/3/2022
 
 * Original close kin Author: Rich Little and Robin Thompson (for RatpackMSE)
 * Expansion and Modifications: Beth Fulton
 *
 *    <b>  Revisions</b>
 *
 *  This set of routines handles passing information to and from R harvest control rulea
 */

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <sjwlib.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <errno.h>
#include <unistd.h>
#include <sjwlib.h>
#include <atlantisboxmodel.h>
#include <atUtilLib.h>

/** Prototypes ***/
FILE * initCloseKinCatchFile(MSEBoxModel *bm);
FILE * initCloseKinSampFile(MSEBoxModel *bm);
FILE * initCloseKinPOFile(MSEBoxModel *bm);
FILE * initCloseKinHSFile(MSEBoxModel *bm);

//void writeCloseKinCatchFile(MSEBoxModel *bm, int year, int sp, FILE *fid);
void writeCloseKinSampFile(MSEBoxModel *bm, int year, int sp, FILE *fid);
void writeCloseKinNcompsfiles(MSEBoxModel *bm, int year, int sp, FILE *fidPO, FILE *fidHS, double *****sim_nkin_HS, double *****sim_nkin_PO, double ****ncomps_HS_yaya, double *****ncomps_PO_syaya);

int poissonRandom(double lambda);

/******************************************************************************

CKsimulator -  main entry for the CK simulator

******************************************************************************/
void CKsimulator(MSEBoxModel *bm, int sp, int year){

    int y = 0, a  = 0, j = 0;
    //int r = 0;
    int s = 0;
    int l;
    
    //double temp = 0;
    //double *tmpVec;

    int last_dy = bm->CloseKinEst[sp].last_sy;
    int nsampy = bm->CloseKinEst[sp].last_sy - bm->CloseKinEst[sp].first_sy + 1;
    //int Asamp_propto = 0; //0 for catch
    //int age_power = 1;
    double sumx = 0, sumsx = 0;
    double log_fecIN_sa = 0;
    //double tempC = 0;
    double cumul_psurv = 0;
    int curr_a = 0, y2 = 0;
    double stocktot = 0;
    int a2 = 0;
    int b1 = 0, b2 = 0;
    int a1_b2 = 0;
    double demog_Pr_HSP = 0, Pr_par1_was_ap_b1 = 0, Pr_par2_is_par1_if_alive = 0;
    int ap_b1;
    int nLinesCompPO = 0, nLinesCompHS = 0;
    double comps = 0, nsamps1 = 0, nsamps2 = 0;
    int nLinesKinPO = 0, nLinesKinHS = 0;
    int this_min = (int)(bm->RBCestimation.RBCspeciesParam[sp][HistYrMin_id]);
    int this_maxage = (int)(bm->RBCestimation.RBCspeciesParam[sp][MaxAge_id]);
        
    /* Leave here for now */
    //TODO: CHECK ORDER OF DIMENSIONS IN THESE DECLARARIONS VS THEIR USE
    double *nsamps_y = Util_Alloc_Init_1D_Double(nsampy, 0.0);
    double ***nsamps_sya = Util_Alloc_Init_3D_Double(bm->RBCestimation.RBCspeciesParam[sp][Nsexes_id], nsampy, bm->RBCestimation.RBCspeciesParam[sp][MaxAge_id], 0.0);

    double *x = Util_Alloc_Init_1D_Double(bm->RBCestimation.RBCspeciesParam[sp][MaxAge_id], 0.0);
    double *sx = Util_Alloc_Init_1D_Double(bm->RBCestimation.RBCspeciesParam[sp][Nsexes_id], 0.0);
    double **inv_totfec_sy = Util_Alloc_Init_2D_Double(bm->RBCestimation.RBCspeciesParam[sp][Nsexes_id], last_dy, 0.0);

    double ***n_sya = Util_Alloc_Init_3D_Double(bm->RBCestimation.RBCspeciesParam[sp][Nsexes_id], last_dy, bm->RBCestimation.RBCspeciesParam[sp][MaxAge_id], 0.0);
    double ****psurv_syay = Util_Alloc_Init_4D_Double(bm->RBCestimation.RBCspeciesParam[sp][Nsexes_id], last_dy, bm->RBCestimation.RBCspeciesParam[sp][MaxAge_id], last_dy, 0.0);
    
    double *****Pr_PO_syaya = Util_Alloc_Init_5D_Double(bm->RBCestimation.RBCspeciesParam[sp][Nsexes_id], bm->CloseKinEst[sp].last_sy, bm->RBCestimation.RBCspeciesParam[sp][MaxAge_id], bm->CloseKinEst[sp].last_sy, bm->RBCestimation.RBCspeciesParam[sp][MaxAge_id], 0.0);
    double *****Pr_HS_syaya = Util_Alloc_Init_5D_Double(bm->RBCestimation.RBCspeciesParam[sp][Nsexes_id], bm->CloseKinEst[sp].last_sy, bm->RBCestimation.RBCspeciesParam[sp][MaxAge_id], bm->CloseKinEst[sp].last_sy, bm->RBCestimation.RBCspeciesParam[sp][MaxAge_id], 0.0);

    double ****Pr_GG_yaya = Util_Alloc_Init_4D_Double(bm->CloseKinEst[sp].last_sy, bm->RBCestimation.RBCspeciesParam[sp][MaxAge_id], bm->CloseKinEst[sp].last_sy, bm->RBCestimation.RBCspeciesParam[sp][MaxAge_id], 0.0);
    double ****Pr_HS_yaya = Util_Alloc_Init_4D_Double(bm->CloseKinEst[sp].last_sy, bm->RBCestimation.RBCspeciesParam[sp][MaxAge_id], bm->CloseKinEst[sp].last_sy, bm->RBCestimation.RBCspeciesParam[sp][MaxAge_id], 0.0);

    double *****ncomps_PO_syaya = Util_Alloc_Init_5D_Double(bm->RBCestimation.RBCspeciesParam[sp][Nsexes_id], bm->CloseKinEst[sp].last_sy, bm->RBCestimation.RBCspeciesParam[sp][MaxAge_id], bm->CloseKinEst[sp].last_sy, bm->RBCestimation.RBCspeciesParam[sp][MaxAge_id], 0.0);
    double ****ncomps_HS_yaya = Util_Alloc_Init_4D_Double(bm->CloseKinEst[sp].last_sy, bm->RBCestimation.RBCspeciesParam[sp][MaxAge_id], bm->CloseKinEst[sp].last_sy, bm->RBCestimation.RBCspeciesParam[sp][MaxAge_id], 0.0);
    
    double *****Ekin_HS_syaya = Util_Alloc_Init_5D_Double(bm->RBCestimation.RBCspeciesParam[sp][Nsexes_id], bm->CloseKinEst[sp].last_sy, bm->RBCestimation.RBCspeciesParam[sp][MaxAge_id], bm->CloseKinEst[sp].last_sy, bm->RBCestimation.RBCspeciesParam[sp][MaxAge_id], 0.0);
    double *****Ekin_PO_syaya = Util_Alloc_Init_5D_Double(bm->RBCestimation.RBCspeciesParam[sp][Nsexes_id], bm->CloseKinEst[sp].last_sy, bm->RBCestimation.RBCspeciesParam[sp][MaxAge_id], bm->CloseKinEst[sp].last_sy, bm->RBCestimation.RBCspeciesParam[sp][MaxAge_id], 0.0);
    
    double *****sim_nkin_PO = Util_Alloc_Init_5D_Double(bm->RBCestimation.RBCspeciesParam[sp][Nsexes_id], bm->CloseKinEst[sp].last_sy, bm->RBCestimation.RBCspeciesParam[sp][MaxAge_id], bm->CloseKinEst[sp].last_sy, bm->RBCestimation.RBCspeciesParam[sp][MaxAge_id], 0.0);
    double *****sim_nkin_HS = Util_Alloc_Init_5D_Double(bm->RBCestimation.RBCspeciesParam[sp][Nsexes_id], bm->CloseKinEst[sp].last_sy, bm->RBCestimation.RBCspeciesParam[sp][MaxAge_id], bm->CloseKinEst[sp].last_sy, bm->RBCestimation.RBCspeciesParam[sp][MaxAge_id], 0.0);

    
    Util_Init_1D_Double(bm->CloseKinEst[sp].fec_expo_s, bm->RBCestimation.RBCspeciesParam[sp][Nsexes_id], 1.0);
    Util_Init_2D_Double(bm->CloseKinEst[sp].fec_sa, bm->RBCestimation.RBCspeciesParam[sp][Nsexes_id], bm->RBCestimation.RBCspeciesParam[sp][MaxAge_id],  0.0);
    Util_Init_1D_Double(bm->CloseKinEst[sp].mature, bm->RBCestimation.RBCspeciesParam[sp][Nlen_id], 0.0);
    Util_Init_1D_Double(bm->CloseKinEst[sp].lengths, bm->RBCestimation.RBCspeciesParam[sp][Nlen_id], 0.0);

    /**
    for (y = bm->CloseKinEst[sp].first_sy; y < bm->CloseKinEst[sp].last_sy; y++) {
        nsamps_y[y] = int(bm->CloseKinEst.nTotSamps / nsampy);
        fprintf(bm->logFile, "Time: %e %s CKsimulator: nsamps_y[%d] : %d\n", bm->dayt, FunctGroupArray[sp].groupCode, y, nsamps_y[y]);
    }
     **/
    
    if (year > bm->RBCestimation.RBCspeciesParam[sp][HistYrMax_id])
        nsampy = bm->CloseKinEst[sp].nCKsamples_proj; //specficiation per year
    else {
        nsampy = bm->CloseKinEst[sp].nCKsamples_hist / (bm->CloseKinEst[sp].last_sy - bm->CloseKinEst[sp].first_sy + 1); //specficiation per year
    }

    for (y = bm->CloseKinEst[sp].first_dy; y < last_dy; y++) {
        nsamps_y[y] = nsampy;
        fprintf(bm->logFile, "Time: %e %s y %d CKsimulator: nsamps_y %e\n", bm->dayt, FunctGroupArray[sp].groupCode, y, nsamps_y[y]);
    }

    for (l = 0 ; l < bm->RBCestimation.RBCspeciesParam[sp][Nlen_id]; l++) {
        bm->CloseKinEst[sp].lengths[l] = bm->RBCestimation.RBCspeciesArray[sp].LoLenBin[l] + (bm->RBCestimation.RBCspeciesArray[sp].HiLenBin[l] - bm->RBCestimation.RBCspeciesArray[sp].LoLenBin[l]) / 2.0;
    }

    for (j = 0; j < FunctGroupArray[sp].numStocks; j++){
    // Maturity as a function of length
        for (l = 0 ; l < bm->RBCestimation.RBCspeciesParam[sp][Nlen_id]; l++) {
            bm->CloseKinEst[sp].mature[l] = 1.0 / (1.0 + exp(bm->RBCestimation.RBCspeciesArray[sp].SlopeMat[j] * (bm->CloseKinEst[sp].lengths[l] - bm->RBCestimation.RBCspeciesArray[sp].Mat50[j])));
        }

        for (s = 0; s < bm->RBCestimation.RBCspeciesParam[sp][Nsexes_id]; s++){
            for (a = 0; a < bm->RBCestimation.RBCspeciesParam[sp][MaxAge_id]; a++) {
                bm->CloseKinEst[sp].fec_sa[s][a] = 0;
                for (int l = 0; l < bm->RBCestimation.RBCspeciesParam[sp][Nlen_id]; l++) {
                    bm->CloseKinEst[sp].fec_sa[s][a] += bm->RBCestimation.RBCspeciesArray[sp].FracLenS[j][s][a][this_min][l] * bm->CloseKinEst[sp].mature[l];
                }
            }
        }
    }

    int amat = 0;
    for (a = 0; a < bm->RBCestimation.RBCspeciesParam[sp][MaxAge_id]; a++) {
        if (bm->CloseKinEst[sp].fec_sa[FEMALE_ID][a] < bm->RBCestimation.RBCspeciesParam[sp][thresh_mat_id]) {
            amat = a;
        }
    }

    // Writeout closekin catch
    if(!bm->closekinCFile) {
        bm->closekinCFile = initCloseKinCatchFile(bm);
    }
    if(!bm->closekinSampFile) {
        bm->closekinSampFile = initCloseKinSampFile(bm);
    }
    //writeCloseKinCatchFile(bm, year, sp, bm->closekinCFile);
    writeCloseKinSampFile(bm, year, sp, bm->closekinSampFile);
    Util_Close_Output_File(bm->closekinCFile);
    Util_Close_Output_File(bm->closekinSampFile);
    
    for (y = bm->CloseKinEst[sp].first_sy; y < bm->CloseKinEst[sp].last_sy; y++){
        sumsx = 0;
        for (s = 0; s < bm->RBCestimation.RBCspeciesParam[sp][Nsexes_id]; s++){
            sumx = 0;
            sx[s] = 0;
            for (a = 0; a < bm->RBCestimation.RBCspeciesParam[sp][MaxAge_id]; a++){
                x[a] = pow( bm->CloseKinEst[sp].samp_prop_to[s][y][a], bm->CloseKinEst[sp].samp_power);
                sumx += x[a];
                sx[s] += pow(bm->CloseKinEst[sp].samp_prop_to[s][y][a], bm->CloseKinEst[sp].samp_power);
                sumsx += pow(bm->CloseKinEst[sp].samp_prop_to[s][y][a], bm->CloseKinEst[sp].samp_power);
            }
            for (a = 0; a < bm->RBCestimation.RBCspeciesParam[sp][MaxAge_id]; a++){
                nsamps_sya[s][y][a] = nsamps_y[y] * (x[a] / sumx);
            }
        }

        // If nsamps_y[y] is across both sexes, scale by the sex ratio
        for (s = 0; s < bm->RBCestimation.RBCspeciesParam[sp][Nsexes_id]; s++){
            fprintf(bm->logFile, "Time: %e %s sex ratio %e\n", bm->dayt, FunctGroupArray[sp].groupCode, sx[s] / sumsx);
            for (a = 0; a < bm->RBCestimation.RBCspeciesParam[sp][MaxAge_id]; a++){
                nsamps_sya[s][y][a] *= sx[s] / sumsx;
            }
        }

        /*
        for (s = 0; s < bm->RBCestimation.RBCspeciesParam[sp][Nsexes_id]; s++){
            fprintf(bm->logFile, "Time: %e %s nsamps_sya ", bm->dayt, s);
            for (a = 0; a < bm->RBCestimation.RBCspeciesParam[sp][MaxAge_id]; a++){
                fprintf(bm->logFile, "%d ", nsamps_sya[s][y][a]);
            }
            fprintf(bm->logFile, "\n");
        }
        */
    }

 
    // Pr_stuff
    //fprintf(bm->logFile, "Time: %e %s CKsimulator: Fecundity ", bm->dayt, FunctGroupArray[sp].groupCode);
    for (s = 0; s < bm->RBCestimation.RBCspeciesParam[sp][Nsexes_id]; s++){
        //fprintf(bm->logFile, "sex ", s);
        for (a = 0;a < bm->RBCestimation.RBCspeciesParam[sp][MaxAge_id]; a++){
            log_fecIN_sa = log(bm->CloseKinEst[sp].fec_sa[s][a] + 0.00001);
            bm->CloseKinEst[sp].fec_sa[s][a] = exp(bm->CloseKinEst[sp].fec_expo_s[s] * log_fecIN_sa);

            //fprintf(bm->logFile, "fec_sa-%d %e ", a, bm->CloseKinEst[sp].fec_sa[s][a]);
        }
        //fprintf(bm->logFile, "\n");
    }

    //fprintf(bm->logFile, "Time: %e %s CKsimulator: ", bm->dayt, FunctGroupArray[sp].groupCode);
    for (y = bm->CloseKinEst[sp].first_dy; y < last_dy; y++) {
        for (s = 0; s < bm->RBCestimation.RBCspeciesParam[sp][Nsexes_id]; s++){
            for (j = 0; j < FunctGroupArray[sp].numStocks; j++){
                for (a = 0; a < bm->RBCestimation.RBCspeciesParam[sp][MaxAge_id]; a++){
                    n_sya[s][y][a] += bm->RBCestimation.RBCspeciesArray[sp].CKSurveyNum[j][s][a][y];
                }
            }
            /*
            fprintf(bm->logFile, "%d ", s);
            for (a = 0; a < bm->RBCestimation.RBCspeciesParam[sp][MaxAge_id]; a++) {
                fprintf(bm->logFile, "Time: %e %s n_sya: %e ", n_sya[s][y][a]);
            }
            fprintf(bm->logFile, "\n");
            */
        }
    }

    //fprintf(bm->logFile, "Time: %e %s CKsimulator: ", bm->dayt, FunctGroupArray[sp].groupCode);
    for (s = 0; s < bm->RBCestimation.RBCspeciesParam[sp][Nsexes_id]; s++){
        for (y = bm->CloseKinEst[sp].first_dy + 1; y < last_dy; y++) {
            //fprintf(bm->logFile, "%d n-ratio: %e ", s, (n[s][y][a] / n[s][y-1][a-1]));

            for (a = 1; a < bm->RBCestimation.RBCspeciesParam[sp][MaxAge_id]; a++){
                  psurv_syay[s][y-1][a-1][y] = n_sya[s][y][a] / n_sya[s][y-1][a-1];
            }

            psurv_syay[s][y-1][this_maxage][y] = n_sya[s][y][this_maxage] / n_sya[s][y- 1][this_maxage-1];

            /*
            for (a = 1;a < bm->RBCestimation.RBCspeciesParam[sp][MaxAge_id]; a++) {
                fprintf(bm->logFile, "surv_syay[s][y-1][a-1][y] ", s, y-1, a-1, y, surv_syay[s][y-1][a-1][y]);
            }
            fprintf(bm->logFile, "\n");
            */
        }
    }

    //multi-year surv probs have already done one year ahead, so just multiply psurv_syay
    //fprintf(bm->logFile, "Time: %e %s CKsimulator: ", bm->dayt, FunctGroupArray[sp].groupCode);
    for (s = 0; s < bm->RBCestimation.RBCspeciesParam[sp][Nsexes_id]; s++){
        //fprintf(bm->logFile, "s %d  ", s);
        for (y = bm->CloseKinEst[sp].first_dy; y < last_dy - 1; y++) {
            /*
            fprintf(bm->logFile, "y %d ", y);
            for (y2 = y; y2 < last_dy; y2++) {
                fprintf(bm->logFile, "y2 %d", y2);
            }
            fprintf(bm->logFile, "\n");
            */
            
            for (a = 0; a < bm->RBCestimation.RBCspeciesParam[sp][MaxAge_id]; a++){
                psurv_syay[s][y][a][y] = 1; //diagonal
                cumul_psurv = 1;
                curr_a = a;
                for (y2 = y+1; y2 < last_dy; y2++) {
                    cumul_psurv *= psurv_syay[s][y2-1][curr_a][y2];
                    psurv_syay[s][y][curr_a][y2] = cumul_psurv;
                    if (curr_a < bm->RBCestimation.RBCspeciesParam[sp][MaxAge_id]) {
                        curr_a++; // ROBIN what does curr_a do?
                    }
                }

                /*
                fprintf(bm->logFile, "age ", a);
                for (y2 = y; y2 < last_dy; y2++) {
                    fprintf(bm->logFile, "psurv_syay[s][y][a][y2]: %e ", s, y, a, y2, psurv_syay[s][y][a][y2]);
                }
                fprintf(bm->logFile, "\n");
                */
            }
        }
    }

    for (s = 0; s < bm->RBCestimation.RBCspeciesParam[sp][Nsexes_id]; s++){
        for (y = bm->CloseKinEst[sp].first_dy + 1; y < last_dy; y++){
            for (a = 0; a < bm->RBCestimation.RBCspeciesParam[sp][MaxAge_id]; a++){
                stocktot += n_sya[s][y][a] * bm->CloseKinEst[sp].fec_sa[s][a];
            }
            if (stocktot > 0) {
                inv_totfec_sy[s][y] = 1 / stocktot;
            } else {
                fprintf(bm->logFile, "Time: %e %s  CKsim stocktot <= 0\n", bm->dayt, FunctGroupArray[sp].groupCode);
                //return;  // If this return is real remember to empty out the arrays
            }
        }
    }

    /*** PO **/
    for (s = 0; s < bm->RBCestimation.RBCspeciesParam[sp][Nsexes_id]; s++){
        for (y = bm->CloseKinEst[sp].first_sy; y < bm->CloseKinEst[sp].last_sy; y++){
            for (a = amat; a < bm->RBCestimation.RBCspeciesParam[sp][MaxAge_id]; a++){
                b1 = y - a;
                for (y2 = bm->CloseKinEst[sp].first_sy; y2 < bm->CloseKinEst[sp].last_sy; y2++) {
                    for (a2 = amat; a2 < bm->RBCestimation.RBCspeciesParam[sp][MaxAge_id]; a2++){
                        b2 = y2 - a2;
                        if (b2 >= bm->CloseKinEst[sp].first_dy){
                            a1_b2 = a - (y - b2);
                            if ((a1_b2 >= amat) && (y > b2)){
                                Pr_PO_syaya[s][y][a][y2][a2] = bm->CloseKinEst[sp].fec_sa[s][a1_b2] * inv_totfec_sy[s][b2];
                            } // alive and bm->CloseKinEst[sp].mature
                        } //if not born before pop
                    }
                }
            }
        }
    }

    /*** HS **/
    for (s = 0; s < bm->RBCestimation.RBCspeciesParam[sp][Nsexes_id]; s++){
        for (y = bm->CloseKinEst[sp].first_sy; y < bm->CloseKinEst[sp].last_sy; y++){
            for (a = amat; a < bm->RBCestimation.RBCspeciesParam[sp][MaxAge_id]; a++){
                b1 = y - a;
                if ((b1 >= bm->CloseKinEst[sp].first_dy) && (b1 < last_dy)) {
                    for (y2 = bm->CloseKinEst[sp].first_sy; y2 < bm->CloseKinEst[sp].last_sy; y2++) {
                        for (a2 = amat;a2 < bm->RBCestimation.RBCspeciesParam[sp][MaxAge_id]; a2++){
                            b2 = y2 - a2;
                            if ((b2 >= bm->CloseKinEst[sp].first_dy) && (b2 < last_dy) )  {
                                // Note, not restricting this to b1>b2, so this is less efficient than it might be
                                demog_Pr_HSP = 0;

                                // *Putative* parental age in b1
                                // Propto numbers-at-age * fec-at-age by sex in b1
                                for (ap_b1 = amat; ap_b1 < bm->RBCestimation.RBCspeciesParam[sp][MaxAge_id]; ap_b1++) {
    
                                    Pr_par1_was_ap_b1 = n_sya[s][b1][ap_b1] * bm->CloseKinEst[sp].fec_sa[s][ap_b1] * inv_totfec_sy[s][b1];
    
                                    // Parent's age in b2, could be in the plus group
                                    int ap_b2 = 0;
                                    if ((ap_b1 + (b2 - b1)) > bm->RBCestimation.RBCspeciesParam[sp][MaxAge_id]) {
                                        ap_b2 = bm->RBCestimation.RBCspeciesParam[sp][MaxAge_id];
                                    } else {
                                        ap_b2 = ap_b1 + (b2 - b1);
                                    }
                                    if( ap_b2 >= amat)  {
                                        // note *no* n_sya term next because it's just 1 candidate animal
                                        Pr_par2_is_par1_if_alive = bm->CloseKinEst[sp].fec_sa[s][ap_b2] * inv_totfec_sy[s][b2];

                                        // Cumul survival prob of parent from b1 to b2, given age ap_b1 in b1
                                        demog_Pr_HSP +=    Pr_par1_was_ap_b1 * Pr_par2_is_par1_if_alive * psurv_syay[s][b1][ap_b1][b2];  // if b2>b1 this will be zero instead of 1; but I think those cases are overwritten anyway

                                    } // for ap_b1
    
                                    Pr_HS_syaya[s][y][a][y2][a2] = demog_Pr_HSP; // *bm->CloseKinEst.bm->CloseKinEst.HSP_false_neg_loss_rate; // moved to Likelihood calcs
                                    Pr_HS_yaya[y][a][y2][a2] += Pr_HS_syaya[s][y][a][y2][a2]; //ROBIN: Pr_HS_yaya[y][a][y2][a2] isn;t really used anywhere. It's a bit complicated in the original. Is it needed?
    
                                } // if parent bm->CloseKinEst[sp].mature in b2

                            } //  if b2 within temporal scope of model
                        }
                    }
                }
            }
        }
    }

    /*************** R - calc_ncomps ***************/
    
    /*** PO **/
    for (s = 0; s < bm->RBCestimation.RBCspeciesParam[sp][Nsexes_id]; s++) {
        for (y = bm->CloseKinEst[sp].first_sy; y < bm->CloseKinEst[sp].last_sy; y++) {
            for (a = amat; a < bm->RBCestimation.RBCspeciesParam[sp][MaxAge_id]; a++) {
                b1 = y - a;
                nsamps1 = nsamps_sya[s][y][a];
                for (y2 = bm->CloseKinEst[sp].first_sy; y2 < bm->CloseKinEst[sp].last_sy; y2++) { //  I could start in year y1, but I've chosen to populate this whole matrix instead
                    for (a2 = 0; a2 < bm->RBCestimation.RBCspeciesParam[sp][MaxAge_id]; a2++) {
                        b2 = y2 - a2;
                        comps = 0;
                        for (int s1 = 1; s1 < bm->RBCestimation.RBCspeciesParam[sp][Nsexes_id]; s1++)
                            comps += nsamps1 * nsamps_sya[s1][y2][a2];  // sex of offspring doesn't matter

                        // PO
                        if ((y > b2) && (b2 >= bm->CloseKinEst[sp].first_dy)) { // 'parent' must have been alive in b2 (assuming lethal sampling); b2 must fall within the demographic span of the model
                            ncomps_PO_syaya[s][y][a][y2][a2] = comps;
                            if ((int) comps) {
                                nLinesCompPO++;
                            }
                        }
                    } // a2
                } // y2
            } // a1
        } // y
    } // s

    /*** HSP  **/
    // animal 1 is older, it might be the grandparent
    for (y = bm->CloseKinEst[sp].first_sy; y < bm->CloseKinEst[sp].last_sy; y++) {
        for (a = amat; a < bm->RBCestimation.RBCspeciesParam[sp][MaxAge_id]; a++) { // NB *no* plus-group comps, ALSO NOTE THAT ONLY bm->CloseKinEst[sp].mature INDIVIDUALS ARE ALLOWED TO BE POTENTIAL GRANDPARENTS, ASSUMING LETHAL SAMPLING
            b1 = y - a;
            nsamps1 = 0;
            for (s = 0; s < bm->RBCestimation.RBCspeciesParam[sp][Nsexes_id]; s++) {
                nsamps1 += nsamps_sya[s][y][a];
            }

            for (y2 = bm->CloseKinEst[sp].first_sy; y2 < bm->CloseKinEst[sp].last_sy; y2++) { //  I could start in year y1, but I've chosen to populate this whole matrix instead
                for (a2 = 0; a2 < bm->RBCestimation.RBCspeciesParam[sp][MaxAge_id]; a2++) {
                    b2 = y2 - a2;
                    // Removing this, because it was meant to compensate for within group comparisons, and those
                    // only occur when ya=y2 AND a1=a2 and those are same cohort animals, so we ignore them anyway
                    if ((y == y2) && (a == a2)) {
                        comps =  ((nsamps1 - 1.0) * nsamps1) / 2.0; // because nsamps1=nsamps2
                    }
                    else  {
                    //if (b1 < b2) { //avoid double counting by only considering cases where the first animal is older
                    // do them all, avoid double counting within the LnL calcs
                        nsamps2 = 0;
                        for (s = 0; s < bm->RBCestimation.RBCspeciesParam[sp][Nsexes_id]; s++) {
                            nsamps2 += nsamps_sya[s][y2][a2];
                        }
                        comps = nsamps1 * nsamps2;

                        if ((y == (bm->CloseKinEst[sp].last_sy - 1)) && (a == (bm->RBCestimation.RBCspeciesParam[sp][MaxAge_id] - 1)) && (y2 == (bm->CloseKinEst[sp].last_sy - 2)) && (a2 == (bm->RBCestimation.RBCspeciesParam[sp][MaxAge_id] - 2))) {
                            fprintf(bm->logFile, "Time: %e %s nsamps_sya[FEMALE][%d][%d]: %e nsamps_sya[MALE]: %e nsamps1: %e nsamps_sya[FEMALE][%d][%d]: %e nsamps_sya[MALE]: %e nsamps2: %e comps: %e\n", bm->dayt, FunctGroupArray[sp].groupCode, y, a, nsamps_sya[FEMALE_ID][y][a], nsamps_sya[MALE_ID][y][a], nsamps1, y2, a2, nsamps_sya[FEMALE_ID][y2][a2], nsamps_sya[MALE_ID][y2][a2], nsamps2, comps);
                        }

                    }
                    ncomps_HS_yaya[y][a][y2][a2] = comps;
                    if ((int) comps) nLinesCompHS++;

                } // a2
            } // y2
        } // a1
    } // y1
        
    /*** calc_exp_kin **/
    // Multiply probabilities by ncomps
    for (y = bm->CloseKinEst[sp].first_sy; y < bm->CloseKinEst[sp].last_sy; y++) {    // NB skip plus-group since its earlier age is uncertain
        for (a = 0; a < bm->RBCestimation.RBCspeciesParam[sp][MaxAge_id]; a++) {
            for (y2 = bm->CloseKinEst[sp].first_sy; y2 < bm->CloseKinEst[sp].last_sy; y2++) {
                for (a2 = 0; a2 < bm->RBCestimation.RBCspeciesParam[sp][MaxAge_id]; a2++) {
                    for (s = 0; s < bm->RBCestimation.RBCspeciesParam[sp][Nsexes_id]; s++) {
                        Ekin_HS_syaya[s][y][a][y2][a2] = Pr_HS_syaya[s][y][a][y2][a2] * ncomps_HS_yaya[y][a][y2][a2]; // s is s_parent
                        Ekin_PO_syaya[s][y][a][y2][a2] = Pr_PO_syaya[s][y][a][y2][a2] * ncomps_PO_syaya[s][y][a][y2][a2];  // s is s1, the candidate parent
                    } // s1
                } // a2
            } // y2
        } // a1
    } // y1

    //fprintf(bm->logFile, "Time: %e %s ", bm->dayt, FunctGroupArray[sp].groupCode);
    for (s = 0; s < bm->RBCestimation.RBCspeciesParam[sp][Nsexes_id]; s++) {
        for (y = bm->CloseKinEst[sp].first_sy; y < bm->CloseKinEst[sp].last_sy; y++) {    // NB skip plus-group since its earlier age is uncertain
            for (a = 0; a < bm->RBCestimation.RBCspeciesParam[sp][MaxAge_id]; a++) {
                for (y2 = bm->CloseKinEst[sp].first_sy; y2 < bm->CloseKinEst[sp].last_sy; y2++) {
                    for (a2 = 0; a2 < bm->RBCestimation.RBCspeciesParam[sp][MaxAge_id]; a2++) {
                        fprintf(bm->logFile, "Ekin_HS_syaya[%d][%d][%d][%d][%d] %e ", s, y, a, y2, a2, Ekin_HS_syaya[s][y][a][y2][a2]);
                    }
                    fprintf(bm->logFile, "\n");
                }
            }
        }
    }
    
    /*** Draw kin pairs **/
    nLinesKinPO = 0;
    nLinesKinHS = 0;
    for (y = bm->CloseKinEst[sp].first_sy; y < bm->CloseKinEst[sp].last_sy; y++) {    // NB skip plus-group since its earlier age is uncertain
        for (a = 0; a < bm->RBCestimation.RBCspeciesParam[sp][MaxAge_id]; a++) {
            for (y2 = bm->CloseKinEst[sp].first_sy; y2 < bm->CloseKinEst[sp].last_sy; y2++) {
                for (a2 = 0; a2 < bm->RBCestimation.RBCspeciesParam[sp][MaxAge_id]; a2++) {
                    for (s = 0; s < bm->RBCestimation.RBCspeciesParam[sp][Nsexes_id]; s++) {
                        sim_nkin_PO[s][y][a][y2][a2] = poissonRandom(Ekin_PO_syaya[s][y][a][y2][a2]);
                        if((int) sim_nkin_PO[s][y][a][y2][a2]) {
                            nLinesKinPO++;
                        }
                        sim_nkin_HS[s][y][a][y2][a2] = poissonRandom(Ekin_HS_syaya[s][y][a][y2][a2]);
                        if((int) sim_nkin_HS[s][y][a][y2][a2]) {
                            nLinesKinHS++;
                        }
                    }
                }
            }
        }
    }

    /** write ncomps to files */
    if(!bm->CloseKinPOFile) {
        bm->CloseKinPOFile = initCloseKinPOFile(bm);
    }
    if(!bm->CloseKinHSFile) {
        bm->CloseKinHSFile = initCloseKinHSFile(bm);
    }
    writeCloseKinNcompsfiles(bm, sp, year, bm->CloseKinPOFile,  bm->CloseKinHSFile, sim_nkin_HS, sim_nkin_PO, ncomps_HS_yaya, ncomps_PO_syaya);
    Util_Close_Output_File(bm->CloseKinPOFile);
    Util_Close_Output_File(bm->CloseKinHSFile);
    
    /* Clean up */
    free1d(nsamps_y);
    free3d(nsamps_sya);
    free1d(x);
    free1d(sx);
    free2d(inv_totfec_sy);
    free3d(n_sya);
    free4d(psurv_syay);
    free5d(Pr_PO_syaya);
    free5d(Pr_HS_syaya);
    free4d(Pr_GG_yaya);
    free4d(Pr_HS_yaya);
    free5d(ncomps_PO_syaya);
    free4d(ncomps_HS_yaya);
    free5d(Ekin_HS_syaya);
    free5d(Ekin_PO_syaya);
    free5d(sim_nkin_PO);
    free5d(sim_nkin_HS);
   
    return;
}

/*********************************************************************************************************************************************
 
 Close kin related files
 
 ******/

FILE * initCloseKinCatchFile(MSEBoxModel *bm) {
    FILE *fp;
    char *fname = "CloseKinCatch.out";

    /* Create file */
    if ((fp = Util_fopen(bm, fname, "w")) == NULL)
        quit("initCloseKinCatchFile: Can't open %s\n", fname);

    /* File content - nothing yet as added in another routine */
    fprintf(fp, "Time year species yr sex fleet age region CatchNumbersAtAge\n");

    /* Return file pointer */
    return (fp);
}

/** Actually needed? If yes then need to iniitalise the matrix and generate the content
void writeCloseKinCatchFile(MSEBoxModel *bm, int year, int sp, FILE *fid) {
    int y, s, f, a, r;

    for (y = bm->CloseKinEst[sp].first_sy; y < bm->CloseKinEst[sp].last_sy; y++){
        for (s = 0; s < bm->RBCestimation.RBCspeciesParam[sp][Nsexes_id]; s++){
            for (f = 0; f < bm->RBCestimation.RBCspeciesParam[sp][NumFisheries_id]; f++){
                for (a = 0; a < bm->RBCestimation.RBCspeciesParam[sp][MaxAge_id]; a++){
                    for ( r = 0; r < bm->RBCestimation.RBCspeciesParam[sp][NumRegions_id]; r++){
                        fprintf(fid, "%e %d %s %d %d %d %d %d %e\n", bm->dayt, year, FunctGroupArray[sp].groupCode, y, s, f, a, r, bm->RBCestimation.RBCspeciesArray[sp].CatchNumbersAtAge[y][r][f][s][a]); // check that year is correct here
                    }
                }
            }
        }
    }
    return;
}
*/

FILE * initCloseKinSampFile(MSEBoxModel *bm) {
    FILE *fp;
    char *fname = "CloseKinSamp.out";

    /* Create file */
    if ((fp = Util_fopen(bm, fname, "w")) == NULL)
        quit("initCloseKinSampFile: Can't open %s\n", fname);

    /* File content - nothing yet as added in another routine */
    fprintf(fp, "Time year species sex yr age samp_prop_to\n");

    /* Return file pointer */
    return (fp);
}

void writeCloseKinSampFile(MSEBoxModel *bm, int year, int sp, FILE *fid) {
    int y, s, a;
    //int f, r;

    for (y = bm->CloseKinEst[sp].first_sy; y < bm->CloseKinEst[sp].last_sy; y++){
        for (s = 0; s < bm->RBCestimation.RBCspeciesParam[sp][Nsexes_id]; s++){
            for (a = 0; a < bm->RBCestimation.RBCspeciesParam[sp][MaxAge_id]; a++){
                /* Redid a bit differnt to Ratpack
                for (f = 0; f < bm->RBCestimation.RBCspeciesParam[sp][NumFisheries_id]; f++){
                    bm->CloseKinEst[sp].samp_prop_to[s][y][a] += bm->RBCestimation.RBCspeciesArray[sp].CatchNumbersAtAge[y][r][f][s][a]; // check that year is corrct here
                }
                */
                bm->CloseKinEst[sp].samp_prop_to[s][y][a] += bm->RBCestimation.RBCspeciesArray[sp].ageprops[s][a] *  FunctGroupArray[sp].speciesParams[samplesize_id];
                if (a <  bm->CloseKinEst[sp].a_ck_thresh) bm->CloseKinEst[sp].samp_prop_to[s][y][a] *= bm->CloseKinEst[sp].a_ck_scalar;
                else if (a < bm->CloseKinEst[sp].a_ck_thresh2) bm->CloseKinEst[sp].samp_prop_to[s][y][a] *= bm->CloseKinEst[sp].a_ck_scalar2;

                fprintf(bm->logFile, "%e %d %s %d %d %d CloseKinEst: %e\n", bm->dayt, year, FunctGroupArray[sp].groupCode, s, y, a, bm->CloseKinEst[sp].samp_prop_to[s][y][a]);
            }
        }
    }
    return;
}

FILE * initCloseKinPOFile(MSEBoxModel *bm) {
    FILE *fp;
    char *fname = "ncomps_PO.out";

    /* Create file */
    if ((fp = Util_fopen(bm, fname, "w")) == NULL)
        quit("initCloseKinPOFile: Can't open %s\n", fname);

    /* File content - nothing yet as added in another routine */
    fprintf(fp, "Time year species \n");

    /* Return file pointer */
    return (fp);
}

FILE * initCloseKinHSFile(MSEBoxModel *bm) {
    FILE *fp;
    char *fname = "ncomps_HS.out";

    /* Create file */
    if ((fp = Util_fopen(bm, fname, "w")) == NULL)
        quit("initCloseKinHSFile: Can't open %s\n", fname);

    /* File content - nothing yet as added in another routine */
    fprintf(fp, "Time year species sex\n");

    /* Return file pointer */
    return (fp);
}

void writeCloseKinNcompsfiles(MSEBoxModel *bm, int year, int sp, FILE *fidPO, FILE *fidHS, double *****sim_nkin_HS, double *****sim_nkin_PO, double ****ncomps_HS_yaya, double *****ncomps_PO_syaya) {
    int y, s, a, y2, a2, b1, b2;
    //int f;

    for (y = bm->CloseKinEst[sp].first_sy; y < bm->CloseKinEst[sp].last_sy; y++) {
        for (a = 0; a < bm->RBCestimation.RBCspeciesParam[sp][MaxAge_id]; a++) {
             for (y2 = bm->CloseKinEst[sp].first_sy; y2 < bm->CloseKinEst[sp].last_sy; y2++) {
                 for (a2 = 0; a2 < bm->RBCestimation.RBCspeciesParam[sp][MaxAge_id]; a2++) {
                     for (s = 0; s < bm->RBCestimation.RBCspeciesParam[sp][Nsexes_id]; s++) {
                         b1 = y - a;
                         b2 = y2 - a2;

                         fprintf(fidPO, "%e %d %s %d %d %d %d %d %d %d\n", bm->dayt, year, FunctGroupArray[sp].groupCode, s, y, a, y2, a2, ((int)ncomps_PO_syaya[s][y][a][y2][a2]), ((int)sim_nkin_PO[s][y][a][y2][a2]));
                     
                         if (b1<b2) {
                             fprintf(fidHS, "%e %d %s %d %d %d %d %d %d %d %e\n", bm->dayt, year, FunctGroupArray[sp].groupCode, s, y, a, y2, a2, ((int)sim_nkin_HS[s][y][a][y2][a2]), ((int)sim_nkin_HS[s][y][a][y2][a2]), ncomps_HS_yaya[y][a][y2][a2]);
                        } else {
                            fprintf(fidHS, "%e %d %s %d %d %d %d %d %d %d %e\n", bm->dayt, year, FunctGroupArray[sp].groupCode, s, y2, a2, y, a, ((int)sim_nkin_HS[s][y][a][y2][a2]), ((int)sim_nkin_HS[s][y][a][y2][a2]), ncomps_HS_yaya[y][a][y2][a2]);
                            
                            //TODO: ASK Robin if this should be ((int)sim_nkin_HS[s][y2][a2][y][a]), ((int)sim_nkin_HS[s][y2][a2][y][a]).....
                        }
                    }
                }
            }
        }
    }

    return;
}

//******************************************************************************
//
// Name:  poissonRandom
// Description:
// calls :
// called by: CKSimulator()
// created  :
// based on: https://en.wikipedia.org/wiki/Poisson_distribution#Generating_Poisson-distributed_random_variables
//
//******************************************************************************
int poissonRandom(double lambda) {
  int n = 0; //counter of iteration
  double limit;
  double x;  //pseudo random number
    
  limit = exp(-lambda);
  x = (double) rand() / RAND_MAX;
  while (x > limit) {
    n++;
    x *= (double) rand() / RAND_MAX;
  }
  return n;
}
