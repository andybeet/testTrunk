/**
 * The purpose of the code in this file is to test the code that writes out SS3 input files, runs SS3 and reads back in the input files.
 * This is designed to be run with dummy values simply to test that the format of the files is correct and can be read in by SS3.
 *
 *	The code sets up dummy array structures that are used by the SS3 code.
 *
 *
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <signal.h>
#include <math.h>
#include <ctype.h>
#include <time.h>
#include <string.h>
#include <sjwlib.h>
#include <atlantisboxmodel.h>
#include <atUtilLib.h>

/**
 * When this function is called it is assumed that all arrays are allocated and zeroed.
 * This function is simply used to set up values that can be used to write out SS3 input files so we can test that the format is correct.
 *
 *
 */
int SS3_Test_Init(MSEBoxModel *bm) {

#ifndef _WIN32
	int groupIndex = 2;
	int sumregion = 0;
	int nYears;
	int AccumAge;

	bm->RBCestimation.RBCspeciesParam[groupIndex][HistYrMin_id] = 1915;
	bm->RBCestimation.RBCspeciesParam[groupIndex][HistYrMax_id] = 2008;

	nYears = bm->RBCestimation.RBCspeciesParam[groupIndex][HistYrMax_id] - bm->RBCestimation.RBCspeciesParam[groupIndex][HistYrMin_id];

	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[0][sumregion][0] = 445.6;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[1][sumregion][0] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[2][sumregion][0] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[3][sumregion][0] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[0][sumregion][1] = 447.4;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[1][sumregion][1] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[2][sumregion][1] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[3][sumregion][1] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[0][sumregion][2] = 518.6;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[1][sumregion][2] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[2][sumregion][2] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[3][sumregion][2] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[0][sumregion][3] = 805.7;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[1][sumregion][3] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[2][sumregion][3] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[3][sumregion][3] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[0][sumregion][4] = 1381.4;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[1][sumregion][4] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[2][sumregion][4] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[3][sumregion][4] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[0][sumregion][5] = 1117.3;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[1][sumregion][5] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[2][sumregion][5] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[3][sumregion][5] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[0][sumregion][6] = 1556.8;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[1][sumregion][6] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[2][sumregion][6] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[3][sumregion][6] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[0][sumregion][7] = 1007.8;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[1][sumregion][7] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[2][sumregion][7] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[3][sumregion][7] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[0][sumregion][8] = 954.6;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[1][sumregion][8] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[2][sumregion][8] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[3][sumregion][8] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[0][sumregion][9] = 1627.7;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[1][sumregion][9] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[2][sumregion][9] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[3][sumregion][9] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[0][sumregion][10] = 2363.3;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[1][sumregion][10] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[2][sumregion][10] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[3][sumregion][10] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[0][sumregion][11] = 2600.1;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[1][sumregion][11] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[2][sumregion][11] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[3][sumregion][11] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[0][sumregion][12] = 3281.5;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[1][sumregion][12] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[2][sumregion][12] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[3][sumregion][12] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[0][sumregion][13] = 3932.9;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[1][sumregion][13] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[2][sumregion][13] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[3][sumregion][13] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[0][sumregion][14] = 4521.7;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[1][sumregion][14] = 121.8;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[2][sumregion][14] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[3][sumregion][14] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[0][sumregion][15] = 3994.4;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[1][sumregion][15] = 396;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[2][sumregion][15] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[3][sumregion][15] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[0][sumregion][16] = 3518.4;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[1][sumregion][16] = 4.7;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[2][sumregion][16] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[3][sumregion][16] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[0][sumregion][17] = 3170.1;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[1][sumregion][17] = 462;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[2][sumregion][17] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[3][sumregion][17] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[0][sumregion][18] = 2946.8;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[1][sumregion][18] = 52.5;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[2][sumregion][18] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[3][sumregion][18] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[0][sumregion][19] = 2733.7;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[1][sumregion][19] = 331.7;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[2][sumregion][19] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[3][sumregion][19] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[0][sumregion][20] = 3016.4;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[1][sumregion][20] = 324.5;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[2][sumregion][20] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[3][sumregion][20] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[0][sumregion][21] = 3254.8;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[1][sumregion][21] = 1046.2;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[2][sumregion][21] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[3][sumregion][21] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[0][sumregion][22] = 3260.6;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[1][sumregion][22] = 764.8;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[2][sumregion][22] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[3][sumregion][22] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[0][sumregion][23] = 2753;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[1][sumregion][23] = 869.5;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[2][sumregion][23] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[3][sumregion][23] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[0][sumregion][24] = 2227.5;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[1][sumregion][24] = 1241.4;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[2][sumregion][24] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[3][sumregion][24] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[0][sumregion][25] = 977.9;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[1][sumregion][25] = 1329.3;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[2][sumregion][25] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[3][sumregion][25] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[0][sumregion][26] = 483;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[1][sumregion][26] = 1506.1;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[2][sumregion][26] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[3][sumregion][26] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[0][sumregion][27] = 200.8;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[1][sumregion][27] = 269.8;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[2][sumregion][27] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[3][sumregion][27] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[0][sumregion][28] = 267.1;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[1][sumregion][28] = 380.7;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[2][sumregion][28] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[3][sumregion][28] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[0][sumregion][29] = 377.9;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[1][sumregion][29] = 3149.1;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[2][sumregion][29] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[3][sumregion][29] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[0][sumregion][30] = 1143.8;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[1][sumregion][30] = 2601.6;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[2][sumregion][30] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[3][sumregion][30] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[0][sumregion][31] = 1305.6;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[1][sumregion][31] = 1611.3;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[2][sumregion][31] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[3][sumregion][31] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[0][sumregion][32] = 974.7;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[1][sumregion][32] = 1009.9;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[2][sumregion][32] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[3][sumregion][32] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[0][sumregion][33] = 812.7;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[1][sumregion][33] = 1064.4;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[2][sumregion][33] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[3][sumregion][33] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[0][sumregion][34] = 396.4;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[1][sumregion][34] = 868.2;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[2][sumregion][34] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[3][sumregion][34] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[0][sumregion][35] = 372.1;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[1][sumregion][35] = 1314.3;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[2][sumregion][35] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[3][sumregion][35] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[0][sumregion][36] = 699.4;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[1][sumregion][36] = 1950.2;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[2][sumregion][36] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[3][sumregion][36] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[0][sumregion][37] = 923.2;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[1][sumregion][37] = 1799;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[2][sumregion][37] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[3][sumregion][37] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[0][sumregion][38] = 620.3;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[1][sumregion][38] = 2682;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[2][sumregion][38] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[3][sumregion][38] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[0][sumregion][39] = 438.9;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[1][sumregion][39] = 2083.9;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[2][sumregion][39] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[3][sumregion][39] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[0][sumregion][40] = 253;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[1][sumregion][40] = 2318.7;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[2][sumregion][40] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[3][sumregion][40] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[0][sumregion][41] = 187.8;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[1][sumregion][41] = 2097.7;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[2][sumregion][41] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[3][sumregion][41] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[0][sumregion][42] = 167.2;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[1][sumregion][42] = 1525.3;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[2][sumregion][42] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[3][sumregion][42] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[0][sumregion][43] = 81.5;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[1][sumregion][43] = 1261.6;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[2][sumregion][43] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[3][sumregion][43] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[0][sumregion][44] = 38.7;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[1][sumregion][44] = 1525.9;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[2][sumregion][44] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[3][sumregion][44] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[0][sumregion][45] = 18.3;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[1][sumregion][45] = 2143.2;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[2][sumregion][45] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[3][sumregion][45] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[0][sumregion][46] = 10.5;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[1][sumregion][46] = 1035.5;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[2][sumregion][46] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[3][sumregion][46] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[0][sumregion][47] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[1][sumregion][47] = 1742;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[2][sumregion][47] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[3][sumregion][47] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[0][sumregion][48] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[1][sumregion][48] = 2956.2;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[2][sumregion][48] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[3][sumregion][48] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[0][sumregion][49] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[1][sumregion][49] = 2793;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[2][sumregion][49] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[3][sumregion][49] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[0][sumregion][50] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[1][sumregion][50] = 2850.6;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[2][sumregion][50] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[3][sumregion][50] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[0][sumregion][51] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[1][sumregion][51] = 2762;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[2][sumregion][51] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[3][sumregion][51] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[0][sumregion][52] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[1][sumregion][52] = 2867;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[2][sumregion][52] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[3][sumregion][52] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[0][sumregion][53] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[1][sumregion][53] = 2355;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[2][sumregion][53] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[3][sumregion][53] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[0][sumregion][54] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[1][sumregion][54] = 3289;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[2][sumregion][54] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[3][sumregion][54] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[0][sumregion][55] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[1][sumregion][55] = 2667;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[2][sumregion][55] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[3][sumregion][55] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[0][sumregion][56] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[1][sumregion][56] = 1792.5;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[2][sumregion][56] = 285.5;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[3][sumregion][56] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[0][sumregion][57] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[1][sumregion][57] = 1981;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[2][sumregion][57] = 491;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[3][sumregion][57] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[0][sumregion][58] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[1][sumregion][58] = 2397.1;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[2][sumregion][58] = 489.9;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[3][sumregion][58] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[0][sumregion][59] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[1][sumregion][59] = 1493.2;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[2][sumregion][59] = 368.8;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[3][sumregion][59] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[0][sumregion][60] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[1][sumregion][60] = 1367.4;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[2][sumregion][60] = 826.6;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[3][sumregion][60] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[0][sumregion][61] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[1][sumregion][61] = 899.5;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[2][sumregion][61] = 711.5;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[3][sumregion][61] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[0][sumregion][62] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[1][sumregion][62] = 976.8;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[2][sumregion][62] = 522.2;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[3][sumregion][62] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[0][sumregion][63] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[1][sumregion][63] = 836;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[2][sumregion][63] = 446;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[3][sumregion][63] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[0][sumregion][64] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[1][sumregion][64] = 927.8;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[2][sumregion][64] = 520.2;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[3][sumregion][64] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[0][sumregion][65] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[1][sumregion][65] = 851;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[2][sumregion][65] = 609;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[3][sumregion][65] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[0][sumregion][66] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[1][sumregion][66] = 418;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[2][sumregion][66] = 877;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[3][sumregion][66] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[0][sumregion][67] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[1][sumregion][67] = 615;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[2][sumregion][67] = 930;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[3][sumregion][67] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[0][sumregion][68] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[1][sumregion][68] = 889;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[2][sumregion][68] = 950;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[3][sumregion][68] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[0][sumregion][69] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[1][sumregion][69] = 890.3;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[2][sumregion][69] = 977.6;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[3][sumregion][69] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[0][sumregion][70] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[1][sumregion][70] = 890.3;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[2][sumregion][70] = 977.6;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[3][sumregion][70] = 30;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[0][sumregion][71] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[1][sumregion][71] = 891.5;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[2][sumregion][71] = 1005.2;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[3][sumregion][71] = 25.9;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[0][sumregion][72] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[1][sumregion][72] = 1357.6;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[2][sumregion][72] = 1108.6;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[3][sumregion][72] = 6;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[0][sumregion][73] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[1][sumregion][73] = 1176.6;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[2][sumregion][73] = 1262.6;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[3][sumregion][73] = 115.7;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[0][sumregion][74] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[1][sumregion][74] = 1189.4;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[2][sumregion][74] = 1318;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[3][sumregion][74] = 128.3;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[0][sumregion][75] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[1][sumregion][75] = 591.4;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[2][sumregion][75] = 1425;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[3][sumregion][75] = 178.2;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[0][sumregion][76] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[1][sumregion][76] = 745.5;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[2][sumregion][76] = 1461.3;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[3][sumregion][76] = 165.6;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[0][sumregion][77] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[1][sumregion][77] = 1019.4;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[2][sumregion][77] = 1080.4;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[3][sumregion][77] = 170;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[0][sumregion][78] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[1][sumregion][78] = 515.6;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[2][sumregion][78] = 961.9;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[3][sumregion][78] = 193.7;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[0][sumregion][79] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[1][sumregion][79] = 626.2;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[2][sumregion][79] = 981.5;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[3][sumregion][79] = 178;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[0][sumregion][80] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[1][sumregion][80] = 563.9;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[2][sumregion][80] = 1188.7;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[3][sumregion][80] = 139.3;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[0][sumregion][81] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[1][sumregion][81] = 711.3;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[2][sumregion][81] = 1264.9;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[3][sumregion][81] = 113.8;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[0][sumregion][82] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[1][sumregion][82] = 1023.3;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[2][sumregion][82] = 1541.8;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[3][sumregion][82] = 174.9;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[0][sumregion][83] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[1][sumregion][83] = 904.9;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[2][sumregion][83] = 1699.7;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[3][sumregion][83] = 185.6;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[0][sumregion][84] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[1][sumregion][84] = 1873.4;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[2][sumregion][84] = 1520.4;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[3][sumregion][84] = 248.2;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[0][sumregion][85] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[1][sumregion][85] = 1286;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[2][sumregion][85] = 2006.3;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[3][sumregion][85] = 202.5;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[0][sumregion][86] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[1][sumregion][86] = 1323.6;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[2][sumregion][86] = 1552.4;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[3][sumregion][86] = 116.1;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[0][sumregion][87] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[1][sumregion][87] = 1357.2;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[2][sumregion][87] = 1677.6;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[3][sumregion][87] = 237.4;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[0][sumregion][88] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[1][sumregion][88] = 1315.3;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[2][sumregion][88] = 2077.3;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[3][sumregion][88] = 276.8;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[0][sumregion][89] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[1][sumregion][89] = 1335.7;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[2][sumregion][89] = 1705.6;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[3][sumregion][89] = 454.5;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[0][sumregion][90] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[1][sumregion][90] = 1413.4;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[2][sumregion][90] = 1420.7;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[3][sumregion][90] = 394;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[0][sumregion][91] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[1][sumregion][91] = 1238.2;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[2][sumregion][91] = 1434.8;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[3][sumregion][91] = 344.3;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[0][sumregion][92] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[1][sumregion][92] = 1489.8;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[2][sumregion][92] = 1357.4;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[3][sumregion][92] = 205.1;

//	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[0][sumregion][maxyr] = 0;
//	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[1][sumregion][maxyr] = 1366.7;
//	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[2][sumregion][maxyr] = 1245.2;
//	bm->RBCestimation.RBCspeciesArray[groupIndex].CatchData[3][sumregion][maxyr] = 188.1;

	bm->RBCestimation.RBCspeciesParam[groupIndex][TAC_old_id] = 2800;

	bm->RBCestimation.RBCspeciesParam[groupIndex][HistYrMin_id] = 1915;
	bm->RBCestimation.RBCspeciesParam[groupIndex][HistYrMax_id] = 2008;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CurrentYear = 2007;

//	bm->RBCestimation.RBCspeciesArray[groupIndex].CPUEcv[0] = 0.5;
//	bm->RBCestimation.RBCspeciesArray[groupIndex].CPUEcv[1] = 0.2;
//	bm->RBCestimation.RBCspeciesArray[groupIndex].CPUEcv[2] = 0.05;
//	bm->RBCestimation.RBCspeciesArray[groupIndex].CPUEcv[3] = 0.4;
//
//	bm->RBCestimation.RBCspeciesArray[groupIndex].CPUEqmu[0] = 0.00004;
//	bm->RBCestimation.RBCspeciesArray[groupIndex].CPUEqmu[1] = 0.00006;
//	bm->RBCestimation.RBCspeciesArray[groupIndex].CPUEqmu[2] = 0.00009;
//	bm->RBCestimation.RBCspeciesArray[groupIndex].CPUEqmu[3] = 0.0005;
//
//	bm->RBCestimation.RBCspeciesArray[groupIndex].CPUEpow[0] = 1;
//	bm->RBCestimation.RBCspeciesArray[groupIndex].CPUEpow[1] = 1;
//	bm->RBCestimation.RBCspeciesArray[groupIndex].CPUEpow[2] = 1;
//	bm->RBCestimation.RBCspeciesArray[groupIndex].CPUEpow[3] = 1;
//
//	bm->RBCestimation.RBCspeciesArray[groupIndex].CPUEvar[0] = 0;
//	bm->RBCestimation.RBCspeciesArray[groupIndex].CPUEvar[1] = 0;
//	bm->RBCestimation.RBCspeciesArray[groupIndex].CPUEvar[2] = 0;
//	bm->RBCestimation.RBCspeciesArray[groupIndex].CPUEvar[3] = 0;
//
//	bm->RBCestimation.RBCspeciesArray[groupIndex].CPUEcorr[0] = 1;
//	bm->RBCestimation.RBCspeciesArray[groupIndex].CPUEcorr[1] = 1;
//	bm->RBCestimation.RBCspeciesArray[groupIndex].CPUEcorr[2] = 1;
//	bm->RBCestimation.RBCspeciesArray[groupIndex].CPUEcorr[3] = 1;

	bm->RBCestimation.RBCspeciesArray[groupIndex].CvLA0[0][0] = 0.103;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CvLAmax[0][0] = 0.103;

	bm->RBCestimation.RBCspeciesArray[groupIndex].CPUEgen[0][sumregion][0] = 0.97;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CPUEgen[0][sumregion][1] = 1.03;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CPUEgen[0][sumregion][2] = 0.5;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CPUEgen[0][sumregion][3] = 0.54;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CPUEgen[0][sumregion][4] = 0.8;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CPUEgen[0][sumregion][5] = 0.89;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CPUEgen[0][sumregion][6] = 1.01;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CPUEgen[0][sumregion][7] = 0.67;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CPUEgen[0][sumregion][8] = 0.98;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CPUEgen[0][sumregion][9] = 1.03;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CPUEgen[0][sumregion][10] = 0.68;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CPUEgen[0][sumregion][11] = 1.22;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CPUEgen[0][sumregion][12] = 0.52;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CPUEgen[0][sumregion][13] = 0.36;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CPUEgen[0][sumregion][14] = 0.37;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CPUEgen[0][sumregion][15] = 0.42;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CPUEgen[0][sumregion][16] = 0.45;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CPUEgen[0][sumregion][17] = 0.36;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CPUEgen[0][sumregion][18] = 0.28;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CPUEgen[0][sumregion][19] = 0.31;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CPUEgen[0][sumregion][20] = 0.4;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CPUEgen[0][sumregion][21] = 0.23;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CPUEgen[0][sumregion][22] = 0.32;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CPUEgen[0][sumregion][23] = 0.19;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CPUEgen[0][sumregion][24] = 0.14;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CPUEgen[0][sumregion][25] = 0.18;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CPUEgen[0][sumregion][26] = 0.08;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CPUEgen[0][sumregion][27] = 0.14;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CPUEgen[0][sumregion][28] = 0.17;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CPUEgen[0][sumregion][29] = 0.41;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CPUEgen[0][sumregion][30] = 0.24;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CPUEgen[0][sumregion][31] = 0.19;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CPUEgen[0][sumregion][32] = 0.19;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CPUEgen[0][sumregion][33] = 0.07;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CPUEgen[0][sumregion][34] = 0.09;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CPUEgen[0][sumregion][35] = 0.18;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CPUEgen[0][sumregion][36] = 0.17;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CPUEgen[0][sumregion][37] = 0.13;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CPUEgen[0][sumregion][38] = 0.29;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CPUEgen[0][sumregion][39] = 0.24;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CPUEgen[0][sumregion][40] = 0.26;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CPUEgen[0][sumregion][41] = 0.32;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CPUEgen[0][sumregion][42] = 0.08;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CPUEgen[0][sumregion][43] = 0.17;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CPUEgen[0][sumregion][44] = 0.25;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CPUEgen[0][sumregion][45] = 0.14;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CPUEgen[0][sumregion][46] = 0.23;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CPUEgen[1][sumregion][71] = 0.43;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CPUEgen[1][sumregion][72] = 0.5;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CPUEgen[1][sumregion][73] = 0.59;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CPUEgen[1][sumregion][74] = 0.87;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CPUEgen[1][sumregion][75] = 0.73;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CPUEgen[1][sumregion][76] = 0.63;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CPUEgen[1][sumregion][77] = 0.46;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CPUEgen[1][sumregion][78] = 0.49;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CPUEgen[1][sumregion][79] = 0.42;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CPUEgen[1][sumregion][80] = 0.5;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CPUEgen[1][sumregion][81] = 0.6;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CPUEgen[1][sumregion][82] = 0.58;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CPUEgen[1][sumregion][83] = 0.46;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CPUEgen[1][sumregion][84] = 0.42;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CPUEgen[1][sumregion][85] = 0.44;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CPUEgen[1][sumregion][86] = 0.47;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CPUEgen[1][sumregion][87] = 0.47;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CPUEgen[1][sumregion][88] = 0.5;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CPUEgen[1][sumregion][89] = 0.31;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CPUEgen[1][sumregion][90] = 0.42;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CPUEgen[1][sumregion][91] = 0.41;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CPUEgen[1][sumregion][92] = 0.48;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CPUEgen[2][sumregion][71] = 0.54;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CPUEgen[2][sumregion][72] = 0.61;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CPUEgen[2][sumregion][73] = 0.68;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CPUEgen[2][sumregion][74] = 0.72;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CPUEgen[2][sumregion][75] = 0.67;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CPUEgen[2][sumregion][76] = 0.7;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CPUEgen[2][sumregion][77] = 0.63;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CPUEgen[2][sumregion][78] = 0.64;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CPUEgen[2][sumregion][79] = 0.66;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CPUEgen[2][sumregion][80] = 0.63;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CPUEgen[2][sumregion][81] = 0.66;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CPUEgen[2][sumregion][82] = 0.59;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CPUEgen[2][sumregion][83] = 0.54;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CPUEgen[2][sumregion][84] = 0.53;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CPUEgen[2][sumregion][85] = 0.48;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CPUEgen[2][sumregion][86] = 0.49;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CPUEgen[2][sumregion][87] = 0.52;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CPUEgen[2][sumregion][88] = 0.48;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CPUEgen[2][sumregion][89] = 0.45;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CPUEgen[2][sumregion][90] = 0.45;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CPUEgen[2][sumregion][91] = 0.46;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CPUEgen[2][sumregion][92] = 0.47;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CPUEgen[3][sumregion][71] = 0.73;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CPUEgen[3][sumregion][72] = 0.58;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CPUEgen[3][sumregion][73] = 0.51;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CPUEgen[3][sumregion][74] = 1.19;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CPUEgen[3][sumregion][75] = 0.83;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CPUEgen[3][sumregion][76] = 0.75;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CPUEgen[3][sumregion][77] = 1.05;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CPUEgen[3][sumregion][78] = 0.85;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CPUEgen[3][sumregion][79] = 0.65;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CPUEgen[3][sumregion][80] = 0.57;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CPUEgen[3][sumregion][81] = 0.66;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CPUEgen[3][sumregion][82] = 0.96;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CPUEgen[3][sumregion][83] = 1.11;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CPUEgen[3][sumregion][84] = 0.85;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CPUEgen[3][sumregion][85] = 0.78;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CPUEgen[3][sumregion][86] = 0.56;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CPUEgen[3][sumregion][87] = 0.68;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CPUEgen[3][sumregion][88] = 0.59;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CPUEgen[3][sumregion][89] = 0.58;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CPUEgen[3][sumregion][90] = 0.51;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CPUEgen[3][sumregion][91] = 1.1;
	bm->RBCestimation.RBCspeciesArray[groupIndex].CPUEgen[3][sumregion][92] = 0.71;

//	bm->RBCestimation.RBCspeciesArray[groupIndex].DiscCV[0] = 0.1;
//	bm->RBCestimation.RBCspeciesArray[groupIndex].DiscCV[1] = 0.1;
//	bm->RBCestimation.RBCspeciesArray[groupIndex].DiscCV[2] = 0.1;
//	bm->RBCestimation.RBCspeciesArray[groupIndex].DiscCV[3] = 0.1;

	bm->RBCestimation.RBCspeciesArray[groupIndex].DiscData[1][sumregion][83][0] = 0.01829;
	bm->RBCestimation.RBCspeciesArray[groupIndex].DiscData[1][sumregion][84][0] = 0.01567;
	bm->RBCestimation.RBCspeciesArray[groupIndex].DiscData[1][sumregion][85][0] = 0.02565;
	bm->RBCestimation.RBCspeciesArray[groupIndex].DiscData[1][sumregion][86][0] = 0.02594;
	bm->RBCestimation.RBCspeciesArray[groupIndex].DiscData[1][sumregion][87][0] = 0.01785;
	bm->RBCestimation.RBCspeciesArray[groupIndex].DiscData[1][sumregion][88][0] = 0.02274;
	bm->RBCestimation.RBCspeciesArray[groupIndex].DiscData[1][sumregion][89][0] = 0.02394;
	bm->RBCestimation.RBCspeciesArray[groupIndex].DiscData[1][sumregion][90][0] = 0.02561;
	bm->RBCestimation.RBCspeciesArray[groupIndex].DiscData[1][sumregion][91][0] = 0.022;
	bm->RBCestimation.RBCspeciesArray[groupIndex].DiscData[1][sumregion][92][0] = 0.02359;
	bm->RBCestimation.RBCspeciesArray[groupIndex].DiscData[2][sumregion][81][0] = 0.02648;
	bm->RBCestimation.RBCspeciesArray[groupIndex].DiscData[2][sumregion][82][0] = 0.02303;
	bm->RBCestimation.RBCspeciesArray[groupIndex].DiscData[2][sumregion][83][0] = 0.02522;
	bm->RBCestimation.RBCspeciesArray[groupIndex].DiscData[2][sumregion][84][0] = 0.02711;
	bm->RBCestimation.RBCspeciesArray[groupIndex].DiscData[2][sumregion][85][0] = 0.03257;
	bm->RBCestimation.RBCspeciesArray[groupIndex].DiscData[2][sumregion][86][0] = 0.03316;
	bm->RBCestimation.RBCspeciesArray[groupIndex].DiscData[2][sumregion][87][0] = 0.03445;
	bm->RBCestimation.RBCspeciesArray[groupIndex].DiscData[2][sumregion][88][0] = 0.03355;
	bm->RBCestimation.RBCspeciesArray[groupIndex].DiscData[2][sumregion][89][0] = 0.03554;
	bm->RBCestimation.RBCspeciesArray[groupIndex].DiscData[2][sumregion][90][0] = 0.03937;
	bm->RBCestimation.RBCspeciesArray[groupIndex].DiscData[2][sumregion][91][0] = 0.04319;
	bm->RBCestimation.RBCspeciesArray[groupIndex].DiscData[2][sumregion][92][0] = 0.03752;
	bm->RBCestimation.RBCspeciesArray[groupIndex].DiscData[3][sumregion][79][0] = 0.01853;
	bm->RBCestimation.RBCspeciesArray[groupIndex].DiscData[3][sumregion][82][0] = 0.01654;
	bm->RBCestimation.RBCspeciesArray[groupIndex].DiscData[3][sumregion][83][0] = 0.02069;
	bm->RBCestimation.RBCspeciesArray[groupIndex].DiscData[3][sumregion][84][0] = 0.02496;
	bm->RBCestimation.RBCspeciesArray[groupIndex].DiscData[3][sumregion][87][0] = 0.03983;
	bm->RBCestimation.RBCspeciesArray[groupIndex].DiscData[3][sumregion][88][0] = 0.03459;
	bm->RBCestimation.RBCspeciesArray[groupIndex].DiscData[3][sumregion][90][0] = 0.05324;
	bm->RBCestimation.RBCspeciesArray[groupIndex].DiscData[3][sumregion][91][0] = 0.05653;

	bm->RBCestimation.RBCspeciesParam[groupIndex][Nlen_id] = 54;
//
//	bm->RBCestimation.RBCspeciesArray[groupIndex].LoLenBin = Util_Alloc_Init_1D_Double((int) bm->RBCestimation.RBCspeciesParam[groupIndex][Nlen_id], 0.0);
//
//	bm->RBCestimation.RBCspeciesArray[groupIndex].LoLenBin[0] = 5;
//	bm->RBCestimation.RBCspeciesArray[groupIndex].LoLenBin[1] = 6;
//	bm->RBCestimation.RBCspeciesArray[groupIndex].LoLenBin[2] = 7;
//	bm->RBCestimation.RBCspeciesArray[groupIndex].LoLenBin[3] = 8;
//	bm->RBCestimation.RBCspeciesArray[groupIndex].LoLenBin[4] = 9;
//	bm->RBCestimation.RBCspeciesArray[groupIndex].LoLenBin[5] = 10;
//	bm->RBCestimation.RBCspeciesArray[groupIndex].LoLenBin[6] = 11;
//	bm->RBCestimation.RBCspeciesArray[groupIndex].LoLenBin[7] = 12;
//	bm->RBCestimation.RBCspeciesArray[groupIndex].LoLenBin[8] = 13;
//	bm->RBCestimation.RBCspeciesArray[groupIndex].LoLenBin[9] = 14;
//	bm->RBCestimation.RBCspeciesArray[groupIndex].LoLenBin[10] = 15;
//	bm->RBCestimation.RBCspeciesArray[groupIndex].LoLenBin[11] = 16;
//	bm->RBCestimation.RBCspeciesArray[groupIndex].LoLenBin[12] = 17;
//	bm->RBCestimation.RBCspeciesArray[groupIndex].LoLenBin[13] = 18;
//	bm->RBCestimation.RBCspeciesArray[groupIndex].LoLenBin[14] = 19;
//	bm->RBCestimation.RBCspeciesArray[groupIndex].LoLenBin[15] = 20;
//	bm->RBCestimation.RBCspeciesArray[groupIndex].LoLenBin[16] = 21;
//	bm->RBCestimation.RBCspeciesArray[groupIndex].LoLenBin[17] = 22;
//	bm->RBCestimation.RBCspeciesArray[groupIndex].LoLenBin[18] = 23;
//	bm->RBCestimation.RBCspeciesArray[groupIndex].LoLenBin[19] = 24;
//	bm->RBCestimation.RBCspeciesArray[groupIndex].LoLenBin[20] = 25;
//	bm->RBCestimation.RBCspeciesArray[groupIndex].LoLenBin[21] = 26;
//	bm->RBCestimation.RBCspeciesArray[groupIndex].LoLenBin[22] = 27;
//	bm->RBCestimation.RBCspeciesArray[groupIndex].LoLenBin[23] = 28;
//	bm->RBCestimation.RBCspeciesArray[groupIndex].LoLenBin[24] = 29;
//	bm->RBCestimation.RBCspeciesArray[groupIndex].LoLenBin[25] = 30;
//	bm->RBCestimation.RBCspeciesArray[groupIndex].LoLenBin[26] = 31;
//	bm->RBCestimation.RBCspeciesArray[groupIndex].LoLenBin[27] = 32;
//	bm->RBCestimation.RBCspeciesArray[groupIndex].LoLenBin[28] = 33;
//	bm->RBCestimation.RBCspeciesArray[groupIndex].LoLenBin[29] = 34;
//	bm->RBCestimation.RBCspeciesArray[groupIndex].LoLenBin[30] = 35;
//	bm->RBCestimation.RBCspeciesArray[groupIndex].LoLenBin[31] = 36;
//	bm->RBCestimation.RBCspeciesArray[groupIndex].LoLenBin[32] = 37;
//	bm->RBCestimation.RBCspeciesArray[groupIndex].LoLenBin[33] = 38;
//	bm->RBCestimation.RBCspeciesArray[groupIndex].LoLenBin[34] = 39;
//	bm->RBCestimation.RBCspeciesArray[groupIndex].LoLenBin[35] = 40;
//	bm->RBCestimation.RBCspeciesArray[groupIndex].LoLenBin[36] = 41;
//	bm->RBCestimation.RBCspeciesArray[groupIndex].LoLenBin[37] = 42;
//	bm->RBCestimation.RBCspeciesArray[groupIndex].LoLenBin[38] = 43;
//	bm->RBCestimation.RBCspeciesArray[groupIndex].LoLenBin[39] = 44;
//	bm->RBCestimation.RBCspeciesArray[groupIndex].LoLenBin[40] = 45;
//	bm->RBCestimation.RBCspeciesArray[groupIndex].LoLenBin[41] = 46;
//	bm->RBCestimation.RBCspeciesArray[groupIndex].LoLenBin[42] = 47;
//	bm->RBCestimation.RBCspeciesArray[groupIndex].LoLenBin[43] = 48;
//	bm->RBCestimation.RBCspeciesArray[groupIndex].LoLenBin[44] = 49;
//	bm->RBCestimation.RBCspeciesArray[groupIndex].LoLenBin[45] = 50;
//	bm->RBCestimation.RBCspeciesArray[groupIndex].LoLenBin[46] = 51;
//	bm->RBCestimation.RBCspeciesArray[groupIndex].LoLenBin[47] = 52;
//	bm->RBCestimation.RBCspeciesArray[groupIndex].LoLenBin[48] = 53;
//	bm->RBCestimation.RBCspeciesArray[groupIndex].LoLenBin[49] = 54;
//	bm->RBCestimation.RBCspeciesArray[groupIndex].LoLenBin[50] = 55;
//	bm->RBCestimation.RBCspeciesArray[groupIndex].LoLenBin[51] = 56;
//	bm->RBCestimation.RBCspeciesArray[groupIndex].LoLenBin[52] = 57;
//	bm->RBCestimation.RBCspeciesArray[groupIndex].LoLenBin[53] = 58;

	// Retained.

	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][0][0] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][1][0] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][2][0] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][3][0] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][0][1] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][1][1] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][2][1] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][3][1] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][0][2] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][1][2] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][2][2] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][3][2] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][0][3] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][1][3] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][2][3] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][3][3] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][0][4] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][1][4] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][2][4] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][3][4] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][0][5] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][1][5] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][2][5] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][3][5] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][0][6] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][1][6] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][2][6] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][3][6] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][0][7] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][1][7] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][2][7] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][3][7] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][0][8] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][1][8] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][2][8] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][3][8] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][0][9] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][1][9] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][2][9] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][3][9] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][0][10] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][1][10] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][2][10] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][3][10] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][0][11] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][1][11] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][2][11] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][3][11] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][0][12] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][1][12] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][2][12] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][3][12] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][0][13] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][1][13] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][2][13] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][3][13] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][0][14] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][1][14] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][2][14] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][3][14] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][0][15] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][1][15] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][2][15] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][3][15] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][0][16] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][1][16] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][2][16] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][3][16] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][0][17] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][1][17] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][2][17] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][3][17] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][0][18] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][1][18] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][2][18] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][3][18] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][0][19] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][1][19] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][2][19] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][3][19] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][0][20] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][1][20] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][2][20] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][3][20] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][0][21] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][1][21] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][2][21] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][3][21] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][0][22] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][1][22] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][2][22] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][3][22] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][0][23] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][1][23] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][2][23] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][3][23] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][0][24] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][1][24] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][2][24] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][3][24] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][0][25] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][1][25] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][2][25] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][3][25] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][0][26] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][1][26] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][2][26] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][3][26] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][0][27] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][1][27] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][2][27] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][3][27] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][0][28] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][1][28] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][2][28] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][3][28] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][0][29] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][1][29] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][2][29] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][3][29] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][0][30] = 1;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][1][30] = 1;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][2][30] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][3][30] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][0][31] = 1;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][1][31] = 1;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][2][31] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][3][31] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][0][32] = 1;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][1][32] = 1;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][2][32] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][3][32] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][0][33] = 1;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][1][33] = 1;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][2][33] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][3][33] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][0][34] = 1;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][1][34] = 1;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][2][34] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][3][34] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][0][35] = 1;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][1][35] = 1;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][2][35] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][3][35] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][0][36] = 1;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][1][36] = 1;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][2][36] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][3][36] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][0][37] = 1;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][1][37] = 1;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][2][37] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][3][37] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][0][38] = 1;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][1][38] = 1;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][2][38] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][3][38] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][0][39] = 1;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][1][39] = 1;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][2][39] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][3][39] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][0][40] = 1;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][1][40] = 1;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][2][40] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][3][40] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][0][41] = 1;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][1][41] = 1;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][2][41] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][3][41] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][0][42] = 1;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][1][42] = 1;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][2][42] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][3][42] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][0][43] = 1;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][1][43] = 1;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][2][43] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][3][43] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][0][44] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][1][44] = 1;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][2][44] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][3][44] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][0][45] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][1][45] = 1;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][2][45] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][3][45] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][0][46] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][1][46] = 1;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][2][46] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][3][46] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][0][47] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][1][47] = 1;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][2][47] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][3][47] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][0][48] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][1][48] = 1;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][2][48] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][3][48] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][0][49] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][1][49] = 1;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][2][49] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][3][49] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][0][50] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][1][50] = 1;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][2][50] = 1;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][3][50] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][0][51] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][1][51] = 1;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][2][51] = 1;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][3][51] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][0][52] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][1][52] = 1;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][2][52] = 1;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][3][52] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][0][53] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][1][53] = 1;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][2][53] = 1;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][3][53] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][0][54] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][1][54] = 1;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][2][54] = 1;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][3][54] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][0][55] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][1][55] = 1;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][2][55] = 1;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][3][55] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][0][56] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][1][56] = 1;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][2][56] = 1;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][3][56] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][0][57] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][1][57] = 1;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][2][57] = 1;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][3][57] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][0][58] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][1][58] = 1;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][2][58] = 1;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][3][58] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][0][59] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][1][59] = 1;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][2][59] = 1;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][3][59] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][0][60] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][1][60] = 1;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][2][60] = 1;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][3][60] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][0][61] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][1][61] = 1;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][2][61] = 1;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][3][61] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][0][62] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][1][62] = 1;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][2][62] = 1;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][3][62] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][0][63] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][1][63] = 1;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][2][63] = 1;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][3][63] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][0][64] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][1][64] = 1;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][2][64] = 1;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][3][64] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][0][65] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][1][65] = 1;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][2][65] = 1;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][3][65] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][0][66] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][1][66] = 1;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][2][66] = 1;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][3][66] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][0][67] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][1][67] = 1;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][2][67] = 1;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][3][67] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][0][68] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][1][68] = 1;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][2][68] = 1;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][3][68] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][0][69] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][1][69] = 1;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][2][69] = 1;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][3][69] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][0][70] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][1][70] = 1;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][2][70] = 1;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][3][70] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][0][71] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][1][71] = 1;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][2][71] = 1;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][3][71] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][0][72] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][1][72] = 1;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][2][72] = 1;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][3][72] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][0][73] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][1][73] = 1;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][2][73] = 1;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][3][73] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][0][74] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][1][74] = 1;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][2][74] = 1;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][3][74] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][0][75] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][1][75] = 1;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][2][75] = 1;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][3][75] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][0][76] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][1][76] = 1;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][2][76] = 1;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][3][76] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][0][77] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][1][77] = 1;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][2][77] = 1;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][3][77] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][0][78] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][1][78] = 1;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][2][78] = 1;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][3][78] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][0][79] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][1][79] = 1;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][2][79] = 1;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][3][79] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][0][80] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][1][80] = 1;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][2][80] = 1;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][3][80] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][0][81] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][1][81] = 1;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][2][81] = 1;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][3][81] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][0][82] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][1][82] = 1;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][2][82] = 1;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][3][82] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][0][83] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][1][83] = 1;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][2][83] = 1;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][3][83] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][0][84] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][1][84] = 1;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][2][84] = 1;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][3][84] = 1;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][0][85] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][1][85] = 1;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][2][85] = 1;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][3][85] = 1;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][0][86] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][1][86] = 1;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][2][86] = 1;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][3][86] = 1;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][0][87] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][1][87] = 1;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][2][87] = 1;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][3][87] = 1;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][0][88] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][1][88] = 1;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][2][88] = 1;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][3][88] = 1;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][0][89] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][1][89] = 1;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][2][89] = 1;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][3][89] = 1;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][0][90] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][1][90] = 1;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][2][90] = 1;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][3][90] = 1;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][0][91] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][1][91] = 1;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][2][91] = 1;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][3][91] = 1;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][0][92] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][1][92] = 1;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][2][92] = 1;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][3][92] = 1;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][0][93] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][1][93] = 1;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][2][93] = 1;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[0][3][93] = 1;

	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][0][0] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][1][0] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][2][0] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][3][0] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][0][1] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][1][1] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][2][1] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][3][1] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][0][2] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][1][2] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][2][2] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][3][2] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][0][3] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][1][3] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][2][3] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][3][3] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][0][4] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][1][4] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][2][4] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][3][4] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][0][5] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][1][5] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][2][5] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][3][5] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][0][6] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][1][6] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][2][6] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][3][6] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][0][7] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][1][7] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][2][7] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][3][7] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][0][8] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][1][8] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][2][8] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][3][8] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][0][9] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][1][9] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][2][9] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][3][9] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][0][10] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][1][10] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][2][10] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][3][10] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][0][11] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][1][11] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][2][11] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][3][11] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][0][12] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][1][12] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][2][12] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][3][12] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][0][13] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][1][13] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][2][13] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][3][13] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][0][14] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][1][14] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][2][14] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][3][14] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][0][15] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][1][15] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][2][15] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][3][15] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][0][16] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][1][16] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][2][16] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][3][16] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][0][17] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][1][17] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][2][17] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][3][17] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][0][18] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][1][18] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][2][18] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][3][18] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][0][19] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][1][19] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][2][19] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][3][19] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][0][20] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][1][20] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][2][20] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][3][20] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][0][21] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][1][21] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][2][21] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][3][21] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][0][22] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][1][22] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][2][22] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][3][22] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][0][23] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][1][23] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][2][23] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][3][23] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][0][24] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][1][24] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][2][24] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][3][24] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][0][25] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][1][25] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][2][25] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][3][25] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][0][26] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][1][26] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][2][26] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][3][26] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][0][27] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][1][27] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][2][27] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][3][27] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][0][28] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][1][28] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][2][28] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][3][28] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][0][29] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][1][29] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][2][29] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][3][29] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][0][30] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][1][30] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][2][30] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][3][30] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][0][31] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][1][31] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][2][31] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][3][31] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][0][32] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][1][32] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][2][32] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][3][32] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][0][33] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][1][33] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][2][33] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][3][33] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][0][34] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][1][34] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][2][34] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][3][34] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][0][35] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][1][35] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][2][35] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][3][35] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][0][36] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][1][36] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][2][36] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][3][36] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][0][37] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][1][37] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][2][37] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][3][37] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][0][38] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][1][38] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][2][38] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][3][38] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][0][39] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][1][39] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][2][39] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][3][39] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][0][40] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][1][40] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][2][40] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][3][40] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][0][41] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][1][41] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][2][41] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][3][41] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][0][42] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][1][42] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][2][42] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][3][42] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][0][43] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][1][43] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][2][43] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][3][43] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][0][44] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][1][44] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][2][44] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][3][44] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][0][45] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][1][45] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][2][45] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][3][45] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][0][46] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][1][46] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][2][46] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][3][46] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][0][47] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][1][47] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][2][47] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][3][47] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][0][48] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][1][48] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][2][48] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][3][48] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][0][49] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][1][49] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][2][49] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][3][49] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][0][50] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][1][50] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][2][50] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][3][50] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][0][51] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][1][51] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][2][51] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][3][51] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][0][52] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][1][52] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][2][52] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][3][52] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][0][53] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][1][53] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][2][53] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][3][53] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][0][54] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][1][54] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][2][54] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][3][54] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][0][55] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][1][55] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][2][55] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][3][55] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][0][56] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][1][56] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][2][56] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][3][56] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][0][57] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][1][57] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][2][57] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][3][57] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][0][58] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][1][58] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][2][58] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][3][58] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][0][59] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][1][59] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][2][59] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][3][59] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][0][60] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][1][60] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][2][60] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][3][60] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][0][61] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][1][61] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][2][61] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][3][61] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][0][62] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][1][62] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][2][62] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][3][62] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][0][63] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][1][63] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][2][63] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][3][63] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][0][64] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][1][64] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][2][64] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][3][64] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][0][65] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][1][65] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][2][65] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][3][65] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][0][66] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][1][66] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][2][66] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][3][66] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][0][67] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][1][67] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][2][67] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][3][67] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][0][68] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][1][68] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][2][68] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][3][68] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][0][69] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][1][69] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][2][69] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][3][69] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][0][70] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][1][70] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][2][70] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][3][70] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][0][71] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][1][71] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][2][71] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][3][71] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][0][72] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][1][72] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][2][72] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][3][72] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][0][73] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][1][73] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][2][73] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][3][73] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][0][74] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][1][74] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][2][74] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][3][74] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][0][75] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][1][75] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][2][75] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][3][75] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][0][76] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][1][76] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][2][76] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][3][76] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][0][77] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][1][77] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][2][77] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][3][77] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][0][78] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][1][78] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][2][78] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][3][78] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][0][79] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][1][79] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][2][79] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][3][79] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][0][80] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][1][80] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][2][80] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][3][80] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][0][81] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][1][81] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][2][81] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][3][81] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][0][82] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][1][82] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][2][82] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][3][82] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][0][83] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][1][83] = 1;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][2][83] = 1;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][3][83] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][0][84] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][1][84] = 1;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][2][84] = 1;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][3][84] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][0][85] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][1][85] = 1;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][2][85] = 1;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][3][85] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][0][86] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][1][86] = 1;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][2][86] = 1;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][3][86] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][0][87] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][1][87] = 1;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][2][87] = 1;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][3][87] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][0][88] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][1][88] = 1;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][2][88] = 1;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][3][88] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][0][89] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][1][89] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][2][89] = 1;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][3][89] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][0][90] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][1][90] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][2][90] = 1;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][3][90] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][0][91] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][1][91] = 1;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][2][91] = 1;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][3][91] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][0][92] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][1][92] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][2][92] = 1;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][3][92] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][0][93] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][1][93] = 0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][2][93] = 1;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthFltYr[2][3][93] = 0;



	bm->RBCestimation.RBCspeciesArray[groupIndex].LenComp = Util_Alloc_Init_5D_Int(bm->RBCestimation.RBCspeciesParam[groupIndex][Nlen_id], 3, nYears,
			bm->K_num_sexes, bm->K_num_fisheries, 0);

	memcpy(bm->RBCestimation.RBCspeciesArray[groupIndex].LenComp[2][1][80][0], (int[]
				)
					{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 4, 4, 3, 9, 8, 12, 19, 16, 32, 16, 25, 22, 27, 18, 12, 14, 14, 7, 10, 6,
							5, 7, 6, 2, 0, 0, 1, 2, 1, 2, 2, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
							0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
			bm->RBCestimation.RBCspeciesParam[groupIndex][Nlen_id] * sizeof(int));
	memcpy(bm->RBCestimation.RBCspeciesArray[groupIndex].LenComp[2][1][81][0], (int[]
				)
					{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 2, 0, 0, 3, 5, 7, 4, 8, 12, 24, 32, 21, 17, 37, 30, 22, 19, 10, 17, 9, 12, 11,
							6, 2, 5, 3, 2, 3, 2, 2, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
							0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
			bm->RBCestimation.RBCspeciesParam[groupIndex][Nlen_id] * sizeof(int));
	memcpy(bm->RBCestimation.RBCspeciesArray[groupIndex].LenComp[2][1][82][0], (int[]
				)
					{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 1, 1, 4, 5, 1, 13, 15, 17, 20, 38, 23, 28, 26, 32, 16, 19, 14, 10, 13, 11,
							14, 8, 12, 5, 5, 2, 2, 0, 2, 1, 2, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
							0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
			bm->RBCestimation.RBCspeciesParam[groupIndex][Nlen_id] * sizeof(int));
	memcpy(bm->RBCestimation.RBCspeciesArray[groupIndex].LenComp[2][1][83][0], (int[]
				)
					{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 1, 3, 3, 6, 7, 13, 15, 11, 17, 24, 32, 23, 26, 14, 27, 10, 18, 10, 13, 7,
							3, 11, 8, 5, 4, 2, 3, 1, 2, 1, 1, 1, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
							0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
			bm->RBCestimation.RBCspeciesParam[groupIndex][Nlen_id] * sizeof(int));
	memcpy(bm->RBCestimation.RBCspeciesArray[groupIndex].LenComp[1][1][83][2], (int[]
				)
					{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 1, 2, 2, 2, 6, 2, 10, 15, 5, 13, 8, 4, 8, 11, 6, 2, 4, 1, 1, 0, 2, 0, 0, 0, 0,
							0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
							0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, bm->RBCestimation.RBCspeciesParam[groupIndex][Nlen_id] * sizeof(int));
	memcpy(bm->RBCestimation.RBCspeciesArray[groupIndex].LenComp[2][1][84][0], (int[]
				)
					{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 1, 0, 0, 2, 4, 7, 8, 13, 17, 26, 28, 29, 42, 37, 25, 36, 27, 34, 29, 20, 18, 13, 16,
							13, 7, 12, 11, 7, 3, 2, 1, 2, 0, 2, 1, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
							0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
			bm->RBCestimation.RBCspeciesParam[groupIndex][Nlen_id] * sizeof(int));
	memcpy(bm->RBCestimation.RBCspeciesArray[groupIndex].LenComp[1][1][84][2], (int[]
				)
					{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 2, 5, 0, 4, 9, 4, 7, 14, 16, 20, 16, 14, 14, 22, 10, 12, 5, 5, 6, 3, 6, 2, 0, 0, 2, 0,
							0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
							0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
			bm->RBCestimation.RBCspeciesParam[groupIndex][Nlen_id] * sizeof(int));
	memcpy(bm->RBCestimation.RBCspeciesArray[groupIndex].LenComp[2][1][85][0], (int[]
				)
					{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 1, 0, 1, 6, 8, 10, 10, 22, 21, 28, 36, 30, 30, 36, 26, 21, 13, 13, 10, 14, 10,
							11, 6, 6, 3, 1, 3, 3, 3, 0, 1, 1, 0, 2, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
							0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
			bm->RBCestimation.RBCspeciesParam[groupIndex][Nlen_id] * sizeof(int));
	memcpy(bm->RBCestimation.RBCspeciesArray[groupIndex].LenComp[1][1][85][2], (int[]
				)
					{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 2, 0, 2, 1, 2, 8, 8, 12, 17, 7, 9, 9, 10, 16, 9, 6, 7, 4, 2, 0, 0, 1, 1, 1, 0, 0, 0, 1,
							1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
							0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, bm->RBCestimation.RBCspeciesParam[groupIndex][Nlen_id] * sizeof(int));
	memcpy(bm->RBCestimation.RBCspeciesArray[groupIndex].LenComp[2][1][86][0], (int[]
				)
					{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 5, 4, 9, 14, 22, 19, 24, 34, 41, 31, 30, 33, 39, 29, 21, 20, 15, 20, 6, 6,
							11, 9, 2, 6, 0, 0, 1, 1, 1, 0, 1, 2, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
							0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
			bm->RBCestimation.RBCspeciesParam[groupIndex][Nlen_id] * sizeof(int));
	memcpy(bm->RBCestimation.RBCspeciesArray[groupIndex].LenComp[1][1][86][2], (int[]
				)
					{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 1, 4, 2, 9, 7, 17, 9, 18, 21, 15, 12, 10, 7, 10, 3, 6, 1, 2, 0, 0, 1, 0, 0, 0, 0,
							0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
							0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, bm->RBCestimation.RBCspeciesParam[groupIndex][Nlen_id] * sizeof(int));
	memcpy(bm->RBCestimation.RBCspeciesArray[groupIndex].LenComp[2][1][87][0], (int[]
				)
					{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 3, 5, 8, 12, 21, 33, 40, 27, 38, 34, 50, 39, 25, 22, 15, 11, 10, 11, 5,
							2, 1, 2, 2, 4, 3, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
							0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
			bm->RBCestimation.RBCspeciesParam[groupIndex][Nlen_id] * sizeof(int));
	memcpy(bm->RBCestimation.RBCspeciesArray[groupIndex].LenComp[1][1][87][2], (int[]
				)
					{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 1, 1, 7, 6, 4, 12, 14, 17, 20, 17, 12, 8, 11, 13, 3, 2, 4, 2, 0, 0, 0, 0, 0,
							0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
							0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
			bm->RBCestimation.RBCspeciesParam[groupIndex][Nlen_id] * sizeof(int));
	memcpy(bm->RBCestimation.RBCspeciesArray[groupIndex].LenComp[2][1][88][0], (int[]
				)
					{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 1, 1, 2, 4, 8, 12, 10, 12, 20, 28, 31, 33, 46, 35, 32, 16, 13, 14, 8, 8, 6, 6,
							3, 2, 1, 2, 2, 1, 2, 1, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
							0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
			bm->RBCestimation.RBCspeciesParam[groupIndex][Nlen_id] * sizeof(int));
	memcpy(bm->RBCestimation.RBCspeciesArray[groupIndex].LenComp[1][1][88][2], (int[]
				)
					{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 5, 3, 10, 9, 11, 8, 13, 11, 14, 14, 13, 10, 6, 4, 3, 3, 1, 0, 1, 0, 1, 0,
							0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
							0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
			bm->RBCestimation.RBCspeciesParam[groupIndex][Nlen_id] * sizeof(int));
	memcpy(bm->RBCestimation.RBCspeciesArray[groupIndex].LenComp[2][1][89][0], (int[]
				)
					{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 1, 5, 2, 5, 11, 17, 15, 23, 30, 33, 39, 31, 27, 25, 18, 13, 11, 8, 6, 4, 4,
							3, 4, 4, 2, 3, 1, 0, 0, 2, 1, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
							0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
			bm->RBCestimation.RBCspeciesParam[groupIndex][Nlen_id] * sizeof(int));
	memcpy(bm->RBCestimation.RBCspeciesArray[groupIndex].LenComp[2][1][90][0], (int[]
				)
					{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 1, 2, 7, 5, 6, 12, 10, 19, 27, 25, 29, 37, 34, 31, 31, 31, 25, 24, 17, 8, 9, 4,
							9, 8, 3, 2, 5, 1, 1, 2, 4, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
							0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
			bm->RBCestimation.RBCspeciesParam[groupIndex][Nlen_id] * sizeof(int));
	memcpy(bm->RBCestimation.RBCspeciesArray[groupIndex].LenComp[2][1][91][0], (int[]
				)
					{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 4, 5, 5, 8, 14, 21, 31, 32, 26, 28, 41, 38, 31, 30, 35, 16, 25, 4, 4, 5, 6,
							4, 3, 3, 2, 1, 1, 0, 2, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
							0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
			bm->RBCestimation.RBCspeciesParam[groupIndex][Nlen_id] * sizeof(int));
	memcpy(bm->RBCestimation.RBCspeciesArray[groupIndex].LenComp[1][1][91][2], (int[]
				)
					{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 2, 0, 3, 3, 4, 8, 9, 15, 5, 9, 13, 15, 7, 15, 8, 3, 6, 3, 0, 3, 0, 0, 0, 0, 0, 0, 0,
							0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
							0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, bm->RBCestimation.RBCspeciesParam[groupIndex][Nlen_id] * sizeof(int));
	memcpy(bm->RBCestimation.RBCspeciesArray[groupIndex].LenComp[2][1][92][0], (int[]
				)
					{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 4, 9, 8, 14, 20, 25, 33, 41, 45, 39, 47, 35, 37, 28, 19, 14, 13, 6, 8,
							5, 5, 6, 2, 1, 2, 2, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
							0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
			bm->RBCestimation.RBCspeciesParam[groupIndex][Nlen_id] * sizeof(int));
	memcpy(bm->RBCestimation.RBCspeciesArray[groupIndex].LenComp[2][2][80][0], (int[]
				)
					{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 2, 2, 7, 12, 9, 19, 15, 42, 26, 32, 40, 41, 51, 38, 38, 39, 37, 28, 24, 12,
							19, 20, 6, 11, 8, 13, 13, 5, 2, 2, 1, 1, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
							0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
			bm->RBCestimation.RBCspeciesParam[groupIndex][Nlen_id] * sizeof(int));
	memcpy(bm->RBCestimation.RBCspeciesArray[groupIndex].LenComp[2][2][81][0], (int[]
				)
					{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 1, 3, 3, 9, 9, 21, 18, 28, 27, 41, 35, 48, 56, 38, 49, 52, 26, 28, 19, 17,
							18, 26, 8, 10, 4, 7, 5, 7, 3, 3, 2, 3, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
							0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
			bm->RBCestimation.RBCspeciesParam[groupIndex][Nlen_id] * sizeof(int));
	memcpy(bm->RBCestimation.RBCspeciesArray[groupIndex].LenComp[2][2][82][0], (int[]
				)
					{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 1, 0, 1, 0, 2, 2, 10, 7, 18, 14, 24, 24, 27, 31, 45, 35, 43, 37, 39, 31, 31, 33, 19,
							20, 14, 6, 9, 12, 8, 6, 6, 2, 3, 1, 2, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
							0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
			bm->RBCestimation.RBCspeciesParam[groupIndex][Nlen_id] * sizeof(int));
	memcpy(bm->RBCestimation.RBCspeciesArray[groupIndex].LenComp[2][2][83][0], (int[]
				)
					{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 3, 4, 10, 7, 4, 18, 14, 20, 26, 41, 37, 40, 36, 48, 39, 28, 32, 28, 28, 27,
							16, 18, 19, 11, 8, 12, 5, 6, 2, 2, 0, 0, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
							0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
			bm->RBCestimation.RBCspeciesParam[groupIndex][Nlen_id] * sizeof(int));
	memcpy(bm->RBCestimation.RBCspeciesArray[groupIndex].LenComp[1][2][83][2], (int[]
				)
					{ 0, 0, 0, 0, 6, 1, 1, 0, 0, 0, 0, 3, 9, 12, 14, 10, 5, 2, 2, 1, 11, 8, 24, 19, 16, 28, 19, 21, 22, 13, 16, 13, 16, 13, 16, 8, 7, 3, 8, 5,
							1, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
							0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
			bm->RBCestimation.RBCspeciesParam[groupIndex][Nlen_id] * sizeof(int));
	memcpy(bm->RBCestimation.RBCspeciesArray[groupIndex].LenComp[2][2][84][0], (int[]
				)
					{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 1, 1, 2, 0, 0, 1, 1, 1, 7, 3, 11, 11, 11, 5, 18, 15, 24, 23, 32, 30, 33, 40, 23, 17, 21, 17, 17,
							14, 18, 6, 8, 3, 3, 4, 1, 3, 0, 3, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
							0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
			bm->RBCestimation.RBCspeciesParam[groupIndex][Nlen_id] * sizeof(int));
	memcpy(bm->RBCestimation.RBCspeciesArray[groupIndex].LenComp[1][2][84][2], (int[]
				)
					{ 0, 0, 0, 1, 2, 3, 0, 0, 0, 0, 0, 6, 6, 10, 7, 6, 3, 6, 5, 6, 8, 12, 21, 16, 19, 9, 14, 14, 14, 17, 17, 9, 5, 4, 8, 6, 6, 1, 1, 1, 0, 0, 0,
							0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
							0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
			bm->RBCestimation.RBCspeciesParam[groupIndex][Nlen_id] * sizeof(int));
	memcpy(bm->RBCestimation.RBCspeciesArray[groupIndex].LenComp[2][2][85][0], (int[]
				)
					{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 1, 0, 0, 3, 2, 9, 11, 12, 11, 18, 18, 24, 27, 45, 34, 34, 38, 38, 42, 28, 20, 16, 24,
							20, 13, 12, 12, 7, 5, 5, 4, 1, 3, 3, 3, 4, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
							0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
			bm->RBCestimation.RBCspeciesParam[groupIndex][Nlen_id] * sizeof(int));
	memcpy(bm->RBCestimation.RBCspeciesArray[groupIndex].LenComp[1][2][85][2], (int[]
				)
					{ 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 2, 7, 10, 12, 4, 3, 2, 4, 9, 12, 10, 19, 27, 28, 27, 27, 20, 22, 17, 14, 13, 11, 11, 11, 4, 4, 0, 2, 2,
							0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
							0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
			bm->RBCestimation.RBCspeciesParam[groupIndex][Nlen_id] * sizeof(int));
	memcpy(bm->RBCestimation.RBCspeciesArray[groupIndex].LenComp[2][2][86][0], (int[]
				)
					{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 9, 6, 7, 10, 13, 13, 28, 22, 37, 31, 29, 44, 33, 42, 39, 22, 25, 23, 15, 8,
							3, 7, 8, 3, 3, 5, 3, 1, 1, 1, 1, 0, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
							0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
			bm->RBCestimation.RBCspeciesParam[groupIndex][Nlen_id] * sizeof(int));
	memcpy(bm->RBCestimation.RBCspeciesArray[groupIndex].LenComp[1][2][86][2], (int[]
				)
					{ 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 1, 0, 3, 9, 6, 6, 4, 1, 1, 9, 10, 18, 17, 27, 29, 24, 20, 20, 21, 16, 23, 8, 19, 1, 4, 7, 4, 0, 2, 2, 0, 0,
							0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
							0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
			bm->RBCestimation.RBCspeciesParam[groupIndex][Nlen_id] * sizeof(int));
	memcpy(bm->RBCestimation.RBCspeciesArray[groupIndex].LenComp[2][2][87][0], (int[]
				)
					{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 2, 0, 7, 2, 8, 6, 14, 9, 24, 32, 20, 47, 34, 39, 49, 40, 37, 26, 23, 23, 14, 7,
							9, 3, 5, 3, 3, 3, 1, 0, 3, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
							0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
			bm->RBCestimation.RBCspeciesParam[groupIndex][Nlen_id] * sizeof(int));
	memcpy(bm->RBCestimation.RBCspeciesArray[groupIndex].LenComp[1][2][87][2], (int[]
				)
					{ 0, 0, 0, 1, 4, 1, 0, 0, 0, 0, 0, 3, 0, 12, 8, 9, 2, 0, 2, 6, 4, 9, 14, 15, 21, 23, 22, 24, 20, 19, 14, 16, 19, 7, 8, 6, 5, 1, 2, 1, 0, 0,
							0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
							0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
			bm->RBCestimation.RBCspeciesParam[groupIndex][Nlen_id] * sizeof(int));
	memcpy(bm->RBCestimation.RBCspeciesArray[groupIndex].LenComp[2][2][88][0], (int[]
				)
					{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 1, 2, 2, 6, 7, 11, 9, 18, 23, 37, 34, 40, 47, 56, 37, 36, 36, 44, 26, 19, 18, 11,
							13, 5, 6, 8, 1, 2, 4, 2, 2, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
							0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
			bm->RBCestimation.RBCspeciesParam[groupIndex][Nlen_id] * sizeof(int));
	memcpy(bm->RBCestimation.RBCspeciesArray[groupIndex].LenComp[1][2][88][2], (int[]
				)
					{ 0, 0, 0, 1, 7, 3, 1, 0, 0, 0, 0, 0, 3, 9, 7, 4, 5, 3, 2, 8, 5, 10, 16, 19, 19, 13, 24, 21, 20, 20, 19, 19, 14, 7, 12, 6, 2, 1, 1, 1, 2, 0,
							0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
							0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
			bm->RBCestimation.RBCspeciesParam[groupIndex][Nlen_id] * sizeof(int));
	memcpy(bm->RBCestimation.RBCspeciesArray[groupIndex].LenComp[2][2][89][0], (int[]
				)
					{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 3, 2, 3, 4, 5, 10, 9, 16, 17, 23, 26, 34, 45, 47, 31, 32, 32, 34, 30, 29, 18, 13,
							9, 11, 2, 1, 0, 2, 1, 2, 2, 3, 2, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
							0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
			bm->RBCestimation.RBCspeciesParam[groupIndex][Nlen_id] * sizeof(int));
	memcpy(bm->RBCestimation.RBCspeciesArray[groupIndex].LenComp[1][2][89][2], (int[]
				)
					{ 0, 0, 0, 1, 1, 2, 0, 0, 0, 1, 0, 3, 5, 12, 7, 8, 6, 1, 1, 6, 12, 6, 12, 20, 14, 22, 16, 19, 19, 16, 16, 11, 16, 11, 5, 5, 4, 2, 1, 1, 1,
							0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
							0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
			bm->RBCestimation.RBCspeciesParam[groupIndex][Nlen_id] * sizeof(int));
	memcpy(bm->RBCestimation.RBCspeciesArray[groupIndex].LenComp[2][2][90][0], (int[]
				)
					{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 1, 2, 2, 5, 8, 12, 17, 13, 23, 16, 25, 28, 39, 33, 28, 40, 25, 24, 16, 12, 10,
							12, 8, 10, 6, 1, 5, 2, 2, 0, 2, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
							0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
			bm->RBCestimation.RBCspeciesParam[groupIndex][Nlen_id] * sizeof(int));
	memcpy(bm->RBCestimation.RBCspeciesArray[groupIndex].LenComp[1][2][90][2], (int[]
				)
					{ 0, 0, 0, 1, 3, 3, 0, 0, 0, 0, 0, 3, 0, 9, 9, 10, 3, 2, 4, 4, 6, 14, 17, 14, 16, 16, 14, 22, 20, 19, 14, 16, 17, 7, 3, 2, 3, 0, 2, 0, 0, 0,
							0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
							0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
			bm->RBCestimation.RBCspeciesParam[groupIndex][Nlen_id] * sizeof(int));
	memcpy(bm->RBCestimation.RBCspeciesArray[groupIndex].LenComp[2][2][91][0], (int[]
				)
					{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 1, 6, 4, 12, 12, 14, 20, 26, 26, 28, 23, 28, 43, 29, 33, 27, 24, 14, 20, 8,
							14, 3, 4, 4, 7, 3, 2, 4, 0, 1, 3, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
							0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
			bm->RBCestimation.RBCspeciesParam[groupIndex][Nlen_id] * sizeof(int));
	memcpy(bm->RBCestimation.RBCspeciesArray[groupIndex].LenComp[1][2][91][2], (int[]
				)
					{ 0, 0, 0, 2, 3, 1, 1, 0, 0, 0, 0, 3, 3, 10, 3, 6, 4, 1, 6, 8, 15, 14, 26, 18, 18, 28, 22, 15, 15, 11, 16, 13, 3, 6, 2, 4, 1, 1, 3, 0, 0, 0,
							0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
							0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
			bm->RBCestimation.RBCspeciesParam[groupIndex][Nlen_id] * sizeof(int));
	memcpy(bm->RBCestimation.RBCspeciesArray[groupIndex].LenComp[2][2][92][0], (int[]
				)
					{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 3, 3, 0, 0, 0, 0, 2, 4, 4, 5, 9, 10, 22, 19, 25, 29, 27, 42, 40, 36, 34, 31, 26, 21, 12, 12, 14,
							9, 5, 7, 0, 6, 3, 3, 1, 1, 4, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
							0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
			bm->RBCestimation.RBCspeciesParam[groupIndex][Nlen_id] * sizeof(int));
	memcpy(bm->RBCestimation.RBCspeciesArray[groupIndex].LenComp[1][2][92][2], (int[]
				)
					{ 0, 0, 0, 0, 3, 1, 0, 0, 0, 0, 0, 1, 4, 8, 7, 3, 6, 3, 3, 5, 8, 10, 13, 15, 19, 19, 14, 29, 13, 19, 19, 18, 9, 7, 3, 0, 2, 2, 0, 3, 1, 0,
							0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
							0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
			bm->RBCestimation.RBCspeciesParam[groupIndex][Nlen_id] * sizeof(int));
	memcpy(bm->RBCestimation.RBCspeciesArray[groupIndex].LenComp[2][3][80][0], (int[]
				)
					{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 0, 4, 2, 1, 4, 3, 8, 3, 4, 4, 6, 1, 4, 4, 4, 3, 2, 4,
							4, 2, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
							0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, bm->RBCestimation.RBCspeciesParam[groupIndex][Nlen_id] * sizeof(int));
	memcpy(bm->RBCestimation.RBCspeciesArray[groupIndex].LenComp[2][3][81][0], (int[]
				)
					{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 3, 0, 1, 2, 1, 0, 1, 3, 2, 2, 0, 0, 4, 1, 1, 2, 2, 3,
							3, 1, 2, 2, 1, 2, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
							0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, bm->RBCestimation.RBCspeciesParam[groupIndex][Nlen_id] * sizeof(int));
	memcpy(bm->RBCestimation.RBCspeciesArray[groupIndex].LenComp[2][3][82][0], (int[]
				)
					{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 1, 0, 1, 2, 0, 4, 6, 5, 4, 4, 4, 2, 3, 6, 5, 3, 2, 4,
							1, 0, 1, 2, 1, 2, 1, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
							0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, bm->RBCestimation.RBCspeciesParam[groupIndex][Nlen_id] * sizeof(int));
	memcpy(bm->RBCestimation.RBCspeciesArray[groupIndex].LenComp[2][3][83][0], (int[]
				)
					{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 2, 1, 0, 3, 2, 7, 1, 1, 9, 4, 5, 13, 3, 4, 6, 1, 3, 2, 3, 3,
							1, 1, 0, 0, 1, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
							0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, bm->RBCestimation.RBCspeciesParam[groupIndex][Nlen_id] * sizeof(int));
	memcpy(bm->RBCestimation.RBCspeciesArray[groupIndex].LenComp[2][3][84][0], (int[]
				)
					{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 3, 3, 5, 1, 1, 3, 9, 5, 4, 2, 3, 4, 4, 4, 2, 4, 1, 0,
							1, 2, 2, 0, 1, 2, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
							0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, bm->RBCestimation.RBCspeciesParam[groupIndex][Nlen_id] * sizeof(int));
	memcpy(bm->RBCestimation.RBCspeciesArray[groupIndex].LenComp[2][3][85][0], (int[]
				)
					{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 3, 3, 4, 3, 5, 2, 4, 2, 0, 5, 2, 2, 3, 2, 3, 1, 2, 1,
							3, 1, 1, 1, 0, 1, 1, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
							0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, bm->RBCestimation.RBCspeciesParam[groupIndex][Nlen_id] * sizeof(int));
	memcpy(bm->RBCestimation.RBCspeciesArray[groupIndex].LenComp[2][3][86][0], (int[]
				)
					{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 3, 4, 1, 3, 4, 2, 2, 2, 3, 2, 2, 0, 1, 2, 3, 1, 2, 0, 1,
							1, 0, 1, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
							0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, bm->RBCestimation.RBCspeciesParam[groupIndex][Nlen_id] * sizeof(int));
	memcpy(bm->RBCestimation.RBCspeciesArray[groupIndex].LenComp[2][3][87][0], (int[]
				)
					{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 1, 5, 0, 8, 6, 6, 10, 8, 7, 4, 2, 1, 1, 2, 2, 0, 2, 0,
							1, 2, 1, 2, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
							0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, bm->RBCestimation.RBCspeciesParam[groupIndex][Nlen_id] * sizeof(int));
	memcpy(bm->RBCestimation.RBCspeciesArray[groupIndex].LenComp[2][3][88][0], (int[]
				)
					{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 3, 1, 1, 4, 2, 8, 6, 3, 3, 6, 4, 3, 2, 0, 6, 1, 2, 1, 0,
							3, 1, 3, 1, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
							0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, bm->RBCestimation.RBCspeciesParam[groupIndex][Nlen_id] * sizeof(int));
	memcpy(bm->RBCestimation.RBCspeciesArray[groupIndex].LenComp[2][3][89][0], (int[]
				)
					{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 4, 5, 5, 9, 7, 9, 11, 10, 15, 11, 10, 5, 6, 4, 6, 5, 5,
							3, 1, 4, 3, 2, 1, 0, 0, 0, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
							0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, bm->RBCestimation.RBCspeciesParam[groupIndex][Nlen_id] * sizeof(int));
	memcpy(bm->RBCestimation.RBCspeciesArray[groupIndex].LenComp[2][3][90][0], (int[]
				)
					{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 3, 10, 5, 8, 7, 10, 6, 10, 8, 7, 6, 14, 11, 3, 4, 2, 4,
							2, 2, 3, 4, 0, 3, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
							0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, bm->RBCestimation.RBCspeciesParam[groupIndex][Nlen_id] * sizeof(int));
	memcpy(bm->RBCestimation.RBCspeciesArray[groupIndex].LenComp[2][3][91][0], (int[]
				)
					{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 4, 4, 4, 7, 9, 6, 10, 13, 11, 6, 12, 2, 5, 5, 9, 4, 2, 1, 2,
							1, 1, 1, 1, 1, 2, 0, 1, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
							0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, bm->RBCestimation.RBCspeciesParam[groupIndex][Nlen_id] * sizeof(int));
	memcpy(bm->RBCestimation.RBCspeciesArray[groupIndex].LenComp[2][3][92][0], (int[]
				)
					{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 2, 5, 0, 5, 1, 5, 5, 3, 3, 3, 3, 2, 0, 5, 1, 2, 3, 0,
							1, 2, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
							0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, bm->RBCestimation.RBCspeciesParam[groupIndex][Nlen_id] * sizeof(int));

	bm->RBCestimation.RBCspeciesParam[groupIndex][MaxAge_id] = 21;

	AccumAge = bm->RBCestimation.RBCspeciesParam[groupIndex][MaxAge_id] + ceil(bm->RBCestimation.RBCspeciesParam[groupIndex][MaxAge_id] / 4.0); // accumulator age used in SS

	bm->RBCestimation.RBCspeciesParam[groupIndex][AccumAge_id] = AccumAge;
	printf("AccumAge = %d\n", AccumAge);

	bm->RBCestimation.RBCspeciesArray[groupIndex].Ageing_error = Util_Alloc_Init_1D_Double(AccumAge, 0.0);
	bm->RBCestimation.RBCspeciesArray[groupIndex].Ageing_error[0] = 0.1965;
	bm->RBCestimation.RBCspeciesArray[groupIndex].Ageing_error[1] = 0.2766;
	bm->RBCestimation.RBCspeciesArray[groupIndex].Ageing_error[2] = 0.3589;
	bm->RBCestimation.RBCspeciesArray[groupIndex].Ageing_error[3] = 0.4435;
	bm->RBCestimation.RBCspeciesArray[groupIndex].Ageing_error[4] = 0.5304;
	bm->RBCestimation.RBCspeciesArray[groupIndex].Ageing_error[5] = 0.6198;
	bm->RBCestimation.RBCspeciesArray[groupIndex].Ageing_error[6] = 0.7116;
	bm->RBCestimation.RBCspeciesArray[groupIndex].Ageing_error[7] = 0.8059;
	bm->RBCestimation.RBCspeciesArray[groupIndex].Ageing_error[8] = 0.9029;
	bm->RBCestimation.RBCspeciesArray[groupIndex].Ageing_error[9] = 1.0025;
	bm->RBCestimation.RBCspeciesArray[groupIndex].Ageing_error[10] = 1.1049;
	bm->RBCestimation.RBCspeciesArray[groupIndex].Ageing_error[11] = 1.2102;
	bm->RBCestimation.RBCspeciesArray[groupIndex].Ageing_error[12] = 1.3183;
	bm->RBCestimation.RBCspeciesArray[groupIndex].Ageing_error[13] = 1.4294;
	bm->RBCestimation.RBCspeciesArray[groupIndex].Ageing_error[14] = 1.5436;
	bm->RBCestimation.RBCspeciesArray[groupIndex].Ageing_error[15] = 1.661;
	bm->RBCestimation.RBCspeciesArray[groupIndex].Ageing_error[16] = 1.7816;
	bm->RBCestimation.RBCspeciesArray[groupIndex].Ageing_error[17] = 1.9056;
	bm->RBCestimation.RBCspeciesArray[groupIndex].Ageing_error[18] = 2.0329;
	bm->RBCestimation.RBCspeciesArray[groupIndex].Ageing_error[19] = 2.1638;
	bm->RBCestimation.RBCspeciesArray[groupIndex].Ageing_error[20] = 2.2983;
	bm->RBCestimation.RBCspeciesArray[groupIndex].Ageing_error[21] = 2.29;
	bm->RBCestimation.RBCspeciesArray[groupIndex].Ageing_error[22] = 2.29;
	bm->RBCestimation.RBCspeciesArray[groupIndex].Ageing_error[23] = 2.29;
	bm->RBCestimation.RBCspeciesArray[groupIndex].Ageing_error[24] = 2.29;
	bm->RBCestimation.RBCspeciesArray[groupIndex].Ageing_error[25] = 2.29;

	//bm->RBCestimation.RBCspeciesArray[groupIndex].AFss = Util_Alloc_Init_4D_Int(4, nYears, bm->K_num_sexes,
	//		bm->RBCestimation.RBCspeciesParam[groupIndex][NumFisheries_id], 0);

	//	bm->RBCestimation.RBCspeciesArray[groupIndex].LFss[1][0][80][0] = 309;

	bm->RBCestimation.RBCspeciesArray[groupIndex].AFss[1][0][83][0] = 176;
	bm->RBCestimation.RBCspeciesArray[groupIndex].AFss[1][0][85][0] = 188;
	bm->RBCestimation.RBCspeciesArray[groupIndex].AFss[1][0][86][0] = 223;
	bm->RBCestimation.RBCspeciesArray[groupIndex].AFss[1][0][87][0] = 200;
	bm->RBCestimation.RBCspeciesArray[groupIndex].AFss[1][0][88][0] = 175;
	bm->RBCestimation.RBCspeciesArray[groupIndex].AFss[1][0][89][0] = 197;
	bm->RBCestimation.RBCspeciesArray[groupIndex].AFss[1][0][90][0] = 237;
	bm->RBCestimation.RBCspeciesArray[groupIndex].AFss[1][0][91][0] = 222;
	bm->RBCestimation.RBCspeciesArray[groupIndex].AFss[1][0][92][0] = 236;
	bm->RBCestimation.RBCspeciesArray[groupIndex].AFss[2][0][83][0] = 297;
	bm->RBCestimation.RBCspeciesArray[groupIndex].AFss[2][0][84][0] = 229;
	bm->RBCestimation.RBCspeciesArray[groupIndex].AFss[2][0][85][0] = 290;
	bm->RBCestimation.RBCspeciesArray[groupIndex].AFss[2][0][86][0] = 256;
	bm->RBCestimation.RBCspeciesArray[groupIndex].AFss[2][0][87][0] = 259;
	bm->RBCestimation.RBCspeciesArray[groupIndex].AFss[2][0][88][0] = 290;
	bm->RBCestimation.RBCspeciesArray[groupIndex].AFss[2][0][89][0] = 236;
	bm->RBCestimation.RBCspeciesArray[groupIndex].AFss[2][0][90][0] = 203;
	bm->RBCestimation.RBCspeciesArray[groupIndex].AFss[2][0][91][0] = 226;
	bm->RBCestimation.RBCspeciesArray[groupIndex].AFss[2][0][92][0] = 232;
	bm->RBCestimation.RBCspeciesArray[groupIndex].AFss[3][0][84][0] = 39;
	bm->RBCestimation.RBCspeciesArray[groupIndex].AFss[3][0][85][0] = 22;
	bm->RBCestimation.RBCspeciesArray[groupIndex].AFss[3][0][87][0] = 41;
	bm->RBCestimation.RBCspeciesArray[groupIndex].AFss[3][0][90][0] = 60;

	bm->RBCestimation.RBCspeciesArray[groupIndex].AgeComp = Util_Alloc_Init_5D_Int(bm->RBCestimation.RBCspeciesParam[groupIndex][MaxAge_id], 3, nYears,
			bm->K_num_sexes, bm->K_num_fisheries, 0);

	memcpy(bm->RBCestimation.RBCspeciesArray[groupIndex].AgeComp[1][0][183][0], (int[]
				)
					{ 0, 0, 20, 31, 32, 24, 21, 18, 17, 1, 4, 4, 2, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
			bm->RBCestimation.RBCspeciesParam[groupIndex][MaxAge_id] * sizeof(int));
	memcpy(bm->RBCestimation.RBCspeciesArray[groupIndex].AgeComp[1][0][185][0], (int[]
				)
					{ 0, 0, 22, 51, 31, 23, 17, 10, 20, 4, 4, 0, 1, 0, 1, 2, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
			bm->RBCestimation.RBCspeciesParam[groupIndex][MaxAge_id] * sizeof(int));
	memcpy(bm->RBCestimation.RBCspeciesArray[groupIndex].AgeComp[1][0][186][0], (int[]
				)
					{ 0, 0, 34, 75, 47, 19, 15, 10, 7, 7, 2, 2, 0, 2, 0, 0, 0, 2, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
			bm->RBCestimation.RBCspeciesParam[groupIndex][MaxAge_id] * sizeof(int));
	memcpy(bm->RBCestimation.RBCspeciesArray[groupIndex].AgeComp[1][0][187][0], (int[]
				)
					{ 0, 0, 15, 69, 52, 24, 7, 13, 5, 3, 8, 1, 2, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
			bm->RBCestimation.RBCspeciesParam[groupIndex][MaxAge_id] * sizeof(int));
	memcpy(bm->RBCestimation.RBCspeciesArray[groupIndex].AgeComp[1][0][188][0], (int[]
				)
					{ 0, 0, 14, 40, 53, 32, 24, 6, 1, 2, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
			bm->RBCestimation.RBCspeciesParam[groupIndex][MaxAge_id] * sizeof(int));
	memcpy(bm->RBCestimation.RBCspeciesArray[groupIndex].AgeComp[1][0][189][0], (int[]
				)
					{ 0, 0, 20, 47, 32, 41, 29, 11, 8, 2, 0, 1, 2, 1, 2, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
			bm->RBCestimation.RBCspeciesParam[groupIndex][MaxAge_id] * sizeof(int));
	memcpy(bm->RBCestimation.RBCspeciesArray[groupIndex].AgeComp[1][0][190][0], (int[]
				)
					{ 0, 0, 30, 73, 44, 23, 32, 16, 8, 1, 3, 3, 2, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
			bm->RBCestimation.RBCspeciesParam[groupIndex][MaxAge_id] * sizeof(int));
	memcpy(bm->RBCestimation.RBCspeciesArray[groupIndex].AgeComp[1][0][191][0], (int[]
				)
					{ 0, 0, 22, 74, 40, 31, 11, 23, 5, 3, 5, 2, 1, 2, 0, 1, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
			bm->RBCestimation.RBCspeciesParam[groupIndex][MaxAge_id] * sizeof(int));
	memcpy(bm->RBCestimation.RBCspeciesArray[groupIndex].AgeComp[1][0][192][0], (int[]
				)
					{ 0, 0, 23, 88, 59, 27, 13, 6, 11, 6, 2, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
			bm->RBCestimation.RBCspeciesParam[groupIndex][MaxAge_id] * sizeof(int));
	memcpy(bm->RBCestimation.RBCspeciesArray[groupIndex].AgeComp[2][0][183][0], (int[]
				)
					{ 0, 0, 18, 39, 56, 41, 46, 25, 31, 3, 3, 7, 5, 10, 4, 6, 0, 1, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
			bm->RBCestimation.RBCspeciesParam[groupIndex][MaxAge_id] * sizeof(int));
	memcpy(bm->RBCestimation.RBCspeciesArray[groupIndex].AgeComp[2][0][184][0], (int[]
				)
					{ 0, 0, 19, 44, 35, 30, 23, 30, 12, 12, 4, 5, 4, 5, 5, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
			bm->RBCestimation.RBCspeciesParam[groupIndex][MaxAge_id] * sizeof(int));
	memcpy(bm->RBCestimation.RBCspeciesArray[groupIndex].AgeComp[2][0][185][0], (int[]
				)
					{ 0, 1, 31, 81, 44, 23, 27, 22, 23, 7, 6, 4, 2, 2, 6, 5, 1, 2, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
			bm->RBCestimation.RBCspeciesParam[groupIndex][MaxAge_id] * sizeof(int));
	memcpy(bm->RBCestimation.RBCspeciesArray[groupIndex].AgeComp[2][0][186][0], (int[]
				)
					{ 0, 0, 35, 81, 52, 33, 11, 10, 11, 6, 5, 4, 0, 0, 1, 1, 3, 1, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
			bm->RBCestimation.RBCspeciesParam[groupIndex][MaxAge_id] * sizeof(int));
	memcpy(bm->RBCestimation.RBCspeciesArray[groupIndex].AgeComp[2][0][187][0], (int[]
				)
					{ 0, 0, 12, 78, 66, 39, 19, 9, 8, 5, 9, 4, 5, 0, 0, 1, 1, 2, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
			bm->RBCestimation.RBCspeciesParam[groupIndex][MaxAge_id] * sizeof(int));
	memcpy(bm->RBCestimation.RBCspeciesArray[groupIndex].AgeComp[2][0][188][0], (int[]
				)
					{ 0, 0, 23, 50, 85, 49, 35, 14, 6, 6, 5, 5, 1, 1, 1, 2, 1, 2, 1, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
			bm->RBCestimation.RBCspeciesParam[groupIndex][MaxAge_id] * sizeof(int));
	memcpy(bm->RBCestimation.RBCspeciesArray[groupIndex].AgeComp[2][0][189][0], (int[]
				)
					{ 0, 1, 17, 52, 38, 65, 29, 13, 6, 2, 1, 3, 2, 2, 1, 1, 0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
			bm->RBCestimation.RBCspeciesParam[groupIndex][MaxAge_id] * sizeof(int));
	memcpy(bm->RBCestimation.RBCspeciesArray[groupIndex].AgeComp[2][0][190][0], (int[]
				)
					{ 0, 3, 20, 44, 54, 21, 33, 8, 6, 3, 3, 2, 2, 1, 0, 0, 0, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
			bm->RBCestimation.RBCspeciesParam[groupIndex][MaxAge_id] * sizeof(int));
	memcpy(bm->RBCestimation.RBCspeciesArray[groupIndex].AgeComp[2][0][191][0], (int[]
				)
					{ 0, 2, 26, 68, 51, 27, 16, 15, 8, 3, 1, 6, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
			bm->RBCestimation.RBCspeciesParam[groupIndex][MaxAge_id] * sizeof(int));
	memcpy(bm->RBCestimation.RBCspeciesArray[groupIndex].AgeComp[2][0][192][0], (int[]
				)
					{ 0, 0, 17, 76, 53, 27, 13, 8, 11, 10, 7, 4, 1, 2, 1, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
			bm->RBCestimation.RBCspeciesParam[groupIndex][MaxAge_id] * sizeof(int));
	memcpy(bm->RBCestimation.RBCspeciesArray[groupIndex].AgeComp[3][0][184][0], (int[]
				)
					{ 0, 0, 1, 2, 3, 10, 5, 7, 6, 2, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
			bm->RBCestimation.RBCspeciesParam[groupIndex][MaxAge_id] * sizeof(int));
	memcpy(bm->RBCestimation.RBCspeciesArray[groupIndex].AgeComp[3][0][185][0], (int[]
				)
					{ 0, 0, 0, 6, 4, 2, 3, 2, 1, 1, 1, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
			bm->RBCestimation.RBCspeciesParam[groupIndex][MaxAge_id] * sizeof(int));
	memcpy(bm->RBCestimation.RBCspeciesArray[groupIndex].AgeComp[3][0][187][0], (int[]
				)
					{ 0, 0, 1, 9, 10, 7, 1, 3, 3, 2, 1, 0, 1, 0, 0, 0, 1, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
			bm->RBCestimation.RBCspeciesParam[groupIndex][MaxAge_id] * sizeof(int));

	memcpy(bm->RBCestimation.RBCspeciesArray[groupIndex].AgeComp[3][0][190][0], (int[]
				)
					{ 0, 0, 3, 15, 6, 8, 8, 10, 3, 1, 2, 1, 1, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
			bm->RBCestimation.RBCspeciesParam[groupIndex][MaxAge_id] * sizeof(int));

	bm->RBCestimation.RBCspeciesArray[groupIndex].AgeFltYr = Util_Alloc_Init_3D_Int(nYears, bm->K_num_fisheries, 4, 0.0);

	memcpy(bm->RBCestimation.RBCspeciesArray[groupIndex].AgeFltYr[0][0], (int[]
				)
					{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
							0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
							0, 0, 0, 0 }, nYears * sizeof(int));

	memcpy(bm->RBCestimation.RBCspeciesArray[groupIndex].AgeFltYr[0][1], (int[]
				)
					{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
							0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 1, 1, 1, 1,
							1, 1, 1, 1 }, nYears * sizeof(int));

	memcpy(bm->RBCestimation.RBCspeciesArray[groupIndex].AgeFltYr[0][1], (int[]
				)
					{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
							0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1,
							1, 1, 1, 1 }, nYears * sizeof(int));

	memcpy(bm->RBCestimation.RBCspeciesArray[groupIndex].AgeFltYr[0][1], (int[]
				)
					{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
							0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 1, 0, 0,
							1, 0, 0, 1 }, nYears * sizeof(int));

	bm->RBCestimation.RBCspeciesArray[groupIndex].SSMort[0][0] = 0.22;
	bm->RBCestimation.RBCspeciesArray[groupIndex].SSMort[0][1] = 0.22;

	bm->RBCestimation.RBCspeciesArray[groupIndex].Start_SelInflect[0] = -31.8673;
	bm->RBCestimation.RBCspeciesArray[groupIndex].Start_SelInflect[1] = 32.4;
	bm->RBCestimation.RBCspeciesArray[groupIndex].Start_SelInflect[2] = 36.5;
	bm->RBCestimation.RBCspeciesArray[groupIndex].Start_SelInflect[3] = 52.3;

	bm->RBCestimation.RBCspeciesArray[groupIndex].Start_SelWidth[0] = -6.50136;
	bm->RBCestimation.RBCspeciesArray[groupIndex].Start_SelWidth[1] = 5.99;
	bm->RBCestimation.RBCspeciesArray[groupIndex].Start_SelWidth[2] = 10.88;
	bm->RBCestimation.RBCspeciesArray[groupIndex].Start_SelWidth[3] = -15.9286;

	bm->RBCestimation.RBCspeciesArray[groupIndex].Start_RetInflect[0] = -25.0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].Start_RetInflect[1] = 16.98;
	bm->RBCestimation.RBCspeciesArray[groupIndex].Start_RetInflect[2] = 22.85;
	bm->RBCestimation.RBCspeciesArray[groupIndex].Start_RetInflect[3] = 31.5;

	//-1.0 5.37 4.44 0.81
	bm->RBCestimation.RBCspeciesArray[groupIndex].Start_RetSlope[0] = -1.0;
	bm->RBCestimation.RBCspeciesArray[groupIndex].Start_RetSlope[1] = 5.37;
	bm->RBCestimation.RBCspeciesArray[groupIndex].Start_RetSlope[2] = 4.44;
	bm->RBCestimation.RBCspeciesArray[groupIndex].Start_RetSlope[2] = 0.81;

	bm->RBCestimation.RBCspeciesArray[groupIndex].Hsteep[0] = 0.672796;

	bm->RBCestimation.RBCspeciesArray[groupIndex].Hsteep[0] = 0.672796;

	bm->RBCestimation.RBCspeciesArray[groupIndex].Varadj_CPUE[0] = 0.2;
	bm->RBCestimation.RBCspeciesArray[groupIndex].Varadj_CPUE[1] = 0.179;
	bm->RBCestimation.RBCspeciesArray[groupIndex].Varadj_CPUE[2] = 0.14;
	bm->RBCestimation.RBCspeciesArray[groupIndex].Varadj_CPUE[3] = 0.23;

	bm->RBCestimation.RBCspeciesArray[groupIndex].Varadj_discard[0] = 0.5;
	bm->RBCestimation.RBCspeciesArray[groupIndex].Varadj_discard[1] = 0.02;
	bm->RBCestimation.RBCspeciesArray[groupIndex].Varadj_discard[2] = 0.5;
	bm->RBCestimation.RBCspeciesArray[groupIndex].Varadj_discard[3] = 0.5;

	bm->RBCestimation.RBCspeciesArray[groupIndex].Varadj_length[0] = 1.7;
	bm->RBCestimation.RBCspeciesArray[groupIndex].Varadj_length[1] = 1.2;
	bm->RBCestimation.RBCspeciesArray[groupIndex].Varadj_length[2] = 2.14;
	bm->RBCestimation.RBCspeciesArray[groupIndex].Varadj_length[3] = 0.6;

	bm->RBCestimation.RBCspeciesArray[groupIndex].Varadj_age[0] = 1;
	bm->RBCestimation.RBCspeciesArray[groupIndex].Varadj_age[1] = 0.51;
	bm->RBCestimation.RBCspeciesArray[groupIndex].Varadj_age[2] = 0.93;
	bm->RBCestimation.RBCspeciesArray[groupIndex].Varadj_age[3] = 0.8;

	bm->RBCestimation.RBCspeciesArray[groupIndex].WtLen = Util_Alloc_Init_3D_Double((int) bm->RBCestimation.RBCspeciesParam[groupIndex][Nlen_id], 2, 2, 0.0);

	bm->RBCestimation.DiscType = 2;

	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthN[0] = 1000;
	bm->RBCestimation.RBCspeciesArray[groupIndex].LengthN[1] = 500;

	bm->RBCestimation.RBCspeciesArray[groupIndex].LenRetain = Util_Alloc_Init_3D_Double(nYears, (int) bm->RBCestimation.RBCspeciesParam[groupIndex][Nlen_id],
			(int) bm->RBCestimation.RBCspeciesParam[groupIndex][NumFisheries_id], 0.0);


	bm->RBCestimation.RBCspeciesArray[groupIndex].SelLen = Util_Alloc_Init_3D_Double(nYears, (int) bm->RBCestimation.RBCspeciesParam[groupIndex][Nlen_id],
			(int) bm->RBCestimation.RBCspeciesParam[groupIndex][NumFisheries_id], 0.0);
#endif
	return 0;
}

