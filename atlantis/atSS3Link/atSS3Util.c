
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
#include <atSS3Util.h>


//******************************************************************************
//
//  Name: GenMnomial
//  Description: generate from a multinomial
//  Parameters :
//       prop  is the expected proportions (dimension (1,ndim))
//       ndim  is the dimension of vectors of proportions  (Gavin uses fixed value of 2000)
//       tc    is given sample size
//       c     is the generated multinomial sample (dimension (1,ndim))
//       iseed is random number seed
//       withrep = true (1) if sample without replacement else false (0)
//
//  Calls :
//  Called by:
//  Created:  translated from Gavin's Fortran by Sally 13/11/07
//
// ******************************************************************************
void GenMnomial(double *prop, int ndim, int tc, double *c, int *iseed, int withrep)
{
	int i,j,icount,np,a,a1;
	double tot,ntot;
	int *dum = Util_Alloc_Init_1D_Int(15000, 0);

 	// check if some catches were made
	if (tc==0) return;

	Util_Init_1D_Double(c, ndim, 0.0);

	// first total the probabilities and then normalise
	tot = 0.0;
	for (i=0;i<ndim;i++)
		tot = tot + prop[i];

	for (i=0;i<ndim;i++)
		prop[i] = prop[i]/tot;

	// check whether this is stochastic
	// don't put this code because not used ?

	// adjust to handle 'without-replacement sampling'
	ntot = 10000.0;
	if (withrep && tot <14000) ntot = tot;

	// do special set up for index matrix
	np = 0;
	for (i=0;i<ndim;i++){
		icount = (int) (prop[i]*ntot+0.5);
		if (icount>0){
			for (j=0;j<icount;j++)
				dum[np+j] = i;
			np = np + icount;
		}
	}

	// do actual generation
	for (i=0;i<tc;i++){
        a = 0;
		do{
			a1 = (int) (ran3(iseed)*(float)np+1);
			a = dum[a1];
		} while (a<=0 && withrep);
		if (a>0) c[a] += 1;   // add if (a>0) -  a can be 0 if proportion=0 8.Apr.08
		if (withrep) dum[a1] = -1;
	}

	free(dum);
}


//******************************************************************************
//
//  Name: INVM
//  Description: Computes the inverse of a matrix
//
//  Created: got from Natalie's file numfns.cpp
//
// ******************************************************************************


void INVM(double **A, double **AI, int n){

    int i,j;
    double D=1;
    int *INDEX;
	double *col;

	INDEX = Util_Alloc_Init_1D_Int(n, 0);
	col = Util_Alloc_Init_1D_Double(n, 0.0);

	for( i=0 ; i < n; i++)
		INDEX[i] = 1;

    // Perform an LU-decomposition
    LUDCMP(A, n, INDEX, &D);

    // Find the inverse
    for( j=0; j < n; j++){
        for( i=0; i < n; i++) {
        	col[i] = 0.0;
        }
		col[j] = 1.0;
		LUBKSB(A, n, INDEX, col);
		for( i=0; i < n; i++)
			AI[i][j] = col[i];
    }

    // Clean up
    i_free1d(INDEX);
    free1d(col);

    return;
}

//******************************************************************************
//
//  Name: MULTMV
//  Description: Multiplies matrix AA (nxm) and vector B (m) to give vector C (n)
//
//  Created: got from Natalie's file numfns.cpp
//
// ******************************************************************************
//

void MULTMV(double **AA, double *B, double *C, int n, int m){

    int i,j;
    for( i=0; i < n; i++) {
        C[i] = 0.0;
        for( j=0; j < m; j++) {
            C[i] += AA[i][j] * B[j];
        }
	}

    return;
}

//******************************************************************************
//
//  Name: LUDCMP
//  Description: LU decomposition for solving linear equations
//
//  Created: got from Natalie's file numfns.cpp
//
// ******************************************************************************
void LUDCMP (double **AA, int n, int *INDX, double *d){
    int i, j, k, imax = 0;
    double tiny = 1.0e-20;
    double big, dum, sum, temp;
    double *VV;

    VV = Util_Alloc_Init_1D_Double(n, 0.0);

    *d = 1.0;
    for( i=0; i < n; i++) {
        big = 0.;
        for( j=0; j < n; j++) {
            if ((temp=fabs(AA[i][j])) > big)
            	big = temp;
        }
        if(big == 0.0) {
            printf("Singular matrix in routine LUDCMP\n");
        }
        VV[i] = 1.0 / big;
    }

    for( j=0; j < n; j++) {
        for( i=0; i < j; i++) {
			sum = AA[i][j];
 			for( k=0; k < i; k++)
 				sum -= AA[i][k] * AA[k][j];
            AA[i][j] = sum;
        }
        big = 0;
        for( i=j; i < n; i++) {
            sum = AA[i][j];
            for( k=0; k < j; k++)
            	sum -= AA[i][k] * AA[k][j];
            AA[i][j] = sum;
            if ((dum = VV[i] * fabs(sum)) >= big) {
                imax = i;
                big = dum;
            }
        }
        if (j != imax) {
            for (k=0; k < n; k++) {
                dum = AA[imax][k];
                AA[imax][k] = AA[j][k];
                AA[j][k] = dum;
            }
            *d = -(*d);
            VV[imax] = VV[j];
        }
        INDX[j] = imax;
		if (AA[j][j] == 0.0)
			AA[j][j] = tiny;
        if (j != n) {
            dum = 1.0 / AA[j][j];
            for (i=j+1; i < n; i++) {
            	AA[i][j] *= dum;
            }
        }
    }
}

//******************************************************************************
//
//  Name: LUBKSB
//  Description: LU back-substitution for solving linear equations
//
//  Created: got from Natalie's file numfns.cpp
//
// ******************************************************************************
void LUBKSB(double **AA,int n, int *INDX, double *B) {
    int i,ii=0,ip,j;
    double sum;

    for ( i=0; i < n; i++) {
        ip = INDX[i];
        sum = B[ip];
        B[ip] = B[i];
        if (ii) {
            for( j=ii; j < (i-1); j++) {
				sum -= AA[i][j] * B[j];
            }
        } else if (sum) {
			ii=i;
        }
        B[i] = sum;
    }

    for (i = (n - 1); i>=0; i--) {
        sum = B[i];
		for( j= (i+1); j < n; j++)
			sum -= AA[i][j] * B[j];
        B[i] = sum / AA[i][i];
    }
}
