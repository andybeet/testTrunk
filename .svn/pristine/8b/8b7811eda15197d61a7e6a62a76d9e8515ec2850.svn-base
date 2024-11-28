/*********************************************************************

    File:           atreadBM.c

    Created:        Mon Mar 06 18:22:48 EST 1995

    Author:         Stephen Walker
                    CSIRO Division of Oceanography

    Purpose:        Subroutine to read all the data associated
                    with a box model output dump (geometry, variables
                    etc)

    Arguments:      See below

    Returns:        void

    Revisions:      8/8/2004 EA Fulton
					Ported across the code from the southeast (sephys) model

					17/11/2004 EA Fulton
					Converted original routine definitions from
					void
					routine_name(blah,blahblah)
					int blah;
					double blahblah;

					to

					void routine_name(int blah, double blahblah)

					22-May-2009 Bec Gorton
					Changed the ncopen calls to use sjw_ncopen. This function checks that the
					netcdf file exits before calling ncopen which will crash if the file is not
					found.
*********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <signal.h>
#include <math.h>
#include <string.h>
#include <sjwlib.h>
#include <netcdf.h>
#include <atlantisboxmodel.h>
#include <atUtilLib.h>

void readMSEBoxModel(char *name, long dump, MSEBoxModel *bm, FILE *llogfp)
/* name           name of input netCDF file
   dump;          dump number in netCDF file
*/
{
    int fid = 0;
    int b = 0;

    /* Clear the box model - this assumes that it doesn't
     * already contain valid data. If so, allocated memory should
     * be freed first by whoever is calling this routine.
     */
    memset(bm, 0, sizeof(MSEBoxModel));

    /* Store input file name and dump number */
    strcpy(bm->ncIfname,name);
    bm->ncIfdump = dump;

    /* Open input netCDF file */
    fid =  Util_ncopen(bm->inputFolder, name,NC_NOWRITE);

    /* Get global attribute strings */
    ncattget(fid,NC_GLOBAL,"title",bm->title);
    ncattget(fid,NC_GLOBAL,"geometry",bm->geomIfname);
    ncattget(fid,NC_GLOBAL,"parameters",bm->params);

    /* Read the model geometry */
    readMSEBoxModelGeom(bm->geomIfname,bm);

    /* Read the model variable info */
    readBMphysInfo(fid,bm);
    readBMTracerInfo(fid, name, bm);
    readBMEpiInfo(fid,bm);

    /* Allocate memory for variables. Note that the casts (double ***)
     * are only needed to supress compiler warnings. Note also that the
     * allocation routines used here have parameter ordering opposite
     * to the way that the arrays are accessed. Thus if we want a
     * 3d array X, accessed X[k][j][i] (i varies fastest), then we say
     * X = ...alloc3d(ni,nj,nk), where ni is the number of elements
     * in the i dimension etc.
     */
    bm->nom_dz = (double **)alloc2d(bm->wcnz+bm->sednz,bm->nbox);
    bm->dz = (double **)alloc2d(bm->wcnz+bm->sednz,bm->nbox);
    bm->wctr = (double ***)alloc3d(bm->ntracer,bm->wcnz,bm->nbox);
    bm->sedtr = (double ***)alloc3d(bm->ntracer,bm->sednz,bm->nbox);
    bm->epi = (double **)alloc2d(bm->nepi,bm->nbox);
    bm->vol = (double **)alloc2d(bm->wcnz+bm->sednz,bm->nbox);

    /* Set up pointers to data from each box */
    for(b=0; b<bm->nbox; b++) {
		bm->boxes[b].tr = bm->wctr[b];
		bm->boxes[b].sm.tr = bm->sedtr[b];
		bm->boxes[b].epi = bm->epi[b];
		bm->boxes[b].volume = bm->vol[b];
		bm->boxes[b].sm.volume = bm->vol[b]+bm->wcnz;
		bm->boxes[b].dz = bm->dz[b];
		bm->boxes[b].sm.dz = bm->dz[b]+bm->wcnz;

		if(bm->terrestrial_on && bm->boxes[b].type == LAND){
			bm->boxes[b].land.tr = bm->landtr[b];
		}
    }

    /* Read data */
    readBMphysData(fid,bm->ncIfdump,bm, llogfp);
    readBMTracerData(fid,bm->ncIfdump,bm);
    readBMEpiData(fid,bm->ncIfdump,bm);

    /* Close input netCDF file */
    ncclose(fid);
}




/**** Doesn't actually free the memory for the MSEBoxModel structure
 **** itself, but it does free all memory allocated for it's arrays.
 ****/
void freeMSEBoxModel(MSEBoxModel *bm)
{
    i_free1d((int*)bm->is_boundary);
	free2d((double**)bm->nom_dz);
	free2d((double**)bm->dz);
	free3d((double***)bm->wctr);
	free3d((double***)bm->sedtr);
	free2d((double**)bm->epi);
	free2d((double**)bm->vol);

    /* Free memory allocated for model variable info */
    freeBMEpiInfo(bm);
    freeBMTracerInfo(bm);
    freeBMphysInfo(bm);

    /* Free memory allocated for the geometry. */
    freeMSEBoxModelGeom(bm);
}
