/*********************************************************************

 File:           atlandIO.c

 Created:        15-12-2011

 Author:         Beth Fulton
				 CSIRO Marine and Atmospheric Research

 Purpose:        Routines which deal with land flora/fauna data from tracer files. These
				 files are in netCDF format and should be read, written
				 and created using the routines below.

 Arguments:      See each routine below

 Returns:        See each routine below

 Revisions:      15/12/2011 Beth Fulton
				 Ported across the code from the atepiIO.c

 *********************************************************************/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <sjwlib.h>
#include <netcdf.h>
#include <atlantisboxmodel.h>
#include <atUtilLib.h>

int Find_WC_Tracer(MSEBoxModel *bm, char *tracerName){

	int i, index;
	char *name;

	name = strstr(tracerName, "Land_");

	if(name != NULL){
		name = tracerName + strlen("Land_");
	}


	for(i = 0; i < bm->ntracer; i++){

		if( strcmp(name, bm->tinfo[i].name) == 0 ){
		  index = bm->tinfo[i].n;
		  return index;
		}
	}

	return -1;
}

/*******************************************************************//**
 Routine to read the model land data from a netCDF file.
 This routine assumes the structured file as in atepiIO.c

 Routine base code copied from atlIO.c - be careful of copy and paste errors

 *********************************************************************/
void readBMLandInfo(int fid, MSEBoxModel *bm) {
	int ndims = 0; /* Number of dimensions */
	int nvars = 0; /* Number of variables */
	int ngatts = 0; /* Number of global attributes */
	int recdim = 0; /* ID of unlimited dimension */
	int t_id = 0; /* Dimension ids */
	int b_id = 0;
	int dims[MAX_NC_DIMS];
	int *idlist;
	long n = 0;
	int i = 0;

	/* Set netCDF library error handling */
	ncopts = NC_VERBOSE | NC_FATAL;

	if (verbose > 1)
		fprintf(stderr, "Entering readBMLandInfo\n");

	/* Inquire about this file */
	ncinquire(fid, &ndims, &nvars, &ngatts, &recdim);
	if (ndims < 2)
		quit("readBMLandInfo: not enough dimensions\n");
	if (nvars < 2)
		quit("readBMLandInfo: not enough variables\n");

	/* Check dimensions are as expected */
	if ((t_id = ncdimid(fid, "t")) == -1)
		quit("readBMLandInfo: no t dimension\n");
	if (t_id != recdim)
		quit("readBMLandInfo: t dimension not unlimited\n");
	if ((b_id = ncdimid(fid, "b")) == -1)
		quit("readBMLandInfo: no b dimension\n");

	/* Get dimension sizes and check against geometry */
	ncdiminq(fid, b_id, NULL, &n);
	if (n != bm->nbox)
		quit("readBMLandInfo: Number of boxes (b) doesn't match geometry\n");

	/* Find and count land related variables */
	dims[0] = t_id;
	dims[1] = b_id;
	idlist = i_alloc1d(nvars);

	n = ncvarfind(fid, 2, dims, "bmtype", "terrestrial", idlist);
	if (n < 1)
		warn("readBMLandInfo: No terrestrial variables\n");
	if (verbose > 1)
		fprintf(stderr, "readBMLandInfo: %ld terrestrial variables\n", n);


	/* Allocate space for land info */
	if ((bm->landinfo = (EpiInfo *) malloc((size_t)n * sizeof(EpiInfo))) == NULL)
		quit("readBMLandInfo: Can't allocate memory for info\n");

	/* Clear tracer info */
	memset(bm->landinfo, 0, (size_t)n * sizeof(EpiInfo));

	/* Loop over the terrestrial (land) variables to read info */
	for (i = 0; i < n; i++) {
		nc_type dt;
		int nd;
		int d[3];
		int na;
		int len;

		/* Land tracer number (must be equal to i) */
		bm->landinfo[i].n = i;

		/* Variable name */
		ncvarinq(fid, idlist[i], bm->landinfo[i].name, &dt, &nd, d, &na);
		/* Variable attributes */


		/* Units */
		ncattget(fid, idlist[i], "units", bm->landinfo[i].units);

		/* Long name */
		ncattget(fid, idlist[i], "long_name", bm->landinfo[i].long_name);

		/* Data type */
		ncattinq(fid, ncvarid(fid, bm->landinfo[i].name), "dtype", &dt, &len);
		if (nctypelen(dt) != sizeof(bm->landinfo[i].dtype))
			quit("readBMLandInfo: dtype attribute wrong type (%s)\n", bm->landinfo[i].name);
		ncattget(fid, idlist[i], "dtype", &bm->landinfo[i].dtype);

		/* Sum type */
		ncattinq(fid, ncvarid(fid, bm->landinfo[i].name), "sumtype", &dt, &len);
		if (nctypelen(dt) != sizeof(bm->landinfo[i].sumtype))
			quit("readBMLandInfo: sumtype attribute wrong type (%s)\n", bm->landinfo[i].name);
		ncattget(fid, idlist[i], "sumtype", &bm->landinfo[i].sumtype);

		bm->landinfo[i].flagid = 1;

		if (verbose > 1)
			fprintf(stderr, "readBMLandInfo: %s info read\n", bm->landinfo[i].name);
	}

	/* Store number of land tracers */
	bm->nland = (int)n;

	/* Free land variable id array */
	i_free1d(idlist);
}

/**** Free any memory allocated during a readBMLandInfo call.
 ****/
void freeBMLandInfo(MSEBoxModel *bm) {
	free(bm->landinfo);
}

/*********//**
 Routine to write land tracer info to general data file
 ***********/
void writeBMLandInfo(int fid, MSEBoxModel *bm, int dtype) {
	int i = 0;
	nc_type dt;
	int dim[2];
	int vid;

	/* set netCDF data type */
	dt = NC_DOUBLE;

	/* Get dimension ids */
	dim[0] = ncdimid(fid, "t");
	dim[1] = ncdimid(fid, "b");

	/* Loop over the epi variables */
	for (i = 0; i < bm->nland; i++) {

		if (dtype == bm->landinfo[i].dtype && bm->landinfo[i].flagid) {

			/* Define the variable */
			vid = ncvardef(fid, bm->landinfo[i].name, dt, 2, dim);

			/* Set the attributes */
			ncattput(fid, vid, "bmtype", NC_CHAR, (int) strlen("terrestrial") + 1, "terrestrial");
			ncattput(fid, vid, "units", NC_CHAR, (int) strlen(bm->landinfo[i].units) + 1, bm->landinfo[i].units);
			ncattput(fid, vid, "long_name", NC_CHAR, (int) strlen(bm->landinfo[i].long_name) + 1, bm->landinfo[i].long_name);
			ncattput(fid, vid, "dtype", NC_LONG, 1, &bm->landinfo[i].dtype);
			ncattput(fid, vid, "sumtype", NC_LONG, 1, &bm->landinfo[i].sumtype);

			if (verbose > 1)
				fprintf(stderr, "writeBMLandInfo: %s info written\n", bm->landinfo[i].name);
		}
	}
}

/*******************************************************************//**
 Routine to read the land tracer data from a netCDF file. This routine
 assumes that the land tracer information in the MSEBoxModel is valid and
 corresponds with the netCDF file. Minimal checking for consistency
 is done here.
 *********************************************************************/
void readBMLandData(int fid, int dump, MSEBoxModel *bm) {

	int i = 0;
	double *val;
	long start[2];
	long count[2];
	long n = 0;
	int b;

	/* Set netCDF library error handling */
	ncopts = NC_VERBOSE | NC_FATAL;

	/* Check that requested dump exists in file */
	ncdiminq(fid, ncdimid(fid, "t"), NULL, &n);
	if (dump >= n)
		quit("readBMLandData: dump %d not in file (%ld records)\n", dump, n);

	/* Allocate temporary storage for one epi variable */
	val = (double *) alloc1d(bm->nbox);

	/* Set indices for reading variables */
	start[0] = dump;
	start[1] = 0;
	count[0] = 1;
	count[1] = bm->nbox;

	/* Loop over each variable */
	for (i = 0; i < bm->nland; i++) {

		/* Read data */
		ncvarread(fid, bm->landinfo[i].name, sizeof(double), start, count, (void *) val);

		/* Move data to model storage */
		for (b = 0; b < bm->nbox; b++) {
			if(bm->boxes[b].type == LAND){
				bm->landtr[b][i] = (double) val[b];
			}
		}
	}
	/* Free temporary storage */
	free1d((double *) val);

}

/*******************************************************************//**
 Routine to write the land tracer data to a netCDF file. This routine
 assumes that the land tracer information in the MSEBoxModel is valid and
 corresponds with the netCDF file. Minimal checking for consistency
 is done here.
 *********************************************************************/
void writeBMLandData(int fid, int dump, MSEBoxModel *bm, int dtype) {

	int i = 0;
	double *val;
	long start[2];
	long count[2];
	int vid;
	int b;
	int n;

	if (verbose > 0)
		fprintf(stderr, "Entering writeBMEpiData\n");

	/* Set netCDF library error handling */
	ncopts = NC_VERBOSE | NC_FATAL;

	/* Allocate temporary storage for one epi variable */
	val = (doubleINPUT *) alloc1dInput(bm->nbox);

	/* Set indices for writing epis */
	start[0] = dump;
	start[1] = 0;
	count[0] = 1;
	count[1] = bm->nbox;

	/* Loop over each land tracer */
	for (i = 0; i < bm->nland; i++) {

		if (dtype == bm->landinfo[i].dtype && bm->landinfo[i].flagid) {
			/* Get netCDF variable id */
			vid = ncvarid(fid, bm->landinfo[i].name);

			/* Check double compatibility */
			if ((n = ncvarsize(fid, vid)) != sizeof(doubleINPUT))
				quit("writeBMLandData: %s has %d bytes per value in file, code compiled with %d\n", bm->landinfo[i].name, n, sizeof(doubleINPUT));

			/* Get land data from model storage */
			for (b = 0; b < bm->nbox; b++){
				val[b] = (double) bm->landtr[b][i];
			}

			/* Write data */
			ncvarput(fid, vid, start, count, val);
		}

	}

	/* Free temporary storage */
	free1dInput((doubleINPUT *) val);

}
