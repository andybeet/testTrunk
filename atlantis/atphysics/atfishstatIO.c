/*********************************************************************

 File:           atfishstatIO.c

 Created:        Mon June 13 13:54:39 EST 2005

 Author:         Beth Fulton
 CSIRO Marine Research

 Purpose:        Routines which deal with fisheries statistics files. These
 files are in netCDF format and should be read, written
 and created using the routines below.

 Arguments:      See each routine below

 Returns:        See each routine below

 Revisions:      13/6/2005 EA Fulton
 Created from atepiIO.c

 1/9/2005 EA Fulton
 Added initialisation

 *********************************************************************/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <sjwlib.h>
#include <netcdf.h>
#include <atlantisboxmodel.h>
#include <atUtilLib.h>

/*********************************************************************
 Routine to read the model fstats from a netCDF file.
 This routine assumes that the file has the following
 format (roughly written in CDL).

 dimensions:

 t = unlimited;    // Time dimension
 b = N;            // Number of boxes
 z = N;	    // Number of layers (water column + sediment)

 variables:

 // A typical fisheries statistic entry
 double Effort(t,b);                      // Effort values
 Effort:bmtype = "statistics";          // This is a fisheries statistics variable
 Effort:units = "Hours";                // Units for Effort
 Effort:long_name = "Fisheries effort"; // Long name of Effort
 Effort:dtype = 0;                      // Flag => general data type (1 if fisheries data)

 // The netCDF file has global attributes as follows
 :title = "Description of this data set";
 :geometry = "Name of geometry file";
 :parameters = "Name of parameter file";
 :wcnz = N;     // Number of water column layers
 :sednz = N;    // Number of sediment layers

 The routine finds out how many fisheries statistics there
 are, allocates memory and reads the information and
 values for each statistic.

 *********************************************************************/
void readBMFisheriesInfo(int fid, MSEBoxModel *bm) {
	int ndims = 0; /* Number of dimensions */
	int nvars = 0; /* Number of variables */
	int ngatts = 0; /* Number of global attributes */
	int recdim = 0; /* ID of unlimited dimension */
	int t_id = 0; /* Dimension ids */
	int b_id = 0;
	int dims[MAX_NC_DIMS];
	int *idlist = 0;
	long n;
	int i, b;

	/* Set netCDF library error handling */
	ncopts = NC_VERBOSE | NC_FATAL;

	if (verbose > 1)
        printf( "Entering readBMFisheriesInfo\n");

	/* Inquire about this file */
	ncinquire(fid, &ndims, &nvars, &ngatts, &recdim);
	if (ndims < 2)
		quit("readBMFisheriesInfo: not enough dimensions\n");
	if (nvars < 2)
		quit("readBMFisheriesInfo: not enough variables\n");

	/* Check dimensions are as expected */
	if ((t_id = ncdimid(fid, "t")) == -1)
		quit("readBMFisheriesInfo: no t dimension\n");
	if (t_id != recdim)
		quit("readBMFisheriesInfo: t dimension not unlimited\n");
	if ((b_id = ncdimid(fid, "b")) == -1)
		quit("readBMFisheriesInfo: no b dimension\n");

	/* Get dimension sizes and check against geometry */
	ncdiminq(fid, b_id, NULL, &n);
	if (n != bm->nbox)
		quit("readBMFisheriesInfo: Number of boxes (b) doesn't match geometry\n");

	/* Find and count fisheries statistics */
	dims[0] = t_id;
	dims[1] = b_id;
	idlist = i_alloc1d(nvars);
	n = ncvarfind(fid, 2, dims, "bmtype", "statistics", idlist);
	if (n < 1)
		warn("readBMFisheriesInfo: No fisheries statistics\n");
	if (verbose > 1)
        printf( "readBMFisheriesInfo: %ld fisheries statistics\n", n);

	/* Allocate space for statistics info */
	if ((bm->finfo = (FstatInfo *) malloc((size_t)n * sizeof(FstatInfo))) == NULL)
		quit("readBMFisheriesInfo: Can't allocate memory for info\n");

	/* Clear tracer info */
	memset(bm->finfo, 0, (size_t)n * sizeof(FstatInfo));

	/* Loop over the fisheries statistics to read info */
	for (i = 0; i < n; i++) {
		nc_type dt;
		int nd;
		int d[3];
		int na;
		int len;

		/* Fstat number (must be equal to i) */
		bm->finfo[i].n = i;
		/* Fstat name */
		ncvarinq(fid, idlist[i], bm->finfo[i].name, &dt, &nd, d, &na);
		/* Fstat attributes */
		ncattget(fid, idlist[i], "units", bm->finfo[i].units);
		ncattget(fid, idlist[i], "long_name", bm->finfo[i].long_name);

		ncattinq(fid, ncvarid(fid, bm->finfo[i].name), "dtype", &dt, &len);
		if (nctypelen(dt) != sizeof(bm->finfo[i].dtype))
			quit("readBMFisheriesInfo: dtype attribute wrong type (%s)\n", bm->finfo[i].name);
		ncattget(fid, idlist[i], "dtype", &bm->finfo[i].dtype);

		ncattinq(fid, ncvarid(fid, bm->finfo[i].name), "sumtype", &dt, &len);
		if (nctypelen(dt) != sizeof(bm->finfo[i].sumtype))
			quit("readBMFisheriesInfo: sumtype attribute wrong type (%s)\n", bm->finfo[i].name);
		ncattget(fid, idlist[i], "sumtype", &bm->finfo[i].sumtype);

		bm->finfo[i].flagid = 1;

		if (verbose > 1)
			fprintf(stderr, "readBMFisheriesInfo: %d entry %s info read\n", i, bm->finfo[i].name);
		/*		fprintf(stderr,"with units: %s\n",bm->finfo[i].units);
		 fprintf(stderr,"and long name: %s\n",bm->finfo[i].long_name);
		 */
	}

	/* Store number of statistics */
	bm->nfstat = (int)n;

	fprintf(stderr, "Read in %d fisheries statistics\n", bm->nfstat);

	/* Allocate memory for fishstatics */
	bm->fishstat = (double **) alloc2d(bm->nfstat, bm->nbox);

	/* Set up pointers to data from each box */
	for (b = 0; b < bm->nbox; b++) {
		bm->boxes[b].fishstat = bm->fishstat[b];
	}

	/* Initialise values */
	for (i = 0; i < bm->nfstat; i++) {
		for (b = 0; b < bm->nbox; b++) {
			bm->fishstat[b][i] = 0.0;
		}
	}

	/* Free fisheries statistics id array */
	i_free1d(idlist);
}

/**** Free any memory allocated during a readBMFisheriesInfo call.
 ****/
void freeBMFisheriesInfo(MSEBoxModel *bm) {
	free(bm->finfo);
}

/***************************************************************************************
 Routine to write fisheries statistics info (finfo) to general data file
 ***********/
void writeBMFisheriesInfo(int fid, MSEBoxModel *bm, int dtype) {
	int i = 0;
	nc_type dt;
	int dim[2];

	/* set netCDF data type */
	dt = NC_DOUBLE;

	/* Get dimension ids */
	dim[0] = ncdimid(fid, "t");
	dim[1] = ncdimid(fid, "b");

	/* Loop over the fisheries statistics */
	for (i = 0; i < bm->nfstat; i++) {
		int vid;

		if (dtype == bm->finfo[i].dtype && bm->finfo[i].flagid) {

			/* Define the variable */
			vid = ncvardef(fid, bm->finfo[i].name, dt, 2, dim);

			/* Set the attributes */
			ncattput(fid, vid, "bmtype", NC_CHAR, (int) strlen("statistics") + 1, "statistics");
			ncattput(fid, vid, "units", NC_CHAR, (int) strlen(bm->finfo[i].units) + 1, bm->finfo[i].units);
			ncattput(fid, vid, "long_name", NC_CHAR, (int) strlen(bm->finfo[i].long_name) + 1, bm->finfo[i].long_name);
			ncattput(fid, vid, "dtype", NC_LONG, 1, &bm->finfo[i].dtype);
			ncattput(fid, vid, "sumtype", NC_LONG, 1, &bm->finfo[i].sumtype);

			if (verbose > 1)
				fprintf(stderr, "writeBMFisheriesInfo: %d entry %s info written, dtype = %d\n", i, bm->finfo[i].name, bm->finfo[i].dtype);
			/*		fprintf(stderr,"with units: %s\n",bm->finfo[i].units);
			 fprintf(stderr,"and long name: %s\n",bm->finfo[i].long_name);
			 */

		}
	}
}

/*********************************************************************
 Routine to write the fisheries statistics to a netCDF file. This routine
 assumes that the finfo information in the MSEBoxModel is valid and
 corresponds with the netCDF file. Minimal checking for consistency
 is done here.
 *********************************************************************/
void writeBMFisheriesData(int fid, int dump, MSEBoxModel *bm, int dtype) {
	int i = 0;
	doubleINPUT *val;
	long start[2];
	long count[2];

	if (verbose > 0)
		fprintf(stderr, "Entering writeBMFisheriesData\n");

	/* Set netCDF library error handling */
	ncopts = NC_VERBOSE | NC_FATAL;

	/* Allocate temporary storage for one statistics */
	val = (doubleINPUT *) alloc1dInput(bm->nbox);

	/* Set indices for writing statistics */
	start[0] = dump;
	start[1] = 0;
	count[0] = 1;
	count[1] = bm->nbox;

	/* Loop over each statistic */
	for (i = 0; i < bm->nfstat; i++) {
		int vid;
		int b;
		int n;

		if (dtype == bm->finfo[i].dtype && bm->finfo[i].flagid) {

			/* Get netCDF variable id */
			vid = ncvarid(fid, bm->finfo[i].name);

			/* Check double compatibility */
			if ((n = ncvarsize(fid, vid)) != sizeof(doubleINPUT))
				quit("writeBMFisheriesData: %s has %d bytes per value in file, code compiled with %d\n", bm->finfo[i].name, n, sizeof(doubleINPUT));

			/* Get fisheries statistics from model storage */
			for (b = 0; b < bm->nbox; b++) {
				val[b] = (doubleINPUT) bm->fishstat[b][i];
			}

			/* Write data */
			ncvarput(fid, vid, start, count, val);
		}
	}

	/* Free temporary storage */
	free1dInput((doubleINPUT *) val);

}
