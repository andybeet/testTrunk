/*********************************************************************

 File:           atdiagIO.c

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

/*******************************************************************//**
 Routine to read the model diagnostics from a netCDF file.
 This routine assumes that the file has the following
 format (roughly written in CDL).

 dimensions:

 t = unlimited;    // Time dimension
 b = N;            // Number of boxes
 z = N;	    // Number of layers (water column + sediment)

 variables:

 // A typical diagnostic entry
 double DinoProd(t,b);                           // DinoProd values
 weed:bmtype = "diagnostic";                   // This is a diagnostic output
 weed:units = "mg N s-1";                      // Units for DinoProd
 weed:long_name = "Dinoflagellate production"; // Long name of DinoProd
 weed:dtype = 0;                               // Flag => general data type (1 if fisheries data)

 // The netCDF file has global attributes as follows
 :title = "Description of this data set";
 :geometry = "Name of geometry file";
 :parameters = "Name of parameter file";
 :wcnz = N;     // Number of water column layers
 :sednz = N;    // Number of sediment layers

 The routine finds out how many diagnostic variables there
 are, allocates memory and reads the information and
 values for each diagnostic.

 *********************************************************************/
//void
//readBMDiagInfo(int fid, MSEBoxModel *bm)
//{
//    int ndims = 0;      /* Number of dimensions */
//    int nvars = 0;      /* Number of variables */
//    int ngatts = 0;     /* Number of global attributes */
//    int recdim = 0;     /* ID of unlimited dimension */
//    int t_id = 0;       /* Dimension ids */
//    int b_id = 0;
//    int dims[MAX_NC_DIMS];
//    int *idlist = 0;
//    long n = 0;
//    int i = 0;//, b = 0;
//
//    /* Set netCDF library error handling */
//    ncopts = NC_VERBOSE | NC_FATAL;
//
//	if( verbose > 1 )
//	fprintf(stderr,"Entering readBMDiagInfo\n");
//
//    /* Inquire about this file */
//    ncinquire(fid,&ndims,&nvars,&ngatts,&recdim);
//    if( ndims < 2 )
//        quit("readBMDiagInfo: not enough dimensions\n");
//    if( nvars < 2 )
//        quit("readBMDiagInfo: not enough variables\n");
//
//    /* Check dimensions are as expected */
//    if( (t_id = ncdimid(fid,"t"))  == -1 )
//        quit("readBMDiagInfo: no t dimension\n");
//    if( t_id != recdim )
//        quit("readBMDiagInfo: t dimension not unlimited\n");
//    if( (b_id = ncdimid(fid,"b"))  == -1 )
//        quit("readBMDiagInfo: no b dimension\n");
//
//    /* Get dimension sizes and check against geometry */
//    ncdiminq(fid,b_id,NULL,&n);
//    if( n != bm->nbox )
//        quit("readBMDiagInfo: Number of boxes (b) doesn't match geometry\n");
//
//    /* Find and count diagnostics */
//    dims[0] = t_id;
//    dims[1] = b_id;
//    idlist = i_alloc1d(nvars);
//    n = ncvarfind(fid,2,dims,"bmtype","diagnostic",idlist);
//    if( n < 1 )
//		warn("readBMDiagInfo: No diagnostics\n");
//    if( verbose > 1 )
//		fprintf(stderr,"readBMDiagInfo: %ld diagnostics\n",n);
//
//    /* Allocate space for diagnostic info */
//    if( (bm->dinfo=(DiagInfo *)malloc(n*sizeof(DiagInfo))) == NULL )
//	    quit("readBMDiagInfo: Can't allocate memory for info\n");
//
//	/* Clear tracer info */
//    memset(bm->dinfo, 0, n*sizeof(DiagInfo));
//
//    /* Loop over the diagnostics to read info */
//    for(i=0; i<n; i++) {
//		nc_type dt;
//		int nd;
//		int d[3];
//		int na;
//		int len;
//
//		/* Diagnostic number (must be equal to i) */
//		bm->dinfo[i].n = i;
//		/* Diagnostic name */
//		ncvarinq(fid,idlist[i],bm->dinfo[i].name,&dt,&nd,d,&na);
//
//		/* Diagnostic attributes */
//		ncattget(fid,idlist[i],"units",bm->dinfo[i].units);
//		ncattget(fid,idlist[i],"long_name",bm->dinfo[i].long_name);
//
//		ncattinq(fid,ncvarid(fid,bm->dinfo[i].name),"dtype",&dt,&len);
//
//
//		if( nctypelen(dt) != sizeof(double))
//			quit("readBMDiagInfo: dtype attribute wrong type (%s)\n",bm->dinfo[i].name);
//		ncattget(fid,idlist[i],"dtype",&bm->dinfo[i].dtype);
//
//		ncattinq(fid,ncvarid(fid,bm->dinfo[i].name),"sumtype",&dt,&len);
//		if( nctypelen(dt) != sizeof(double))
//			quit("readBMDiagInfo: sumtype attribute wrong type (%s)\n",bm->dinfo[i].name);
//		ncattget(fid,idlist[i],"sumtype",&bm->dinfo[i].sumtype);
//
//		bm->dinfo[i].flagid = 1;
//
//		if( verbose > 1 )
//			fprintf(stderr,"readBMDiagInfo: %d entry %s info read\n",i, bm->dinfo[i].name);
//	/*		fprintf(stderr,"with units: %s\n",bm->dinfo[i].units);
//			fprintf(stderr,"and long name: %s\n",bm->dinfo[i].long_name);
//	*/
//    }
//
//    /* Store number of diagnostics */
////    bm->ndiag = n;
////
////	fprintf(stderr, "Read in %d diagnostics\n", bm->ndiag);
////
////	/* Allocate memory for diagnostics */
////    switch( sizeof(double) ) {
////	case sizeof(float):
////	    bm->diagnost = (double **)f_alloc2d(bm->ndiag,bm->nbox);
////	    break;
////	case sizeof(double):
////	    bm->diagnost = (double **)alloc2d(bm->ndiag,bm->nbox);
////	    break;
////	default:
////	    quit("setupMSEBoxModel: Unknown size for double\n");
////    }
////
////    /* Set up pointers to data from each box */
////    for(b=0; b<bm->nbox; b++) {
////		bm->boxes[b].diagnost = bm->diagnost[b];
////	}
////
////	/* Initialise values */
////    for(i=0; i<bm->ndiag; i++) {
////		for(b=0; b<bm->nbox; b++) {
////			bm->diagnost[b][i] = 0.0;
////		}
////	}
//
//
//    /* Free diagnostic id array */
//    i_free1d(idlist);
//}

/**//** Free any memory allocated during a readBMDiagInfo call.
 ****/
void freeBMDiagInfo(MSEBoxModel *bm) {
	free(bm->dinfo);
}

/*********//**
 Routine to write diagnostic info to general data file
 ***********/
void writeBMDiagInfo(int fid, MSEBoxModel *bm, int dtype) {
	int i = 0;
	nc_type dt;
	int dim[2];

	/* set netCDF data type */
	dt = NC_DOUBLE;

	/* Get dimension ids */
	dim[0] = ncdimid(fid, "t");
	dim[1] = ncdimid(fid, "b");

	/* Loop over the diagnostics */
	for (i = 0; i < bm->ndiag; i++) {
		int vid;

		if (dtype == bm->dinfo[i].dtype && bm->dinfo[i].flagid) {

			/* Define the variable */
			vid = ncvardef(fid, bm->dinfo[i].name, dt, 2, dim);

			/* Set the attributes */
			ncattput(fid, vid, "bmtype", NC_CHAR, (int) strlen("diagnostic") + 1, "diagnostic");
			ncattput(fid, vid, "units", NC_CHAR, (int) strlen(bm->dinfo[i].units) + 1, bm->dinfo[i].units);
			ncattput(fid, vid, "long_name", NC_CHAR, (int) strlen(bm->dinfo[i].long_name) + 1, bm->dinfo[i].long_name);
			ncattput(fid, vid, "dtype", NC_LONG, 1, &bm->dinfo[i].dtype);
			ncattput(fid, vid, "sumtype", NC_LONG, 1, &bm->dinfo[i].sumtype);

			if (verbose > 1)
				fprintf(stderr, "writeBMDiagInfo: %d entry %s info written\n", i, bm->dinfo[i].name);
			/*		fprintf(stderr,"with units: %s\n",bm->dinfo[i].units);
			 fprintf(stderr,"and long name: %s\n",bm->dinfo[i].long_name);
			 */

		}
	}

}

/*******************************************************************//**
 Routine to write the diagnostics to a netCDF file. This routine
 assumes that the dinfo in the MSEBoxModel is valid and
 corresponds with the netCDF file. Minimal checking for consistency
 is done here.
 *********************************************************************/
void writeBMDiagData(int fid, int dump, MSEBoxModel *bm, int dtype) {
	int i = 0;
	doubleINPUT *val;
	long start[2];
	long count[2];

	if (verbose > 0)
		printf("Entering writeBMDiagData\n");

   /* Set netCDF library error handling */
	ncopts = NC_VERBOSE | NC_FATAL;

	/* Allocate temporary storage for one diagnostic */
	val = (doubleINPUT *) alloc1dInput(bm->nbox);

	/* Set indices for writing diagnostics */
	start[0] = dump;
	start[1] = 0;
	count[0] = 1;
	count[1] = bm->nbox;
    
    /* Loop over each diagnostic */
	for (i = 0; i < bm->ndiag; i++) {
		int vid;
		int b;
		int n;

		if (dtype == bm->dinfo[i].dtype && bm->dinfo[i].flagid) {
			/* Get netCDF variable id */
			vid = ncvarid(fid, bm->dinfo[i].name);

			/* Check double compatibility */
			if ((n = ncvarsize(fid, vid)) != sizeof(doubleINPUT))
				quit("writeBMDiagData: %s has %d bytes per value in file, code compiled with %d\n", bm->dinfo[i].name, n, sizeof(doubleINPUT));

			/* Get diagnostic from model storage */
			for (b = 0; b < bm->nbox; b++) {
				val[b] = (doubleINPUT) bm->diagnost[b][i];
			}

			/* Write data */
			ncvarput(fid, vid, start, count, val);
		}
	}

	/* Free temporary storage */
	free1dInput((doubleINPUT *) val);

}
