/*

    File:           netcdf.c

    Created:        Fri Nov 25 10:16:56 EST 1994

    Author:         Stephen Walker
                    CSIRO Division of Oceanography

    Purpose:        Useful netCDF related routines

    Arguments:      See each routine below

    Returns:        void

	 Revisions:      22-May-2009 Bec Gorton
					Added the sjw_ncopen function to warp the ncopen function. This calls fopen
					to check that the file exists before ncopen is called. The ncopen function
					now crashes if we try to open a file that does not exist. It should just return NULL.

					26-May-2011 Bec Gorton
					Removed the sjw_ncopen function - its been replaced by the Util_ncopen function in the util lib.

    $Id: netcdf.c 3226 2012-06-21 02:24:47Z gor171 $


    Changes:


    21-06-2012 Bec Gorton
    Added functions to check the values read in from the netcdfs are not undefined.

*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <netcdf.h>
#include "sjwlib.h"

/** Locate a netCDF dimension name and return the dimension id.
  * A replacement function to mimic the NetCDF 2.X function ncdimid.
  * @param fid file descriptor of an open netcdf file.
  * @param name netcdf dimension name.
  * @return dimension identifier.
  */
int nc_dimid(int fid, const char* name)
{
int id = -1;

	if (nc_inq_dimid(fid, name, &id) != NC_NOERR)
	   id = -1;

	return id;
}


/** Locate a netCDF variable name and return the variable id.
  * A replacement function to mimic the NetCDF 2.X function ncvarid.
  * @param fid file descriptor of an open netcdf file.
  * @param name netcdf variable name.
  * @return variable identifier.
  */
int nc_varid(int fid, const char* name)
{
int id = -1;

	if (nc_inq_varid(fid, name, &id) != NC_NOERR)
	   id = -1;

	return id;
}


/** Finds which variables in a netCDF file have a certain
  * set of dimensions and a particular attribute with a particular
  * value. Arguments are as follows:
  *
  * @param fid file descriptor of an open netcdf file.
  * @param nvdims number of dimensions for matching variables. Only
  *               checked if > 0.
  * @param vdims list of dimension ids for matching variables. Only
  *              checked if not NULL.
  * @param attr name of attribute. Only checked if not NULL.
  * @param attval value of attribute (must be a string).
  *               Only checked if not NULL.
  * @param list Output list of variable ids which have the above
  * 	        properties
  * @return number of matching variables found.
  */
int
ncvarfind(int fid, int nvdims, int *vdims, char *attr, char *attval, int *list)
{
    int ndims;
    int nvars;
    int ngatts;
    int recdim;
    int id;
    int i;
    int n;

    /* Find out file properties */
    nc_inq(fid,&ndims,&nvars,&ngatts,&recdim);

    if( ndims < nvdims )
	return(0);

    /* Loop over the variables */
    n = 0;
    for(id=0; id<nvars; id++) {
		char name[MAXLINELEN];
		nc_type dt;
		int nd;
		int d[NC_MAX_VAR_DIMS];
		int na;
		char aval[MAXLINELEN];
		int attnum;
		int match;

		match = 1;
		memset(aval,0,MAXLINELEN);
		/* Check number of dimensions */
	        nc_inq_var(fid,id,name,&dt,&nd,d,&na);
		if( (nvdims>0) && (nd!=nvdims) ) match = 0;

		/* Check each dimension */
		for(i=0; vdims!=NULL && match && i<nd; i++)
		    if( (vdims[i]>=0) && (d[i]!=vdims[i]) ) match = 0;

		/* Check attribute existence */
		if( match && (attr != NULL) &&
		    (nc_inq_attid(fid,id,attr, &attnum) != NC_NOERR) )
		    match = 0;

		/* Check attribute value */
		if( match && (attr != NULL) && (attval != NULL)) {
		    if (nc_get_att_text(fid,id,attr,aval) == NC_NOERR) {
		       if (strcmp(aval,attval) != 0 )
		          match = 0;
		    }
		    else
		       match = 0;
		}

	        /* Store id if we have a match */
		if( match ) {
		    list[n] = id;
		    n++;
		}
    }

    return(n);
}


/** Determine the number of bytes needed to store a
  * single value of a specified netCDF variable.
  * @param fid file descriptor of an open netcdf file.
  * @param vid id of the variable.
  * @return number of bytes.
  */
int ncvarsize(int fid, int vid)
{
    nc_type datatype;
    size_t size = 0;

    nc_inq_vartype(fid,vid,&datatype);
    switch (datatype) {
    case NC_BYTE:
    case NC_CHAR:
	size=sizeof(char);
	break;

    case NC_SHORT:
	size=sizeof(short);
	break;

    case NC_INT:
	size=sizeof(int);
	break;

    case NC_FLOAT:
	size=sizeof(float);
	break;

    case NC_DOUBLE:
	//size=sizeof(double);
    	size = sizeof(double);
	break;

    case NC_NAT:
	abort();
    }

    return((int)size);
}

/**
 * Check that the data read in from the netcdf file is not the default fill value. If it is general some
 * nice error messages and quit.
 *
 */
void checkNetCDFData2D( char *functionName, char *varName, double **data, int size1, int size2){
	int b, i;
	for (b = 0; b < size1; b++) {
		for (i = 0; i < size2; i++) {
			if(data[b][i] == NC_FILL_DOUBLE){

				fprintf(stderr, "\n\nERROR: %s - Undefined data for tracer %s in box %d, layer %d in your initial conditions file.\n\n",
						functionName, varName, b, i);
				fprintf(stderr, "Please make sure all variables in your initial conditions file either have a fill value or have defined data values\n");
				quit("");
			}
		}
	}
}

/**
 * Check that the data read in from the netcdf file is not the default fill value. If it is general some
 * nice error messages and quit.
 *
 */
void checkNetCDFData1D(char *functionName, char *varName, double *data, int size1){
	int b;
	for (b = 0; b < size1; b++) {
		if(data[b] == NC_FILL_DOUBLE){

			fprintf(stderr, "\n\nERROR: %s - Undefined data for tracer %s in box %d in your initial conditions file.\n\n",
					functionName, varName, b);
			fprintf(stderr, "Please make sure all variables in your initial conditions file either have a fill value or have defined data values\n");
			quit("");
		}
	}
}

/**
 * Check that the data read in from the netcdf file is not the default fill value. If it is general some
 * nice error messages and quit.
 *
 */
void checkNetCDFData1DShort(char *functionName, char *varName, short *data, int size1){
	int b;
	for (b = 0; b < size1; b++) {
		if(data[b] == NC_FILL_SHORT){

			fprintf(stderr, "\n\nERROR: %s - Undefined data for tracer %s in box %d in your initial conditions file.\n\n",
					functionName, varName, b);
			fprintf(stderr, "Please make sure all variables in your initial conditions file either have a fill value or have defined data values\n");
			quit("");
		}
	}
}


/** Reads a hyperslab of a single netcdf variable, checking that
  * the variable type matches a specified size.
  *
  * This routine mimics the functionality available below netCDF 3.
  * It is best to avoid this function in any new code.
  *
  * @param fid file descriptor of an open netcdf file.
  * @param name netcdf variable name.
  * @param size variable size.
  * @param start pointer to array of hyperslab start positions.
  * @param count pointer to array of hyperslab sizes.
  * @param buf pointer to the array in which the values will be read.
  */
void ncvarread(int fid, char *name, int size, long int *start, long int *count, void *buf)
{
    int n;
    nc_type datatype;
    int vid = nc_varid(fid, name);
    int returnValue = 0;

    if( (n=ncvarsize(fid,vid)) != size )
        quit("ncvarread: %s has %d bytes per value in file, code expects %d\n",n,size);
    nc_inq_vartype(fid,vid,&datatype);

    switch (datatype) {
    case NC_BYTE:
    case NC_CHAR:
    	returnValue = nc_get_vara_schar(fid,vid,(size_t*)start, (size_t*)count,(signed char*)buf);
    	break;

    case NC_SHORT:
    	returnValue = nc_get_vara_short(fid,vid,(size_t*)start, (size_t*)count,(short*)buf);
    	break;

    case NC_INT:
    	returnValue = nc_get_vara_int(fid,vid,(size_t*)start, (size_t*)count,(int*)buf);
    	break;

    case NC_FLOAT:
    	returnValue = nc_get_vara_float(fid,vid,(size_t*)start, (size_t*)count,(float*)buf);
    	break;

    case NC_DOUBLE:
        //nc_get_vara_double(fid,vid,(size_t*)start, (size_t*)count,(double *)buf);
    	returnValue = nc_get_vara_double(fid,vid,(size_t*)start, (size_t*)count,(double *)buf);
    	break;

    case NC_NAT:
    	abort();
    }
    if(returnValue != NC_NOERR){
    	quit("error reading in netcdf data\n");
    }
}


