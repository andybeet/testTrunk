/*********************************************************************

    File:           ptrack.c
    
    Created:        Mon Jan 15 15:47:40 EST 1996
    
    Author:         Stephen Walker
                    CSIRO Division of Oceanography
    
    Purpose:        Routines for particle tracking
    
    Arguments:      See below

    Returns:        void
    
    Revisions:      none

    $Id: ptrack.c 3369 2012-08-29 06:16:46Z gor171 $

*********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <netcdf.h>
#include <sjwlib.h>


/** Read an array of particles from a netCDF particle file
  * at a specified record.
  *
  * @param name netCDF particle filename.
  * @param rec record number.
  * @param np returned number of particles.
  * @param p returned pointer to a particle structure.
  * @param t returned time for specified record (NULL if not required).
  * @param t_units returned time units (NULL if not required).
  * @param ndump returned number of dumps in particle file (NULL if not required).
  */
void
readparticles(char *name, int rec, long int *np, Particle **p,
	      double *t, char *t_units, int *ndump)
{
    int fid;
    int ndims;       /* Number of dimensions */
    int nvars;       /* Number of variables */
    int natts;       /* Number of attributes */
    int recdim;      /* ID of unlimited dimension */
    nc_type dt;      /* type of variable */
    int t_did;       /* Time dimension id */
    int p_did = -1;  /* p dimension id */
    int dims[NC_MAX_DIMS];
    size_t start[2];
    size_t count[2];
    //size_t n;
    int n;
    int t_vid;
    int x_vid;
    int y_vid;
    int z_vid;
    int f_vid;
    doubleINPUT *b;
    short *f;
    size_t nrec;

    /* Open the netCDF file */
    if( (nc_open(name,NC_NOWRITE,&fid)) != NC_NOERR )
    	quit("readparticles: Can't open %s\n",name);
    
    /* Inquire about this file */
    if( nc_inq(fid,&ndims,&nvars,&natts,&recdim) == -1 )
    	quit("readparticles: nc_inq failed\n");
    if( ndims != 2 )
    	quit("readparticles: not enough dimensions\n");
    if( nvars < 5 )
    	quit("readparticles: not enough variables\n");

    /* Check dimensions are as expected */
    if( (nc_inq_dimid(fid,"t", &t_did)) != NC_NOERR )
    	quit("readparticles: no t dimension");
    if( t_did != recdim )
    	quit("readparticles: t dimension not unlimited\n");
    if( (nc_inq_dimid(fid,"n",&p_did)) != NC_NOERR )
    	quit("readparticles: no n dimension");
    /* Get number of particles */
    nc_inq_dimlen(fid,p_did,(size_t *)&n);

    /* Check that there is a variable called "t" */
    if( (t_vid=nc_varid(fid,"t")) < 0 )
	quit("readparticles: No t variable\n");
    nc_inq_var(fid,t_vid,NULL,&dt,&ndims,dims,&natts);
    if( ndims != 1 || dims[0] != t_did )
	quit("readparticles: t variable has wrong dimensions\n");
    if( dt != NC_DOUBLE )
	quit("readparticles: t variable must have type NC_DOUBLE\n");
    /* Get number of records and check that requested record exists */
    nc_inq_dimlen(fid,t_did,&nrec);
    if( rec >= (int)nrec )
    	quit("readparticles: Record %d not in %s (only %d records)\n",rec,name,nrec);
    if( ndump )
	*ndump = (int)nrec;

    /* Check that there is a variable called "x" */
    if( (x_vid=nc_varid(fid,"x")) < 0 )
	quit("readparticles: No x variable\n");
    nc_inq_var(fid,x_vid,NULL,&dt,&ndims,dims,&natts);
    if( ndims != 2 || dims[0] != t_did || dims[1] != p_did )
	quit("readparticles: x variable has wrong dimensions\n");
    if( dt != NC_DOUBLE )
	quit("readparticles: x variable must have type NC_DOUBLE\n");

    /* Check that there is a variable called "y" */
    if( (y_vid=nc_varid(fid,"y")) < 0 )
	quit("readparticles: No y variable\n");
    nc_inq_var(fid,y_vid,NULL,&dt,&ndims,dims,&natts);
    if( ndims != 2 || dims[0] != t_did || dims[1] != p_did )
	quit("readparticles: y variable has wrong dimensions\n");
    if( dt != NC_DOUBLE )
	quit("readparticles: y variable must have type NC_DOUBLE\n");

    /* Check that there is a variable called "z" */
    if( (z_vid=nc_varid(fid,"z")) < 0 )
	quit("readparticles: No z variable\n");
    nc_inq_var(fid,z_vid,NULL,&dt,&ndims,dims,&natts);
    if( ndims != 2 || dims[0] != t_did || dims[1] != p_did )
	quit("readparticles: z variable has wrong dimensions\n");
    if( dt != NC_DOUBLE )
	quit("readparticles: z variable must have type NC_DOUBLE\n");

    /* Check that there is a variable called "flag" */
    if( (f_vid=nc_varid(fid,"flag")) < 0 )
	quit("readparticles: No flag variable\n");
    nc_inq_var(fid,f_vid,NULL,&dt,&ndims,dims,&natts);
    if( ndims != 2 || dims[0] != t_did || dims[1] != p_did )
	quit("readparticles: flag variable has wrong dimensions\n");
    if( dt != NC_SHORT )
	quit("readparticles: flag variable must have type NC_SHORT\n");

    /* Allocate space, if not already done */
    if( *p == NULL ) {
	*np = (long)n;
	if( (*p=(Particle *)malloc((size_t)(*np)*sizeof(Particle))) == NULL )
	    quit("readparticles: not enough memory for particles\n");
    }
    else if( *np != (int)n )
    	quit("readparticles: Number of particles doesn't match space already allocated\n");
    
    b = alloc1dInput(*np);

    /* Read time information if required */
    start[0] = (size_t)rec;
    count[0] = 1;
    if( t )
	nc_get_vara_double(fid,t_vid,start,count, (doubleINPUT *)b);
    for(n = 0; n < rec; n++){
    	t[n] = (double)b[n];
    }
    if( t_units )
	nc_get_att_text(fid,t_vid,"units",t_units);

    /* Read the particle data */
   // b = alloc1dInput(*np);
    f = s_alloc1d(*np);
    start[0] = (size_t)rec; start[1] = 0;
    count[0] = 1; count[1] = (size_t)*np;
    nc_get_vara_double(fid,x_vid,start,count, (doubleINPUT *)b);
    for(n=0; n<*np; n++)
	(*p)[n].e1 = (double)b[n];
    nc_get_vara_double(fid,y_vid,start,count, (doubleINPUT *)b);
    for(n=0; n<*np; n++)
	(*p)[n].e2 = (double)b[n];
    nc_get_vara_double(fid,z_vid,start,count,(doubleINPUT *)b);
    for(n=0; n<*np; n++)
	(*p)[n].e3 =(double) b[n];
    nc_get_vara_short(fid,f_vid,start,count,f);
    for(n=0; n<*np; n++)
	(*p)[n].flag = f[n];
    free1dInput(b);
    s_free1d(f);
    
    /* Close the file */
    nc_close(fid);
}

/** Create particle tracking output file.
  *
  * @param name particle filename.
  * @param np number of particles per record dump.
  * @param t_units time units of each dump.
  * @return netCDF file descriptor.
  */
int createptfile(char *name, long int np, char *t_units)
{
   int  ncid;			/* netCDF id */
   /* dimension ids */
   int  t_dim, n_dim;
   /* variable ids */
   int  t_id, x_id, y_id, z_id, flag_id;
   /* variable shapes */
   int dims[2];

   /* enter define mode */
   if (nc_create(name, NC_NOCLOBBER, &ncid) != NC_NOERR)
      quit("createptfile: Unable to create the particle file '%s'\n", name);

   /* define dimensions */
   nc_def_dim(ncid, "t", NC_UNLIMITED, &t_dim);
   nc_def_dim(ncid, "n", (size_t)np, &n_dim);

   /* define variables */
   dims[0] = t_dim;
   nc_def_var (ncid, "t", NC_DOUBLE, 1, dims, &t_id);
   dims[1] = n_dim;
   nc_def_var (ncid, "x", NC_DOUBLE, 2, dims, &x_id);
   nc_def_var (ncid, "y", NC_DOUBLE, 2, dims, &y_id);
   nc_def_var (ncid, "z", NC_DOUBLE, 2, dims, &z_id);
   nc_def_var (ncid, "flag", NC_SHORT, 2, dims, &flag_id);

   /* Time variable units attribute */
   nc_put_att_text(ncid,t_id,"units",strlen(t_units),t_units);

   /* leave define mode */
   nc_enddef(ncid);

   return(ncid);
}


/** Write particles at a known record into the particle file.
  *
  * @param fid file descriptor to open and writable netCDF
  *            particle file.
  * @param rec record number to write.
  * @param t time at specified record.
  * @param t np number of particles to write.
  * @param p array of particles to write.
  */
void
writeparticles(int fid, int rec, double t, long int np, Particle *p)
{
	doubleINPUT *b;
    short *f;
    size_t start[2];
    size_t count[2];
    long n;

    /* Write time value */
    start[0] = (size_t)rec;
    count[0] = 1;
    nc_put_vara_double(fid,nc_varid(fid,"t"),start,count,(void *)&t);

    start[1] = 0;
    count[1] = (size_t)np;

    /* Allocate buffers */
    b = alloc1dInput(np);
    f = s_alloc1d(np);

    /* Transfer and write data */
    for(n=0; n<np; n++)
	b[n] = p[n].e1;
    nc_put_vara_double(fid,nc_varid(fid,"x"),start,count, (doubleINPUT *)b);

    for(n=0; n<np; n++)
	b[n] = p[n].e2;
    nc_put_vara_double(fid,nc_varid(fid,"y"),start,count,  (doubleINPUT *)b);

    for(n=0; n<np; n++)
	b[n] = p[n].e3;
    nc_put_vara_double(fid,nc_varid(fid,"z"),start,count,  (doubleINPUT *)b);

    for(n=0; n<np; n++)
	f[n] = p[n].flag;
    nc_put_vara_short(fid,nc_varid(fid,"flag"),start,count,f);

    d_free1d(b);
    s_free1d(f);

    nc_sync(fid);
}
