/*********************************************************************

    File:           timeseries.h
    
    Created:        Wed Dec 14 12:14:42 EST 1994
    
    Author:         Stephen Walker
                    CSIRO Division of Oceanography
    
    Purpose:        Include file for routines which deal with
		    time series data.
    
    Revisions:      11/06/97 JRW
		    Support for polar grids added.

		    06/08/97 JRW
		    Timeseries is now a convenience facility for
		    Datafile.

    $Id: timeseries.h 2761 2011-07-01 04:35:25Z gor171 $

*********************************************************************/

#ifndef _TIMESERIES_H
#define _TIMESERIES_H




/* Time series types */
#define        TS_ASCII     0  /* Ascii file */
#define        TS_NC0       1  /* netCDF with 0 spatial dimensions */
#define        TS_NC1       2  /* netCDF with 1 spatial dimension */
#define        TS_NC2       3  /* netCDF with 2 spatial dimensions */
#define        TS_NC3       4  /* netCDF with 3 spatial dimensions */


/* In converting timeseries to datafile, only a subset of the more
 * commonly accessed structure variables have been preserved.
 * this facilitates some backwards compatibility */
typedef struct {
    int	 type;		/* Timeseries type */
    char *name;         /* Time series name */
    long nt;            /* Number of time values */
    int ti;             /* Index of time variable */
    char *t_units;      /* Time units string */
    double *t;          /* time values */
    int nv;             /* Number of variables */
    char **varname;     /* Variable names */
    char **varunit;     /* Variable units */
    Datafile *df;	/* Pointer to the data file */
} TimeSeries;


/* Prototypes */
void    tsRead(char *folderName, char *name, TimeSeries *ts);
void    tsSetGeoType(TimeSeries *ts, char *type);
void    tsSetDefaultGeoType(char *type);
double  tsEval(TimeSeries *ts, int varid, double t);
double  tsEvalEx(TimeSeries *ts, int varid, double t);
double  tsEvalR(TimeSeries *ts, int varid, double t, int rwindid);
double  tsEvalXY(TimeSeries *ts, int varid, double t, double x, double y);
double  tsEvalXYR(TimeSeries *ts, int varid, double t, double x, double y, int rwindid);
double  tsEvalXYZ(TimeSeries *ts, int varid, double t, double x, double y, double z);
double  tsEvalXYZR(TimeSeries *ts, int id, double t, double x, double y, double z, int rwindid);
int  tsMultiCheck(int ntsfiles, TimeSeries *tsfiles, char *var,
			double tstart, double tstop);
double  tsMultiEvalXYZ(int ntsfiles, TimeSeries *tsfiles, char *var,
			 double t, double x, double y, double z);
int     tsIndex(TimeSeries *ts, char *varname);
void    tsNewTimeUnits(TimeSeries *ts, char *newunits);
void    tsFindTime(TimeSeries *ts, double t, int *before, int *after, double *frac);
void    tsFree(TimeSeries *ts);
void    tsCheckTimes(TimeSeries *ts);
void    tsPrintInfo(TimeSeries *ts, FILE *fp);
void    ts_ascii_write(FILE *fp, TimeSeries *ts);


#endif /* _TIMESERIES_H */

