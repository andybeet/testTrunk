/*********************************************************************

    File:           pointsourcesink.h

    Created:        Tue Jan 2 08:35:49 EST 1996

    Author:         Stephen Walker
                    CSIRO Division of Oceanography

    Purpose:        Include file for point source/sink code

    Revisions:      none

    $Id: pointsourcesink.h 2761 2011-07-01 04:35:25Z gor171 $


	Changes:


	04-05-2010 Bec Gorton
	Added structure to allow each variable in each point source to be scaled
	independently. So there are now scaling arrays stored within the point source structured instead of a single
	array to cover all point source scaling.
*********************************************************************/


#ifndef _POINTSOURCESINK_H
#define _POINTSOURCESINK_H

#include "timeseries.h"

/*********************************************************************
The point source/sink structure
*********************************************************************/
typedef struct
{
    char       name[MAXLINELEN]; /* Name of source/sink */
    double     x;            /* Location x coordinate */
    double     y;            /* Location y coordinate */
    double     z;            /* Location z coordinate */
    double     zlow;         /* Low range of z coordinate */
    double     zhigh;        /* High range of z coordinate */
    TimeSeries *loc;         /* Pointer to location time series */
    int        x_id;         /* Index of x variable in location time series */
    int        y_id;         /* Index of y variable in location time series */
    int        zl_id;        /* Index of z variable in location time series */
    int        zh_id;        /* Index of z variable in location time series */
    int        e1;           /* Integer model grid coordinate */
    int        e2;           /* Integer model grid coordinate */
    int        e3;           /* Integer model grid coordinate */
    int (*xyzijk) (void *, double, double, double, int *, int *, int *);
                             /* Pointer to index routine */
    void       *model_data;  /* Pointer to data needed by index routines */
    TimeSeries ts;          /* Data for this source/sink */
    int       *vid;         /* Box model indices for each time series
			     * variable
			     */
    int        watertsid;   /* Time series index of water variable (ditto) */
	int        rewindid;    /* Flag indicating whether looping of forcing files is on (1=on) */
	int			numPssChanges;
	int 		**pssChange;
	int 		**pssStart;
	int 		**pssPeriod;
	double		**pssMult;
} PointSourceSink;

/* Prototypes */
void    readPointSourceSink(char *folderPath, char *name, char *t_units, PointSourceSink **pss,
	 int *np,  void *data,
	 int (*xyzijk)(void *,double, double, double, int *, int *, int*),
	 int (*trI)(void *, char *));

void    locatePointSourceSink(PointSourceSink *p, double t);

#endif /* _POINTSOURCESINK_H */

