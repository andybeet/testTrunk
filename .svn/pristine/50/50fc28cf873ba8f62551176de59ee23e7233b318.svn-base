/*********************************************************************

    File:           mapproj.h

    Created:        Thu Aug 27 08:55:52 EST 1998

    Author:         Jason R. Waring
                    CSIRO Marine Research

    Purpose:        Include file for map projection routines.

    Revisions:      none

    $Id: mapproj.h 2761 2011-07-01 04:35:25Z gor171 $

*********************************************************************/

#ifndef MAP_PROJ_H
#define MAP_PROJ_H

/* This structure provides the data and functions required to
 * transform from eastings/northings to latitude/longitude and
 * vice versa, for any number of arbitary projections.
 */
typedef struct MapProjection {
    double ellip_major;		/* Ellipsoidal major axis. */
    double ellip_flat;		/* Ellipsoidal flattening. */
    double falsex;		/* False eastings. */
    double falsey;		/* False northings. */
    void *private_data;		/* Data private to the map projection */
    void *(*init)(struct MapProjection *mp, int nargs, char *args[]);
				/* Projection initialisation */
    void (*free)(void *data);	/* Free the private data */
    void (*forward)(struct MapProjection *mp, double lat, double lon,
			double *east, double *north);
				/* Convert lat/lon to eastings/northings */
    void (*inverse)(struct MapProjection *mp, double east, double north,
			double *lat, double *lon);
				/* Convert easting/northings to lat/lon */
} MapProjection;


/* Prototypes */
MapProjection *mpInitialise(int nargs, char *args[]);
void mpFree(MapProjection *mp);
void mpForward(MapProjection *mp, double  lat, double lon,
                        double *east, double *north);
void mpInverse(MapProjection *mp, double east, double north,
                        double *lat, double *lon);

#endif /* MAP_PROJ_H */
