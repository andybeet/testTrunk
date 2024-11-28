/*********************************************************************

    File:           gridlib.h
    
    Created:        Thu Aug 13 19:47:54 EST 1992
    
    Author:         Stephen Walker
                    CSIRO Division of Oceanography
    
    Purpose:        Include file for library of routines 
		    which calculate grid coordinates
    
    Revisions:      None

    $Id: gridlib.h 2761 2011-07-01 04:35:25Z gor171 $

*********************************************************************/
    
void rect_coord(
    double **x,
    double **y,
    double **h1,
    double **h2,
    double **a1,
    double **a2,
    long nce1,
    long nce2,
    double x00,
    double y00,
    double rotn,
    double xinc,
    double yinc
);
void polar_coord(
    double **x,
    double **y,
    double **h1,
    double **h2,
    double **a1,
    double **a2,
    long nce1,
    long nce2,
    double x00,
    double y00,
    double rotn,
    double arc,
    double rmin
);
void ellipt_coord(
    double **x,
    double **y,
    double **h1,
    double **h2,
    double **a1,
    double **a2,
    long nce1,
    long nce2,
    double x00,
    double y00,
    double rotn,
    double ella,
    double taumax,
    double taumin,
    long nsimm
);
void gridmetric(
    double **x,
    double **y,
    long nce1,
    long nce2,
    double **h1,
    double **h2
);
void gridangle(
    double **x,
    double **y,
    long nce1,
    long nce2,
    double **a1,
    double **a2
);
void CentreToCorner(
    int nce1,
    int nce2,
    double **cx,
    double **cy,
    double **gx,
    double **gy
);

void geog_gridmetric(double **x, double **y, int nce1, int nce2,
                      double **h1, double **h2);
void geog_gridangle(double **x, double **y, int nce1, int nce2,
                      double **a1, double **a2);
