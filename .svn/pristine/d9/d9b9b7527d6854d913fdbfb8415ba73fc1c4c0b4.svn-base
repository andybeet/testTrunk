/*********************************************************************

    File:           ptrack.h
    
    Created:        Mon Jan 15 16:02:42 EST 1996
    
    Author:         Stephen Walker
                    CSIRO Division of Oceanography
    
    Purpose:        Include file for particle tracking routines
    
    Revisions:      none

    $Id: ptrack.h 2761 2011-07-01 04:35:25Z gor171 $

*********************************************************************/

#ifndef SJW_PTRACK
#define SJW_PTRACK 1

typedef struct
{
    double e1;
    double e2;
    double e3;
    short flag;
} Particle;

#define PT_ACTIVE 0x001
#define PT_LOST   0x002

int     createptfile(char *name, long np, char *t_units);
void    readparticles(char *name, int rec, long *np, Particle **p, double *t, char *t_units, int *ndump);
void    writeparticles(int fid, int rec, double t, long np, Particle *p);

#endif /* SJW_PTRACK */
