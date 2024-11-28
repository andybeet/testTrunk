/*********************************************************************

    File:           txt_param.h
    
    Created:        Fri Oct 29 09:26:59 EST 1993
    
    Author:         Stephen Walker
                    CSIRO Division of Oceanography
    
    Purpose:        Include file for parameter routines
    
    Revisions:      none

    $Id: txt_param.h 2761 2011-07-01 04:35:25Z gor171 $

*********************************************************************/
    
void	param_fatal(int flag);
long	dparam(FILE *fp, char *label, double *p, double min, double max);
long	fparam(FILE *fp, char *label, float *p, double min, double max);
long	lparam(FILE *fp, char *label, long *p, long min, long max);
long	iparam(FILE *fp, char *label, int *p, int min, int max);
long	sparam(FILE *fp, char *label, short *p, int min, int max);
long	cparam(FILE *fp, char *label, char *p, char *chset);
long	strparam(FILE *fp, char *label, char *p, long maxlen);
long	dparamarray(FILE *fp, char *label, double **p, double min, double max,long *size);
