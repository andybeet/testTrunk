/*

    File:           vgrid.c
    
    Created:        Thu Mar 18 10:32:59 EST 1993
    
    Author:         Stephen Walker
                    CSIRO Division of Oceanography
    
    Purpose:        Calculate vertical grid spacing for anm model
    
    Arguments:      none

    Returns:        void
    
    Revisions:      none

    $Id: vgrid.c 2761 2011-07-01 04:35:25Z gor171 $

*/

#include <stdio.h>
#include <math.h>

int main(void)
{
    int n;
    int nl;
    int k;
    double d1;
    double a;
    double b;
    double c;
    double z;


    fprintf(stderr,"Number of layers > ");
    scanf("%d",&nl);
    n = nl-1;
    fprintf(stderr,"Max depth (+ve)> ");
    scanf("%e",&c);
    c = -fabs(c);
    fprintf(stderr,"Depth of bottom of surface layer (+ve)> ");
    scanf("%e",&d1);
    d1 = fabs(d1);

/*
    a = (d1*(n+1)+c)/(n*n);
    b = d1 - 2.0*a*n;
*/
    a = -c/(n+1) + (d1+c)/n;
    b = -(d1+c)/n - a*n;

    for(k=0; k<=nl; k++) {
	z = a*k*k + b*k + c;
	printf("%d %e\n",k,z);
    }
	return 0;
}



