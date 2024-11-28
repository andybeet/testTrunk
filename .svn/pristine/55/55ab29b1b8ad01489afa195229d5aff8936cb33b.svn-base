/*

    File:           waterprops.c
    
    Created:        Fri Jan 5 15:17:32 EST 1996
    
    Author:         Stephen Walker
                    CSIRO Division of Oceanography
    
    Purpose:        Routines to calculate properties of water
    
    Arguments:      See each routine below

    Returns:        void
    
    Revisions:      none

    $Id: waterprops.c 2761 2011-07-01 04:35:25Z gor171 $

*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sjwlib.h>

/*
 * Water vapour pressure.
 * 
 * Data spline interpolated from table on page 107 of
 * Physical Chemistry, Walter J Moore, Longmans 1963
 * 
 */

/* Table of pressures from 0 to 100 degrees C (inclusive)
 * Units here are Pa
 */
double vp[] = 
{
610.4816,
663.1198,
716.3088,
770.5997,
826.5433,
884.6905,
945.5923,
1009.7995,
1077.8630,
1150.3337,
1227.7627,
1310.6425,
1399.2334,
1493.7374,
1594.3566,
1701.2928,
1814.7483,
1934.9249,
2062.0247,
2196.2497,
2337.8020,
2486.9256,
2644.0329,
2809.5784,
2984.0165,
3167.8018,
3361.3886,
3565.2314,
3779.7847,
4005.5029,
4242.8406,
4492.2728,
4754.3571,
5029.6719,
5318.7956,
5622.3064,
5940.7826,
6274.8027,
6624.9448,
6991.7873,
7375.9085,
7777.9258,
8198.6120,
8638.7792,
9099.2394,
9580.8044,
10084.2863,
10610.4971,
11160.2486,
11734.3529,
12333.6219,
12958.8943,
13611.1151,
14291.2563,
15000.2896,
15739.1868,
16508.9199,
17310.4605,
18144.7806,
19012.8519,
19915.6463,
20854.2258,
21830.0132,
22844.5215,
23899.2638,
24995.7531,
26135.5025,
27320.0249,
28550.8334,
29829.4410,
31157.3607,
32535.9668,
33966.0774,
35448.3720,
36983.5300,
38572.2308,
40215.1539,
41912.9786,
43666.3843,
45476.0504,
47342.6564,
49267.6872,
51255.8497,
53312.6565,
55443.6202,
57654.2531,
59950.0679,
62336.5771,
64819.2931,
67403.7285,
70095.3958,
72897.0096,
75800.0935,
78793.3729,
81865.5733,
85005.4205,
88201.6398,
91442.9569,
94718.0973,
98015.7866,
101324.7504
};


/** 
  * Calculates the Water vapour pressure. Uses spline
  * interpolated data from table on page 107 of Physical 
  * Chemistry, Walter J Moore, Longmans 1963.
  * 
  * @param t Temperature at which to evaluate table.
  * @returnwater vapour.
  */
double water_vapourp(double t)
{
    double f;
    int i;

    if( t < 0 || t > 100 )
	warn("water_vapourp: temperature %e outside table range\n",t);
    if( t <= 0 )
	return(vp[0]);
    if( t >= 100 )
	return(vp[100]);
    i = (int)t;
    f = t-i;
    return( (1-f)*vp[i] + f*vp[i+1] );
}

