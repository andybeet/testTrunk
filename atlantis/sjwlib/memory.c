/*

 File:           memory.c

 Created:        Mon May 24 14:52:04 EST 1993

 Author:         Stephen Walker
 CSIRO Division of Oceanography

 Purpose:        Routines for allocation of 1, 2 and 3
 dimensional arrays of various types


 Revisions:      none

 $Id: memory.c 2906 2011-10-11 05:59:06Z gor171 $


 02-02-2009 Bec Gorton
 Added a new checkMemory variable. This is used to check if memory is NULL before its freed.
 If set to TRUE it will not quit with an error message before trying to free null memory - the porgram running
 will throw an exception which can then be used to find the offending array. Unless debugging the model
 shutdown code this should be set to FALSE.

 */

#include <stdio.h>
#include <stdlib.h>
#include "sjwlib.h"
#include <string.h>

int checkMemory = 1;

/** Allocate and clear a 1d array of double values.
 *
 * @param n1 size of array.
 * @return an array n1 of double value.
 */
long double *
d_alloc1longd(long int n1) {
	long double *p;
	size_t size;
	long i;

	if (n1 < 0)
		quit("alloc1d: you have asked to allocate an array of size %d\n", n1);

	/* allocate storage */
	size = (size_t)n1 * sizeof(long double);
	if ((p = (long double *) malloc(size)) == NULL)
		quit("alloc1d: Not enough memory size = %d\n", size);
	for (i = 0; i < n1; i++)
		p[i] = 0.0;
	return ((long double *) p);
}

/** Deallocates the memory for a 1d array of doubles.
 * Was allocated using alloc1d.
 *
 * @param p pointer to aray of doubles.
 */
void d_free1longd(long double *p) {
	if (p == NULL)
		quit("free1d: attempt to free NULL pointer\n");
	free((void *) p);
}

/** Allocate and clear a 2d array. Note that the storage for the values
 * is one large piece of memory so that it can be read or written
 * easily. The i,j value is accessed by saying array[j][i], because
 * array[j] is a pointer to the jth row of the matrix.
 * n1 is the number of columns, n2 is the number of rows.
 *
 * @param n1 size of first dimension of the array.
 * @param n2 size of second dimension of the array.
 * @return a 2d array n1*n2 double values.
 */
long double **
d_alloc2longd(long int n1, long int n2) {
	long double *p;
	long double **pp;
	size_t size;
	long i;

	if (n1 < 0 || n2 < 0)
		quit("alloc2longd: you have asked to allocate an array of size %d, %d\n", n1, n2);

	/* first allocate main storage */
	size = (size_t)n1 * (size_t)n2 * sizeof(long double);
	if ((p = (long double *) malloc(size)) == NULL)
		quit("alloc2d: Not enough memory\n");
	for (i = 0; i < n1 * n2; i++)
		p[i] = 0.0;

	/* now allocate row pointers */
	size = (size_t)n2 * sizeof(long double *);
	if ((pp = (long double **) malloc(size)) == NULL)
		quit("alloc2d: Can't allocate row pointers\n");
	/* point to the rows */
	for (i = 0; i < n2; i++)
		pp[i] = p + i * n1;
	return (pp);
}

/** Deallocates the memory for a 2d array of long double.
 * Was allocated using alloc2d.
 *
 * @param p pointer to aray of long double.
 */
void d_free2longd(long double **pp) {
	long double *p;

	if (checkMemory) {
		if (pp == NULL)
			quit("free2d: attempt to free NULL pointer\n");
	}
	/* extract pointer to main storage */
	p = pp[0];
	/* free row pointer space */
	free((void *) pp);

	/* free main storage */
	if (checkMemory) {
		if (p == NULL)
			quit("free2d: attempt to free NULL pointer\n");
	}
	free((void *) p);
}

/** Allocate a 3d array. Note that the storage for the values
 * is one large piece of memory so that it can be read or written
 * easily. The i,j,k value is accessed by saying array[k][j][i].
 * array[k] is a pointer to the kth plane of the array. array[k][j]
 * is a pointer to the jth row in the kth plane of the array.
 *
 * @param n1 size of first dimension of the array.
 * @param n2 size of second dimension of the array.
 * @param n3 size of third dimension of the array.
 * @return a 3d array n1*n2*n3 long double values.
 */
long double ***
d_alloc3longd(long int n1, long int n2, long int n3) {
	long double *p;
	long double **pp;
	long double ***ppp;
	size_t size;
	long i;

	if (n1 < 0 || n2 < 0 || n3 < 0)
		quit("alloc3longd: you have asked to allocate an array of size %d, %d, %d\n", n1, n2, n3);

	/* first allocate main storage */
	size = (size_t)(n1 * n2 * n3) * sizeof(long double);
	if ((p = (long double *) malloc(size)) == NULL)
		quit("alloc3d: Not enough memory\n");
	for (i = 0; i < n1 * n2 * n3; i++)
		p[i] = 0.0;

	/* now allocate row pointers */
	size = (size_t)(n2 * n3) * sizeof(long double *);
	if ((pp = (long double **) malloc(size)) == NULL)
		quit("alloc3d: Can't allocate row pointers\n");
	/* point to the rows */
	for (i = 0; i < n2 * n3; i++)
		pp[i] = p + i * n1;

	/* now allocate the plane pointers */
	size = (size_t)n3 * sizeof(long double **);
	if ((ppp = (long double ***) malloc(size)) == NULL)
		quit("alloc3d: Can't allocate plane pointers\n");
	/* point to the planes */
	for (i = 0; i < n3; i++)
		ppp[i] = pp + i * n2;
	return (ppp);
}

/** Deallocates the memory for a 3d array of long double.
 * Was allocated using alloc3d.
 *
 * @param p pointer to aray of long double.
 */
void d_free3longd(long double ***ppp) {
	long double **pp;
	long double *p;

	if (checkMemory) {
		if (ppp == NULL)
			quit("free3d: attempt to free NULL pointer\n");
	}
	/* extract pointer to main storage */
	p = ppp[0][0];
	/* extract pointer to row pointer storage */
	pp = ppp[0];
	/* free plane pointer space */
	free((void *) ppp);

	/* free row pointer space */
	if (checkMemory) {
		if (pp == NULL)
			quit("free3d: attempt to free NULL pointer\n");
	}
	free((void *) pp);

	/* free main storage */
	if (checkMemory) {
		if (p == NULL)
			quit("free3d: attempt to free NULL pointer\n");
	}
	free((void *) p);
}

/** Allocate a 4d array. Note that the storage for the values
 * is one large piece of memory so that it can be read or written
 * easily. The i,j,k,l value is accessed by saying array[l][k][j][i].
 *
 * @param n1 size of first dimension of the array.
 * @param n2 size of second dimension of the array.
 * @param n3 size of third dimension of the array.
 * @param n4 size of third dimension of the array.
 * @return a 4d array n1*n2*n3*n4 long double values.
 */
long double ****
d_alloc4longd(long int n1, long int n2, long int n3, long int n4) {
	long double *p;
	long double **pp;
	long double ***ppp;
	long double ****pppp;
	size_t size;
	long i;

	if (n1 < 0 || n2 < 0 || n3 < 0 || n4 < 0)
		quit("alloc4longd: you have asked to allocate an array of size %d, %d, %d, %d\n", n1, n2, n3, n4);

	/* first allocate main storage */
	size = (size_t)(n1 * n2 * n3 * n4) * sizeof(long double);
	if ((p = (long double *) malloc(size)) == NULL)
		quit("alloc4d: Not enough memory\n");
	for (i = 0; i < n1 * n2 * n3 * n4; i++)
		p[i] = 0.0;

	/* now allocate row pointers */
	size = (size_t)(n2 * n3 * n4) * sizeof(long double *);
	if ((pp = (long double **) malloc(size)) == NULL)
		quit("alloc4d: Can't allocate row pointers\n");
	/* point to the rows */
	for (i = 0; i < n2 * n3 * n4; i++)
		pp[i] = p + i * n1;

	/* now allocate the plane pointers */
	size = (size_t)(n3 * n4) * sizeof(long double **);
	if ((ppp = (long double ***) malloc(size)) == NULL)
		quit("alloc4d: Can't allocate plane pointers\n");
	/* point to the planes */
	for (i = 0; i < n3 * n4; i++)
		ppp[i] = pp + i * n2;

	/* now allocate the cube pointers */
	size = (size_t)n4 * sizeof(long double *);
	if ((pppp = (long double ****) malloc(size)) == NULL)
		quit("alloc4d: Can't allocate cube pointers\n");
	/* point to the cubes */
	for (i = 0; i < n4; i++)
		pppp[i] = ppp + i * n3;

	return (pppp);
}

/** Deallocates the memory for a 4d array of double.
 * Was allocated using alloc4d.
 *
 * @param p pointer to aray of double.
 */
void d_free4longd(long double ****pppp) {
	long double ***ppp;
	long double **pp;
	long double *p;

	if (checkMemory) {
		if (pppp == NULL)
			quit("free4d: attempt to free NULL pointer\n");
	}
	/* extract pointer to main storage */
	p = pppp[0][0][0];
	/* extract pointer to row pointer storage */
	pp = pppp[0][0];
	/* extract pointer to plane pointer storage */
	ppp = pppp[0];
	/* free cube pointer space */
	free((void *) pppp);

	/* free plane pointer space */
	if (checkMemory) {
		if (ppp == NULL)
			quit("free4d: attempt to free NULL pointer\n");
	}
	free((void *) ppp);

	/* free row storage */
	if (checkMemory) {
		if (pp == NULL)
			quit("free4d: attempt to free NULL pointer\n");
	}
	free((void *) pp);

	/* free main storage */
	if (checkMemory) {
		if (p == NULL)
			quit("free3d: attempt to free NULL pointer\n");
	}
	free((void *) p);
}


/** Allocate a 5longd array. Note that the storage for the values
 * is one large piece of memory so that it can be read or written
 * easily. The h,i,j,k,l value is accessed by saying array[l][k][j][i][h].
 *
 * @param n1 size of first dimension of the array.
 * @param n2 size of second dimension of the array.
 * @param n3 size of third dimension of the array.
 * @param n4 size of fourth dimension of the array.
 * @param n5 size of fifth dimension of the array
 * @return a 5longd array n1*n2*n3*n4*n5 long double values.
 */
long double *****
d_alloc5longd(long int n1, long int n2, long int n3, long int n4, long int n5) {
	long double *p;
	long double **pp;
	long double ***ppp;
	long double ****pppp;
	long double *****ppppp;
	size_t size;
	long i;

	if (n1 < 0 || n2 < 0 || n3 < 0 || n4 < 0 || n5 < 0)
		quit("alloc5longd: you have asked to allocate an array of size %d, %d, %d, %d, %d\n", n1, n2, n3, n4, n5);

	/* first allocate main storage */
	size = (size_t)(n1 * n2 * n3 * n4 * n5) * sizeof(long double);
	if ((p = (long double *) malloc(size)) == NULL)
		quit("alloc5longd: Not enough memory\n");
	for (i = 0; i < n1 * n2 * n3 * n4 * n5; i++)
		p[i] = 0.0;

	/* now allocate row pointers */
	size = (size_t)(n2 * n3 * n4 * n5) * sizeof(long double *);
	if ((pp = (long double **) malloc(size)) == NULL)
		quit("alloc5longd: Can't allocate row pointers\n");
	/* point to the rows */
	for (i = 0; i < n2 * n3 * n4 * n5; i++)
		pp[i] = p + i * n1;

	/* now allocate the plane pointers */
	size = (size_t)(n3 * n4 * n5) * sizeof(long double **);
	if ((ppp = (long double ***) malloc(size)) == NULL)
		quit("alloc5longd: Can't allocate plane pointers\n");
	/* point to the planes */
	for (i = 0; i < n3 * n4 * n5; i++)
		ppp[i] = pp + i * n2;

	/* now allocate the cube pointers */
	size = (size_t)(n4 * n5) * sizeof(long double *);
	if ((pppp = (long double ****) malloc(size)) == NULL)
		quit("alloc5longd: Can't allocate cube pointers\n");
	/* point to the cubes */
	for (i = 0; i < n4 * n5; i++)
		pppp[i] = ppp + i * n3;

	/* now allocate the hypercube pointers */
	size = (size_t)n5 * sizeof(long double *);
	if ((ppppp = (long double *****) malloc(size)) == NULL)
		quit("alloc5longd: Can't allocate hypercube pointers\n");
	/* point to hypercubes */
	for (i = 0; i < n5; i++)
		ppppp[i] = pppp + i * n4;

	return (ppppp);
}

/** Deallocates the memory for a 5longd array of long doubles.
 * Was allocated using alloc5longd.
 *
 * @param p pointer to aray of long doubles.
 */
void d_free5longd(long double *****ppppp) {
	long double ****pppp;
	long double ***ppp;
	long double **pp;
	long double *p;

	if (checkMemory)
		if (ppppp == NULL)
			quit("free5longd: attempt to free NULL pointer\n");
	/* extract pointer to main storage */
	p = ppppp[0][0][0][0];
	/* extract pointer to row pointer storage */
	pp = ppppp[0][0][0];
	/* extract pointer to plane pointer storage */
	ppp = ppppp[0][0];
	/* free cube pointer space */
	pppp = ppppp[0];
	/* free hypercube pointer space */
	free((void *) ppppp);

	/* free cube pointer space */
	if (checkMemory)
		if (pppp == NULL)
			quit("free5longd: attempt to free NULL pointer\n");
	free((void *) pppp);

	/* free plane pointer space */
	if (checkMemory)
		if (ppp == NULL)
			quit("free5longd: attempt to free NULL pointer\n");
	free((void *) ppp);

	/* free row storage */
	if (checkMemory)
		if (pp == NULL)
			quit("free5longd: attempt to free NULL pointer\n");
	free((void *) pp);

	/* free main storage */
	if (checkMemory)
		if (p == NULL)
			quit("free5longd: attempt to free NULL pointer\n");
	free((void *) p);

	return;
}

/** Allocate a 6longd array. Note that the storage for the values
 * is one large piece of memory so that it can be read or written
 * easily. The h,i,j,k,l,m value is accessed by saying array[m][l][k][j][i][h].
 *
 * @param n1 size of first dimension of the array.
 * @param n2 size of second dimension of the array.
 * @param n3 size of third dimension of the array.
 * @param n4 size of fourth dimension of the array.
 * @param n5 size of fifth dimension of the array
 * @param n6 size of sixth dimension of the array
 * @return a 6longd array n1*n2*n3*n4*n5*n6 long double values.
 */
long double ******
d_alloc6longd(long int n1, long int n2, long int n3, long int n4, long int n5, long int n6) {
	long double *p;
	long double **pp;
	long double ***ppp;
	long double ****pppp;
	long double *****ppppp;
	long double ******pppppp;
	size_t size;
	long i;

	if (n1 < 0 || n2 < 0 || n3 < 0 || n4 < 0 || n5 < 0 ||n6 < 0)
		quit("alloc6longd: you have asked to allocate an array of size %d, %d, %d, %d, %d, %d\n", n1, n2, n3, n4, n5, n6);

	/* first allocate main storage */
	size = (size_t)( n1 * n2 * n3 * n4 * n5 * n6) * sizeof(long double);
	if ((p = (long double *) malloc(size)) == NULL)
		quit("alloc6longd: Not enough memory\n");
	for (i = 0; i < n1 * n2 * n3 * n4 * n5 * n6; i++)
		p[i] = 0.0;

	/* now allocate row pointers */
	size = (size_t)( n2 * n3 * n4 * n5 * n6) * sizeof(long double *);
	if ((pp = (long double **) malloc(size)) == NULL)
		quit("alloc6longd: Can't allocate row pointers\n");
	/* point to the rows */
	for (i = 0; i < n2 * n3 * n4 * n5 * n6; i++)
		pp[i] = p + i * n1;

	/* now allocate the plane pointers */
	size = (size_t)( n3 * n4 * n5 * n6) * sizeof(long double **);
	if ((ppp = (long double ***) malloc(size)) == NULL)
		quit("alloc6longd: Can't allocate plane pointers\n");
	/* point to the planes */
	for (i = 0; i < n3 * n4 * n5 * n6; i++)
		ppp[i] = pp + i * n2;

	/* now allocate the cube pointers */
	size = (size_t)( n4 * n5 * n6) * sizeof(long double ***);
	if ((pppp = (long double ****) malloc(size)) == NULL)
		quit("alloc6longd: Can't allocate cube pointers\n");
	/* point to the cubes */
	for (i = 0; i < n4 * n5 * n6; i++)
		pppp[i] = ppp + i * n3;

	/* now allocate the hypercube pointers */
	size = (size_t)( n5 * n6) * sizeof(long double **);
	if ((ppppp = (long double *****) malloc(size)) == NULL)
		quit("alloc6longd: Can't allocate hypercube pointers\n");
	/* point to hypercubes */
	for (i = 0; i < n5 * n6; i++)
		ppppp[i] = pppp + i * n4;

	/* now allocate the hypercube pointers */
	size = (size_t)( n6) * sizeof(long double *);
	if ((pppppp = (long double ******) malloc(size)) == NULL)
		quit("alloc6longd: Can't allocate higher hypercube pointers\n");
	/* point to hypercubes */
	for (i = 0; i < n6; i++)
		pppppp[i] = ppppp + i * n5;

	return (pppppp);
}

/** Deallocates the memory for a 6longd array of long doubles.
 * Was allocated using alloc6longd.
 *
 * @param p pointer to aray of long doubles.
 */
void d_free6longd(long double ******pppppp) {
	long double *****ppppp;
	long double ****pppp;
	long double ***ppp;
	long double **pp;
	long double *p;

	if (pppppp == NULL)
		quit("free6longd: attempt to free NULL pointer\n");
	/* extract pointer to main storage */
	p = pppppp[0][0][0][0][0];
	/* extract pointer to row pointer storage */
	pp = pppppp[0][0][0][0];
	/* extract pointer to plane pointer storage */
	ppp = pppppp[0][0][0];
	/* free cube pointer space */
	pppp = pppppp[0][0];
	/* free hypercube pointer space */
	ppppp = pppppp[0];
	/* free higher hypercube pointer space */
	free((void *) pppppp);

	/* free hypercube pointer space */
	if (ppppp == NULL)
		quit("free6longd: attempt to free NULL pointer\n");
	free((void *) ppppp);

	/* free cube pointer space */
	if (pppp == NULL)
		quit("free6longd: attempt to free NULL pointer\n");
	free((void *) pppp);

	/* free plane pointer space */
	if (ppp == NULL)
		quit("free6longd: attempt to free NULL pointer\n");
	free((void *) ppp);

	/* free row storage */
	if (pp == NULL)
		quit("free6longd: attempt to free NULL pointer\n");
	free((void *) pp);

	/* free main storage */
	if (p == NULL)
		quit("free6longd: attempt to free NULL pointer\n");
	free((void *) p);

	return;
}
/** Allocate and clear a 1d array of double values.
 *
 * @param n1 size of array.
 * @return an array n1 of double value.
 */
double *
d_alloc1d(long int n1) {
	double *p;
	size_t size;
	long i;

	if (n1 < 0)
		quit("alloc1d: you have asked to allocate an array of size %d\n", n1);

	/* allocate storage */
	size = (size_t)(n1) *sizeof(double);
	if ((p = (double *) malloc(size)) == NULL)
		quit("alloc1d: Not enough memory size = %d\n", size);
	for (i = 0; i < n1; i++)
		p[i] = 0.0;
	return ((double *) p);
}

/** Deallocates the memory for a 1d array of double.
 * Was allocated using alloc1d.
 *
 * @param p pointer to aray of double.
 */
void d_free1d(double *p) {
	if (checkMemory) {
		if (p == NULL)
			quit("free1d: attempt to free NULL pointer\n");
	}
	free((void *) p);
}

/** Allocate and clear a 2d array. Note that the storage for the values
 * is one large piece of memory so that it can be read or written
 * easily. The i,j value is accessed by saying array[j][i], because
 * array[j] is a pointer to the jth row of the matrix.
 * n1 is the number of columns, n2 is the number of rows.
 *
 * @param n1 size of first dimension of the array.
 * @param n2 size of second dimension of the array.
 * @return a 2d array n1*n2 double values.
 */
double **
d_alloc2d(long int n1, long int n2) {
	double *p;
	double **pp;
	size_t size;
	long i;

	if (n1 < 0 || n2 < 0)
		quit("alloc2d: you have asked to allocate an array of size %d, %d\n", n1, n2);

	/* first allocate main storage */
	size = (size_t)(n1 * n2) *sizeof(double);
	if ((p = (double *) malloc(size)) == NULL)
		quit("alloc2d: Not enough memory\n");
	for (i = 0; i < n1 * n2; i++)
		p[i] = 0.0;

	/* now allocate row pointers */
	size = (size_t)(n2) *sizeof(double *);
	if ((pp = (double **) malloc(size)) == NULL)
		quit("alloc2d: Can't allocate row pointers\n");
	/* point to the rows */
	for (i = 0; i < n2; i++)
		pp[i] = p + i * n1;
	return (pp);
}

/** Deallocates the memory for a 2d array of double.
 * Was allocated using alloc2d.
 *
 * @param p pointer to aray of double.
 */
void d_free2d(double **pp) {
	double *p;

	if (checkMemory) {
		if (pp == NULL)
			quit("free2d: attempt to free NULL pointer\n");
	}
	/* extract pointer to main storage */
	p = pp[0];
	/* free row pointer space */
	free((void *) pp);

	/* free main storage */
	if (checkMemory) {
		if (p == NULL)
			quit("free2d: attempt to free NULL pointer\n");
	}
	free((void *) p);
}

/** Allocate a 3d array. Note that the storage for the values
 * is one large piece of memory so that it can be read or written
 * easily. The i,j,k value is accessed by saying array[k][j][i].
 * array[k] is a pointer to the kth plane of the array. array[k][j]
 * is a pointer to the jth row in the kth plane of the array.
 *
 * @param n1 size of first dimension of the array.
 * @param n2 size of second dimension of the array.
 * @param n3 size of third dimension of the array.
 * @return a 3d array n1*n2*n3 double values.
 */
double ***
d_alloc3d(long int n1, long int n2, long int n3) {
	double *p;
	double **pp;
	double ***ppp;
	size_t size;
	long i;

	if (n1 < 0 || n2 < 0 || n3 < 0)
		quit("alloc3d: you have asked to allocate an array of size %d, %d, %d\n", n1, n2, n3);

	/* first allocate main storage */
	size = (size_t)(n1 * n2 * n3) *sizeof(double);
	if ((p = (double *) malloc(size)) == NULL)
		quit("alloc3d: Not enough memory\n");
	for (i = 0; i < n1 * n2 * n3; i++)
		p[i] = 0.0;

	/* now allocate row pointers */
	size = (size_t)(n2 * n3) *sizeof(double *);
	if ((pp = (double **) malloc(size)) == NULL)
		quit("alloc3d: Can't allocate row pointers\n");
	/* point to the rows */
	for (i = 0; i < n2 * n3; i++)
		pp[i] = p + i * n1;

	/* now allocate the plane pointers */
	size = (size_t)(n3) *sizeof(double **);
	if ((ppp = (double ***) malloc(size)) == NULL)
		quit("alloc3d: Can't allocate plane pointers\n");
	/* point to the planes */
	for (i = 0; i < n3; i++)
		ppp[i] = pp + i * n2;
	return (ppp);
}

/** Deallocates the memory for a 3d array of double.
 * Was allocated using alloc3d.
 *
 * @param p pointer to aray of double.
 */
void d_free3d(double ***ppp) {
	double **pp;
	double *p;

	if (checkMemory) {
		if (ppp == NULL)
			quit("free3d: attempt to free NULL pointer\n");
	}
	/* extract pointer to main storage */
	p = ppp[0][0];
	/* extract pointer to row pointer storage */
	pp = ppp[0];
	/* free plane pointer space */
	free((void *) ppp);

	/* free row pointer space */
	if (checkMemory) {
		if (pp == NULL)
			quit("free3d: attempt to free NULL pointer\n");
	}
	free((void *) pp);

	/* free main storage */
	if (checkMemory) {
		if (p == NULL)
			quit("free3d: attempt to free NULL pointer\n");
	}
	free((void *) p);
}

/** Allocate a 4d array. Note that the storage for the values
 * is one large piece of memory so that it can be read or written
 * easily. The i,j,k,l value is accessed by saying array[l][k][j][i].
 *
 * @param n1 size of first dimension of the array.
 * @param n2 size of second dimension of the array.
 * @param n3 size of third dimension of the array.
 * @param n4 size of third dimension of the array.
 * @return a 4d array n1*n2*n3*n4 double values.
 */
double ****
d_alloc4d(long int n1, long int n2, long int n3, long int n4) {
	double *p;
	double **pp;
	double ***ppp;
	double ****pppp;
	size_t size;
	long i;

	if (n1 < 0 || n2 < 0 || n3 < 0 || n4 < 0)
		quit("alloc4d: you have asked to allocate an array of size %d, %d, %d, %d\n", n1, n2, n3, n4);

	/* first allocate main storage */
	size = (size_t)(n1 * n2 * n3 * n4) *sizeof(double);
	if ((p = (double *) malloc(size)) == NULL)
		quit("alloc4d: Not enough memory\n");
	for (i = 0; i < n1 * n2 * n3 * n4; i++)
		p[i] = 0.0;

	/* now allocate row pointers */
	size = (size_t)(n2 * n3 * n4) *sizeof(double *);
	if ((pp = (double **) malloc(size)) == NULL)
		quit("alloc4d: Can't allocate row pointers\n");
	/* point to the rows */
	for (i = 0; i < n2 * n3 * n4; i++)
		pp[i] = p + i * n1;

	/* now allocate the plane pointers */
	size = (size_t)(n3 * n4) *sizeof(double **);
	if ((ppp = (double ***) malloc(size)) == NULL)
		quit("alloc4d: Can't allocate plane pointers\n");
	/* point to the planes */
	for (i = 0; i < n3 * n4; i++)
		ppp[i] = pp + i * n2;

	/* now allocate the cube pointers */
	size = (size_t)(n4) *sizeof(double *);
	if ((pppp = (double ****) malloc(size)) == NULL)
		quit("alloc4d: Can't allocate cube pointers\n");
	/* point to the cubes */
	for (i = 0; i < n4; i++)
		pppp[i] = ppp + i * n3;

	return (pppp);
}

/** Deallocates the memory for a 4d array of double.
 * Was allocated using alloc4d.
 *
 * @param p pointer to aray of double.
 */
void d_free4d(double ****pppp) {
	double ***ppp;
	double **pp;
	double *p;

	if (checkMemory) {
		if (pppp == NULL)
			quit("free4d: attempt to free NULL pointer\n");
	}
	/* extract pointer to main storage */
	p = pppp[0][0][0];
	/* extract pointer to row pointer storage */
	pp = pppp[0][0];
	/* extract pointer to plane pointer storage */
	ppp = pppp[0];
	/* free cube pointer space */
	free((void *) pppp);

	/* free plane pointer space */
	if (checkMemory) {
		if (ppp == NULL)
			quit("free4d: attempt to free NULL pointer\n");
	}
	free((void *) ppp);

	/* free row storage */
	if (checkMemory) {
		if (pp == NULL)
			quit("free4d: attempt to free NULL pointer\n");
	}
	free((void *) pp);

	/* free main storage */
	if (checkMemory) {
		if (p == NULL)
			quit("free4d: attempt to free NULL pointer\n");
	}
	free((void *) p);
}

/** Allocate a 5d array. Note that the storage for the values
 * is one large piece of memory so that it can be read or written
 * easily. The h,i,j,k,l value is accessed by saying array[l][k][j][i][h].
 *
 * @param n1 size of first dimension of the array.
 * @param n2 size of second dimension of the array.
 * @param n3 size of third dimension of the array.
 * @param n4 size of fourth dimension of the array.
 * @param n5 size of fifth dimension of the array
 * @return a 5d array n1*n2*n3*n4*n5 double values.
 */
double *****
d_alloc5d(long int n1, long int n2, long int n3, long int n4, long int n5) {
	double *p;
	double **pp;
	double ***ppp;
	double ****pppp;
	double *****ppppp;
	size_t size;
	long i;

	if (n1 < 0 || n2 < 0 || n3 < 0 || n4 < 0 || n5 < 0)
		quit("alloc5d: you have asked to allocate an array of size %d, %d, %d, %d, %d\n", n1, n2, n3, n4, n5);

	/* first allocate main storage */
	size = (size_t)(n1 * n2 * n3 * n4 * n5) * sizeof(double);
	if ((p = (double *) malloc(size)) == NULL)
		quit("alloc5d: Not enough memory\n");
	for (i = 0; i < n1 * n2 * n3 * n4 * n5; i++)
		p[i] = 0.0;

	/* now allocate row pointers */
	size = (size_t)(n2 * n3 * n4 * n5) * sizeof(double *);
	if ((pp = (double **) malloc(size)) == NULL)
		quit("alloc5d: Can't allocate row pointers\n");
	/* point to the rows */
	for (i = 0; i < n2 * n3 * n4 * n5; i++)
		pp[i] = p + i * n1;

	/* now allocate the plane pointers */
	size = (size_t)(n3 * n4 * n5) * sizeof(double **);
	if ((ppp = (double ***) malloc(size)) == NULL)
		quit("alloc5d: Can't allocate plane pointers\n");
	/* point to the planes */
	for (i = 0; i < n3 * n4 * n5; i++)
		ppp[i] = pp + i * n2;

	/* now allocate the cube pointers */
	size = (size_t)(n4 * n5) * sizeof(double *);
	if ((pppp = (double ****) malloc(size)) == NULL)
		quit("alloc5d: Can't allocate cube pointers\n");
	/* point to the cubes */
	for (i = 0; i < n4 * n5; i++)
		pppp[i] = ppp + i * n3;

	/* now allocate the hypercube pointers */
	size = (size_t)(n5) * sizeof(double *);
	if ((ppppp = (double *****) malloc(size)) == NULL)
		quit("alloc5d: Can't allocate hypercube pointers\n");
	/* point to hypercubes */
	for (i = 0; i < n5; i++)
		ppppp[i] = pppp + i * n4;

	return (ppppp);
}

/** Deallocates the memory for a 5d array of doubles.
 * Was allocated using alloc5d.
 *
 * @param p pointer to aray of doubles.
 */
void d_free5d(double *****ppppp) {
	double ****pppp;
	double ***ppp;
	double **pp;
	double *p;

	if (checkMemory)
		if (ppppp == NULL)
			quit("free5d: attempt to free NULL pointer\n");
	/* extract pointer to main storage */
	p = ppppp[0][0][0][0];
	/* extract pointer to row pointer storage */
	pp = ppppp[0][0][0];
	/* extract pointer to plane pointer storage */
	ppp = ppppp[0][0];
	/* free cube pointer space */
	pppp = ppppp[0];
	/* free hypercube pointer space */
	free((void *) ppppp);

	/* free cube pointer space */
	if (checkMemory)
		if (pppp == NULL)
			quit("free5d: attempt to free NULL pointer\n");
	free((void *) pppp);

	/* free plane pointer space */
	if (checkMemory)
		if (ppp == NULL)
			quit("free5d: attempt to free NULL pointer\n");
	free((void *) ppp);

	/* free row storage */
	if (checkMemory)
		if (pp == NULL)
			quit("free5d: attempt to free NULL pointer\n");
	free((void *) pp);

	/* free main storage */
	if (checkMemory)
		if (p == NULL)
			quit("free5d: attempt to free NULL pointer\n");
	free((void *) p);

	return;
}

/** Allocate a 6d array. Note that the storage for the values
 * is one large piece of memory so that it can be read or written
 * easily. The h,i,j,k,l,m value is accessed by saying array[m][l][k][j][i][h].
 *
 * @param n1 size of first dimension of the array.
 * @param n2 size of second dimension of the array.
 * @param n3 size of third dimension of the array.
 * @param n4 size of fourth dimension of the array.
 * @param n5 size of fifth dimension of the array
 * @param n6 size of sixth dimension of the array
 * @return a 6d array n1*n2*n3*n4*n5*n6 double values.
 */
double ******
d_alloc6d(long int n1, long int n2, long int n3, long int n4, long int n5, long int n6) {
	double *p;
	double **pp;
	double ***ppp;
	double ****pppp;
	double *****ppppp;
	double ******pppppp;
	size_t size;
	long i;

	if (n1 < 0 || n2 < 0 || n3 < 0 || n4 < 0 || n5 < 0 ||n6 < 0)
		quit("alloc6d: you have asked to allocate an array of size %d, %d, %d, %d, %d, %d\n", n1, n2, n3, n4, n5, n6);

	/* first allocate main storage */
	size = (size_t)( n1 * n2 * n3 * n4 * n5 * n6) * sizeof(double);
	if ((p = (double *) malloc(size)) == NULL)
		quit("alloc6d: Not enough memory\n");
	for (i = 0; i < n1 * n2 * n3 * n4 * n5 * n6; i++)
		p[i] = 0.0;

	/* now allocate row pointers */
	size = (size_t)( n2 * n3 * n4 * n5 * n6) * sizeof(double *);
	if ((pp = (double **) malloc(size)) == NULL)
		quit("alloc6d: Can't allocate row pointers\n");
	/* point to the rows */
	for (i = 0; i < n2 * n3 * n4 * n5 * n6; i++)
		pp[i] = p + i * n1;

	/* now allocate the plane pointers */
	size = (size_t)( n3 * n4 * n5 * n6) * sizeof(double **);
	if ((ppp = (double ***) malloc(size)) == NULL)
		quit("alloc6d: Can't allocate plane pointers\n");
	/* point to the planes */
	for (i = 0; i < n3 * n4 * n5 * n6; i++)
		ppp[i] = pp + i * n2;

	/* now allocate the cube pointers */
	size = (size_t)( n4 * n5 * n6) * sizeof(double ***);
	if ((pppp = (double ****) malloc(size)) == NULL)
		quit("alloc6d: Can't allocate cube pointers\n");
	/* point to the cubes */
	for (i = 0; i < n4 * n5 * n6; i++)
		pppp[i] = ppp + i * n3;

	/* now allocate the hypercube pointers */
	size = (size_t)( n5 * n6) * sizeof(double **);
	if ((ppppp = (double *****) malloc(size)) == NULL)
		quit("alloc6d: Can't allocate hypercube pointers\n");
	/* point to hypercubes */
	for (i = 0; i < n5 * n6; i++)
		ppppp[i] = pppp + i * n4;

	/* now allocate the hypercube pointers */
	size = (size_t)( n6) * sizeof(double *);
	if ((pppppp = (double ******) malloc(size)) == NULL)
		quit("alloc6d: Can't allocate higher hypercube pointers\n");
	/* point to hypercubes */
	for (i = 0; i < n6; i++)
		pppppp[i] = ppppp + i * n5;

	return (pppppp);
}

/** Deallocates the memory for a 6d array of doubles.
 * Was allocated using alloc6d.
 *
 * @param p pointer to aray of doubles.
 */
void d_free6d(double ******pppppp) {
	double *****ppppp;
	double ****pppp;
	double ***ppp;
	double **pp;
	double *p;

	if (pppppp == NULL)
		quit("free6d: attempt to free NULL pointer\n");
	/* extract pointer to main storage */
	p = pppppp[0][0][0][0][0];
	/* extract pointer to row pointer storage */
	pp = pppppp[0][0][0][0];
	/* extract pointer to plane pointer storage */
	ppp = pppppp[0][0][0];
	/* free cube pointer space */
	pppp = pppppp[0][0];
	/* free hypercube pointer space */
	ppppp = pppppp[0];
	/* free higher hypercube pointer space */
	free((void *) pppppp);

	/* free hypercube pointer space */
	if (ppppp == NULL)
		quit("free6d: attempt to free NULL pointer\n");
	free((void *) ppppp);

	/* free cube pointer space */
	if (pppp == NULL)
		quit("free6d: attempt to free NULL pointer\n");
	free((void *) pppp);

	/* free plane pointer space */
	if (ppp == NULL)
		quit("free6d: attempt to free NULL pointer\n");
	free((void *) ppp);

	/* free row storage */
	if (pp == NULL)
		quit("free6d: attempt to free NULL pointer\n");
	free((void *) pp);

	/* free main storage */
	if (p == NULL)
		quit("free6d: attempt to free NULL pointer\n");
	free((void *) p);

	return;
}

/** Allocate and clear a 1d array of float values.
 *
 * @param n1 size of array.
 * @return an array n1 of float value.
 */
float *
f_alloc1d(long int n1) {
	float *p;
	size_t size;
	long i;

	if (n1 < 0)
		quit("f_alloc1d: you have asked to allocate an array of size %d\n", n1);

	/* allocate storage */
	size = (size_t)(n1) *sizeof(float);
	if ((p = (float *) malloc(size)) == NULL)
		quit("f_alloc1d: Not enough memory\n");
	for (i = 0; i < n1; i++)
		p[i] = 0;
	return ((float *) p);
}

/** Deallocates the memory for a 1d array of floats.
 * Was allocated using f_alloc1d.
 *
 * @param p pointer to aray of floats.
 */
void f_free1d(float *p) {
	if (checkMemory) {
		if (p == NULL)
			quit("free1d: attempt to free NULL pointer\n");
	}
	free((void *) p);
}

/** Allocate and clear a 2d array. Note that the storage for the values
 * is one large piece of memory so that it can be read or written
 * easily. The i,j value is accessed by saying array[j][i], because
 * array[j] is a pointer to the jth row of the matrix.
 * n1 is the number of columns, n2 is the number of rows.
 *
 * @param n1 size of first dimension of the array.
 * @param n2 size of second dimension of the array.
 * @return a 2d array n1*n2 float values.
 */
float **f_alloc2d(long int n1, long int n2) {
	float *p;
	float **pp;
	size_t size;
	long i;

	if (n1 < 0 || n2 < 0)
		quit("f_alloc2d: you have asked to allocate an array of size %d, %d\n", n1, n2);

	/* first allocate main storage */
	size = (size_t)(n1 * n2) *sizeof(float);
	if ((p = (float *) malloc(size)) == NULL)
		quit("f_alloc2d: Not enough memory\n");
	for (i = 0; i < n1 * n2; i++)
		p[i] = 0;

	/* now allocate row pointers */
	size = (size_t)(n2) *sizeof(float *);
	if ((pp = (float **) malloc(size)) == NULL)
		quit("f_alloc2d: Can't allocate row pointers\n");
	/* point to the rows */
	for (i = 0; i < n2; i++)
		pp[i] = p + i * n1;
	return (pp);
}

/** Deallocates the memory for a 2d array of floats.
 * Was allocated using f_alloc2d.
 *
 * @param p pointer to aray of floats.
 */
void f_free2d(float **pp) {
	float *p;

	if (checkMemory) {
		if (pp == NULL)
			quit("free2d: attempt to free NULL pointer\n");
	}

	/* extract pointer to main storage */
	p = pp[0];
	/* free row pointer space */
	free((void *) pp);

	/* free main storage */
	if (checkMemory) {
		if (p == NULL)
			quit("free2d: attempt to free NULL pointer\n");
	}
	free((void *) p);
}

/** Allocate a 3d array. Note that the storage for the values
 * is one large piece of memory so that it can be read or written
 * easily. The i,j,k value is accessed by saying array[k][j][i].
 * array[k] is a pointer to the kth plane of the array. array[k][j]
 * is a pointer to the jth row in the kth plane of the array.
 *
 * @param n1 size of first dimension of the array.
 * @param n2 size of second dimension of the array.
 * @param n3 size of third dimension of the array.
 * @return a 3d array n1*n2*n3 float values.
 */
float ***
f_alloc3d(long int n1, long int n2, long int n3) {
	float *p;
	float **pp;
	float ***ppp;
	size_t size;
	long i;

	if (n1 < 0 || n2 < 0 || n3 < 0)
		quit("f_alloc3d: you have asked to allocate an array of size %d, %d, %d\n", n1, n2, n3);

	/* first allocate main storage */
	size = (size_t)(n1 * n2 * n3) *sizeof(float);
	if ((p = (float *) malloc(size)) == NULL)
		quit("f_alloc3d: Not enough memory\n");
	for (i = 0; i < n1 * n2 * n3; i++)
		p[i] = 0;

	/* now allocate row pointers */
	size = (size_t)(n2 * n3) *sizeof(float *);
	if ((pp = (float **) malloc(size)) == NULL)
		quit("f_alloc3d: Can't allocate row pointers\n");
	/* point to the rows */
	for (i = 0; i < n2 * n3; i++)
		pp[i] = p + i * n1;

	/* now allocate the plane pointers */
	size = (size_t)(n3) *sizeof(float **);
	if ((ppp = (float ***) malloc(size)) == NULL)
		quit("f_alloc3d: Can't allocate plane pointers\n");
	/* point to the planes */
	for (i = 0; i < n3; i++)
		ppp[i] = pp + i * n2;
	return (ppp);
}

/** Deallocates the memory for a 3d array of floats.
 * Was allocated using f_alloc3d.
 *
 * @param p pointer to aray of floats.
 */
void f_free3d(float ***ppp) {
	float **pp;
	float *p;
	if (checkMemory) {
		if (ppp == NULL)
			quit("free3d: attempt to free NULL pointer\n");
	}
	/* extract pointer to main storage */
	p = ppp[0][0];
	/* extract pointer to row pointer storage */
	pp = ppp[0];
	/* free plane pointer space */
	free((void *) ppp);

	/* free row pointer space */
	if (checkMemory) {
		if (pp == NULL)
			quit("free3d: attempt to free NULL pointer\n");
	}
	free((void *) pp);

	/* free main storage */
	if (checkMemory) {
		if (p == NULL)
			quit("free3d: attempt to free NULL pointer\n");
	}
	free((void *) p);
}

/** Allocate a 4d array. Note that the storage for the values
 * is one large piece of memory so that it can be read or written
 * easily. The i,j,k,l value is accessed by saying array[l][k][j][i].
 *
 * @param n1 size of first dimension of the array.
 * @param n2 size of second dimension of the array.
 * @param n3 size of third dimension of the array.
 * @param n4 size of third dimension of the array.
 * @return a 4d array n1*n2*n3*n4 float values.
 */
float ****
f_alloc4d(long int n1, long int n2, long int n3, long int n4) {
	float *p;
	float **pp;
	float ***ppp;
	float ****pppp;
	size_t size;
	long i;

	if (n1 < 0 || n2 < 0 || n3 < 0 || n4 < 0)
		quit("f_alloc4d: you have asked to allocate an array of size %d, %d, %d, %d\n", n1, n2, n3, n4);

	/* first allocate main storage */
	size = (size_t)(n1 * n2 * n3 * n4) *sizeof(float);
	if ((p = (float *) malloc(size)) == NULL)
		quit("f_alloc4d: Not enough memory\n");
	for (i = 0; i < n1 * n2 * n3 * n4; i++)
		p[i] = 0;

	/* now allocate row pointers */
	size = (size_t)(n2 * n3 * n4) *sizeof(float *);
	if ((pp = (float **) malloc(size)) == NULL)
		quit("f_alloc4d: Can't allocate row pointers\n");
	/* point to the rows */
	for (i = 0; i < n2 * n3 * n4; i++)
		pp[i] = p + i * n1;

	/* now allocate the plane pointers */
	size = (size_t)(n3 * n4) *sizeof(float **);
	if ((ppp = (float ***) malloc(size)) == NULL)
		quit("f_alloc4d: Can't allocate plane pointers\n");
	/* point to the planes */
	for (i = 0; i < n3 * n4; i++)
		ppp[i] = pp + i * n2;

	/* now allocate the cube pointers */
	size = (size_t)(n4) *sizeof(float *);
	if ((pppp = (float ****) malloc(size)) == NULL)
		quit("f_alloc4d: Can't allocate cube pointers\n");
	/* point to the cubes */
	for (i = 0; i < n4; i++)
		pppp[i] = ppp + i * n3;

	return (pppp);
}

/** Deallocates the memory for a 4d array of floats.
 * Was allocated using f_alloc4d.
 *
 * @param p pointer to aray of floats.
 */
void f_free4d(float ****pppp) {
	float ***ppp;
	float **pp;
	float *p;

	if (checkMemory) {
		if (pppp == NULL)
			quit("free4d: attempt to free NULL pointer\n");
	}
	/* extract pointer to main storage */
	p = pppp[0][0][0];
	/* extract pointer to row pointer storage */
	pp = pppp[0][0];
	/* extract pointer to plane pointer storage */
	ppp = pppp[0];
	/* free cube pointer space */
	free((void *) pppp);

	/* free plane pointer space */
	if (checkMemory) {
		if (ppp == NULL)
			quit("free4d: attempt to free NULL pointer\n");
	}
	free((void *) ppp);

	/* free row storage */
	if (checkMemory) {
		if (pp == NULL)
			quit("free4d: attempt to free NULL pointer\n");
	}
	free((void *) pp);

	/* free main storage */
	if (checkMemory) {
		if (p == NULL)
			quit("free4d: attempt to free NULL pointer\n");
	}
	free((void *) p);
}

/** Allocate and clear a 1d array of long values.
 *
 * @param n1 size of array.
 * @return an array n1 of long values.
 */
long *
l_alloc1d(long int n1) {
	long *p;
	size_t size;
	long i;

	if (n1 < 0)
		quit("l_alloc1d: you have asked to allocate an array of size %d\n", n1);

	/* allocate storage */
	size = (size_t)(n1) *sizeof(long);
	if ((p = (long *) malloc(size)) == NULL)
		quit("l_alloc1d: Not enough memory\n");
	for (i = 0; i < n1; i++)
		p[i] = 0;
	return ((long *) p);
}

/** Deallocates the memory for a 1d array of longs.
 * Was allocated using l_alloc1d.
 *
 * @param p pointer to aray of longs.
 */
void l_free1d(long int *p) {
	if (checkMemory) {
		if (p == NULL)
			quit("free1d: attempt to free NULL pointer\n");
	}
	free((void *) p);
}

/** Allocate and clear a 2d array. Note that the storage for the values
 * is one large piece of memory so that it can be read or written
 * easily. The i,j value is accessed by saying array[j][i], because
 * array[j] is a pointer to the jth row of the matrix.
 * n1 is the number of columns, n2 is the number of rows.
 *
 * @param n1 size of first dimension of the array.
 * @param n2 size of second dimension of the array.
 * @return a 2d array n1*n2 long values.
 */
long **
l_alloc2d(long int n1, long int n2) {
	long *p;
	long **pp;
	size_t size;
	long i;

	if (n1 < 0 || n2 < 0)
		quit("l_alloc2d: you have asked to allocate an array of size %d, %d\n", n1, n2);

	/* first allocate main storage */
	size = (size_t)(n1 * n2) *sizeof(long);
	if ((p = (long *) malloc(size)) == NULL)
		quit("l_alloc2d: Not enough memory\n");
	for (i = 0; i < n1 * n2; i++)
		p[i] = 0;

	/* now allocate row pointers */
	size = (size_t)(n2) *sizeof(long *);
	if ((pp = (long **) malloc(size)) == NULL)
		quit("l_alloc2d: Can't allocate row pointers\n");
	/* point to the rows */
	for (i = 0; i < n2; i++)
		pp[i] = p + i * n1;
	return (pp);
}

/** Deallocates the memory for a 2d array of longs.
 * Was allocated using l_alloc2d.
 *
 * @param p pointer to aray of longs.
 */
void l_free2d(long int **pp) {
	long *p;

	if (checkMemory) {
		if (pp == NULL)
			quit("free2d: attempt to free NULL pointer\n");
	}
	/* extract pointer to main storage */
	p = pp[0];
	/* free row pointer space */
	free((void *) pp);

	/* free main storage */
	if (checkMemory) {
		if (p == NULL)
			quit("free2d: attempt to free NULL pointer\n");
	}
	free((void *) p);
}

/** Allocate a 3d array. Note that the storage for the values
 * is one large piece of memory so that it can be read or written
 * easily. The i,j,k value is accessed by saying array[k][j][i].
 * array[k] is a pointer to the kth plane of the array. array[k][j]
 * is a pointer to the jth row in the kth plane of the array.
 *
 * @param n1 size of first dimension of the array.
 * @param n2 size of second dimension of the array.
 * @param n3 size of third dimension of the array.
 * @return a 3d array n1*n2*n3 long values.
 */
long ***
l_alloc3d(long int n1, long int n2, long int n3) {
	long *p;
	long **pp;
	long ***ppp;
	size_t size;
	long i;

	if (n1 < 0 || n2 < 0 || n3 < 0)
		quit("l_alloc3d: you have asked to allocate an array of size %d, %d, %d\n", n1, n2, n3);

	/* first allocate main storage */
	size = (size_t)(n1 * n2 * n3) *sizeof(long);
	if ((p = (long *) malloc(size)) == NULL)
		quit("l_alloc3d: Not enough memory\n");
	for (i = 0; i < n1 * n2 * n3; i++)
		p[i] = 0;

	/* now allocate row pointers */
	size = (size_t)(n2 * n3) *sizeof(long *);
	if ((pp = (long **) malloc(size)) == NULL)
		quit("l_alloc3d: Can't allocate row pointers\n");
	/* point to the rows */
	for (i = 0; i < n2 * n3; i++)
		pp[i] = p + i * n1;

	/* now allocate the plane pointers */
	size = (size_t)(n3) *sizeof(long **);
	if ((ppp = (long ***) malloc(size)) == NULL)
		quit("l_alloc3d: Can't allocate plane pointers\n");
	/* point to the planes */
	for (i = 0; i < n3; i++)
		ppp[i] = pp + i * n2;
	return (ppp);
}

/** Deallocates the memory for a 3d array of longs.
 * Was allocated using l_alloc3d.
 *
 * @param p pointer to aray of longs.
 */
void l_free3d(long int ***ppp) {
	long **pp;
	long *p;

	if (checkMemory) {
		if (ppp == NULL)
			quit("free3d: attempt to free NULL pointer\n");
	}
	/* extract pointer to main storage */
	p = ppp[0][0];
	/* extract pointer to row pointer storage */
	pp = ppp[0];
	/* free plane pointer space */
	free((void *) ppp);

	/* free row pointer space */
	if (checkMemory) {
		if (pp == NULL)
			quit("free3d: attempt to free NULL pointer\n");
	}
	free((void *) pp);

	/* free main storage */
	if (checkMemory) {
		if (p == NULL)
			quit("free3d: attempt to free NULL pointer\n");
	}
	free((void *) p);
}

/** Allocate and clear a 1d array of integer values.
 *
 * @param n1 size of array.
 * @return an array n1 of integer values.
 */
int *
i_alloc1d(long int n1) {
	int *p;
	size_t size;
	long i;

	if (n1 < 0)
		quit("i_alloc1d: you have asked to allocate an array of size %d\n", n1);

	/* allocate storage */
	size = (size_t)(n1) *sizeof(int);
	if ((p = (int *) malloc(size)) == NULL)
		quit("i_alloc1d: Not enough memory\n");
	for (i = 0; i < n1; i++)
		p[i] = 0;
	return ((int *) p);
}

/** Deallocates the memory for a 1d array of integers.
 * Was allocated using i_alloc1d.
 *
 * @param p pointer to aray of integers.
 */
/* free a 1d array */
void i_free1d(int *p) {
	if (checkMemory) {
		if (p == NULL)
			quit("free1d: attempt to free NULL pointer\n");
	}
	free((void *) p);
}

/** Allocate and clear a 2d array. Note that the storage for the values
 * is one large piece of memory so that it can be read or written
 * easily. The i,j value is accessed by saying array[j][i], because
 * array[j] is a pointer to the jth row of the matrix.
 * n1 is the number of columns, n2 is the number of rows.
 *
 * @param n1 size of first dimension of the array.
 * @param n2 size of second dimension of the array.
 * @return a 2d array n1*n2 integer values.
 */
int **i_alloc2d(long int n1, long int n2) {
	int *p;
	int **pp;
	size_t size;
	long i;

	if (n1 < 0 || n2 < 0)
		quit("i_alloc2d: you have asked to allocate an array of size %d, %d\n", n1, n2);

	/* first allocate main storage */
	size = (size_t)(n1 * n2) *sizeof(int);
	if ((p = (int *) malloc(size)) == NULL)
		quit("i_alloc2d: Not enough memory\n");
	for (i = 0; i < n1 * n2; i++)
		p[i] = 0;

	/* now allocate row pointers */
	size = (size_t)(n2) *sizeof(int *);
	if ((pp = (int **) malloc(size)) == NULL)
		quit("i_alloc2d: Can't allocate row pointers\n");
	/* point to the rows */
	for (i = 0; i < n2; i++)
		pp[i] = p + i * n1;
	return (pp);
}

/** Deallocates the memory for a 2d array of integers.
 * Was allocated using i_alloc2d.
 *
 * @param p pointer to aray of integers.
 */
void i_free2d(int **pp) {
	int *p;
	if (checkMemory) {
		if (pp == NULL)
			quit("free2d: attempt to free NULL pointer\n");
	}
	/* extract pointer to main storage */
	p = pp[0];
	/* free row pointer space */
	free((void *) pp);

	/* free main storage */
	if (checkMemory) {
		if (p == NULL)
			quit("free2d: attempt to free NULL pointer\n");
	}
	free((void *) p);
}

/** Allocate a 3d array. Note that the storage for the values
 * is one large piece of memory so that it can be read or written
 * easily. The i,j,k value is accessed by saying array[k][j][i].
 * array[k] is a pointer to the kth plane of the array. array[k][j]
 * is a pointer to the jth row in the kth plane of the array.
 *
 * @param n1 size of first dimension of the array.
 * @param n2 size of second dimension of the array.
 * @param n3 size of third dimension of the array.
 * @return a 3d array n1*n2*n3 integer values.
 */
int ***
i_alloc3d(long int n1, long int n2, long int n3) {
	int *p;
	int **pp;
	int ***ppp;
	size_t size;
	long i;

	if (n1 < 0 || n2 < 0 || n3 < 0)
		quit("i_alloc3d: you have asked to allocate an array of size %d, %d, %d\n", n1, n2, n3);

	/* first allocate main storage */
	size = (size_t)(n1 * n2 * n3) *sizeof(int);
	if ((p = (int *) malloc(size)) == NULL)
		quit("i_alloc3d: Not enough memory\n");
	for (i = 0; i < n1 * n2 * n3; i++)
		p[i] = 0;

	/* now allocate row pointers */
	size = (size_t)(n2 * n3) *sizeof(int *);
	if ((pp = (int **) malloc(size)) == NULL)
		quit("i_alloc3d: Can't allocate row pointers\n");
	/* point to the rows */
	for (i = 0; i < n2 * n3; i++)
		pp[i] = p + i * n1;

	/* now allocate the plane pointers */
	size = (size_t)(n3) *sizeof(int **);
	if ((ppp = (int ***) malloc(size)) == NULL)
		quit("alloc3d: Can't allocate plane pointers\n");
	/* point to the planes */
	for (i = 0; i < n3; i++)
		ppp[i] = pp + i * n2;
	return (ppp);
}

/** Deallocates the memory for a 3d array of integers.
 * Was allocated using i_alloc3d.
 *
 * @param p pointer to aray of integers.
 */
void i_free3d(int ***ppp) {
	int **pp;
	int *p;

	if (checkMemory) {
		if (ppp == NULL)
			quit("free3d: attempt to free NULL pointer\n");
	}
	/* extract pointer to main storage */
	p = ppp[0][0];
	/* extract pointer to row pointer storage */
	pp = ppp[0];
	/* free plane pointer space */
	free((void *) ppp);

	/* free row pointer space */
	if (checkMemory) {
		if (pp == NULL)
			quit("free3d: attempt to free NULL pointer\n");
	}
	free((void *) pp);

	/* free main storage */
	if (checkMemory) {
		if (p == NULL)
			quit("free3d: attempt to free NULL pointer\n");
	}
	free((void *) p);
}

/** Allocate a 4d array. Note that the storage for the values
 * is one large piece of memory so that it can be read or written
 * easily. The i,j,k,l value is accessed by saying array[l][k][j][i].
 *
 * @param n1 size of first dimension of the array.
 * @param n2 size of second dimension of the array.
 * @param n3 size of third dimension of the array.
 * @param n4 size of third dimension of the array.
 * @return a 4d array n1*n2*n3*n4 integer values.
 */
int ****
i_alloc4d(long int n1, long int n2, long int n3, long int n4) {
	int *p;
	int **pp;
	int ***ppp;
	int ****pppp;
	size_t size;
	long i;

	if (n1 < 0 || n2 < 0 || n3 < 0 || n4 < 0)
		quit("i_alloc4d: you have asked to allocate an array of size %d, %d, %d, %d\n", n1, n2, n3, n4);

	/* first allocate main storage */
	size = (size_t)(n1 * n2 * n3 * n4) *sizeof(int);
	if ((p = (int *) malloc(size)) == NULL)
		quit("i_alloc4d: Not enough memory\n");
	for (i = 0; i < n1 * n2 * n3 * n4; i++)
		p[i] = 0;

	/* now allocate row pointers */
	size = (size_t)(n2 * n3 * n4) *sizeof(int *);
	if ((pp = (int **) malloc(size)) == NULL)
		quit("i_alloc4d: Can't allocate row pointers\n");
	/* point to the rows */
	for (i = 0; i < n2 * n3 * n4; i++)
		pp[i] = p + i * n1;

	/* now allocate the plane pointers */
	size = (size_t)(n3 * n4) *sizeof(int **);
	if ((ppp = (int ***) malloc(size)) == NULL)
		quit("i_alloc4d: Can't allocate plane pointers\n");
	/* point to the planes */
	for (i = 0; i < n3 * n4; i++)
		ppp[i] = pp + i * n2;

	/* now allocate the cube pointers */
	size = (size_t)(n4) *sizeof(int *);
	if ((pppp = (int ****) malloc(size)) == NULL)
		quit("i_alloc4d: Can't allocate cube pointers\n");
	/* point to the cubes */
	for (i = 0; i < n4; i++)
		pppp[i] = ppp + i * n3;

	return (pppp);
}

/** Deallocates the memory for a 4d array of integers.
 * Was allocated using i_alloc4d.
 *
 * @param p pointer to aray of integers.
 */
void i_free4d(int ****pppp) {
	int ***ppp;
	int **pp;
	int *p;

	if (checkMemory) {
		if (pppp == NULL)
			quit("free4d: attempt to free NULL pointer\n");
	}
	/* extract pointer to main storage */
	p = pppp[0][0][0];
	/* extract pointer to row pointer storage */
	pp = pppp[0][0];
	/* extract pointer to plane pointer storage */
	ppp = pppp[0];
	/* free cube pointer space */
	free((void *) pppp);

	/* free plane pointer space */
	if (checkMemory) {
		if (ppp == NULL)
			quit("free4d: attempt to free NULL pointer\n");
	}
	free((void *) ppp);

	/* free row storage */
	if (checkMemory) {
		if (pp == NULL)
			quit("free4d: attempt to free NULL pointer\n");
	}
	free((void *) pp);

	/* free main storage */
	if (checkMemory) {
		if (p == NULL)
			quit("free3d: attempt to free NULL pointer\n");
	}

	free((void *) p);
}

/** Allocate and clear a 1d array of short values.
 *
 * @param n1 size of array.
 * @return an array n1 of short values.
 */
short *
s_alloc1d(long int n1) {
	short *p;
	size_t size;
	long i;

	if (n1 < 0)
		quit("s_alloc1d: you have asked to allocate an array of size %d\n", n1);

	/* allocate storage */
	size = (size_t)(n1) *sizeof(short);
	if ((p = (short *) malloc(size)) == NULL)
		quit("s_alloc1d: Not enough memory\n");
	for (i = 0; i < n1; i++)
		p[i] = 0;
	return ((short *) p);
}

/** Deallocates the memory for a 1d array of shorts.
 * Was allocated using s_alloc1d.
 *
 * @param p pointer to aray of shorts.
 */
void s_free1d(short int *p) {
	if (p == NULL)
		quit("s_free1d: attempt to free NULL pointer\n");
	free((void *) p);
}

/** Allocate and clear a 2d array. Note that the storage for the values
 * is one large piece of memory so that it can be read or written
 * easily. The i,j value is accessed by saying array[j][i], because
 * array[j] is a pointer to the jth row of the matrix.
 * n1 is the number of columns, n2 is the number of rows.
 *
 * @param n1 size of first dimension of the array.
 * @param n2 size of second dimension of the array.
 * @return a 2d array n1*n2 short values.
 */
short **
s_alloc2d(long int n1, long int n2) {
	short *p;
	short **pp;
	size_t size;
	long i;

	if (n1 < 0 || n2 < 0)
		quit("s_alloc2d: you have asked to allocate an array of size %d, %d\n", n1, n2);

	/* first allocate main storage */
	size = (size_t)(n1 * n2) *sizeof(short);
	if ((p = (short *) malloc(size)) == NULL)
		quit("s_alloc2d: Not enough memory\n");
	for (i = 0; i < n1 * n2; i++)
		p[i] = 0;

	/* now allocate row pointers */
	size = (size_t)(n2) *sizeof(short *);
	if ((pp = (short **) malloc(size)) == NULL)
		quit("s_alloc2d: Can't allocate row pointers\n");
	/* point to the rows */
	for (i = 0; i < n2; i++)
		pp[i] = p + i * n1;
	return (pp);
}

/** Deallocates the memory for a 2d array of shorts.
 * Was allocated using s_alloc2d.
 *
 * @param p pointer to aray of shorts.
 */
void s_free2d(short int **pp) {
	short *p;
	if (checkMemory) {
		if (pp == NULL)
			quit("free2d: attempt to free NULL pointer\n");
	}
	/* extract pointer to main storage */
	p = pp[0];
	/* free row pointer space */
	free((void *) pp);

	/* free main storage */
	if (checkMemory) {
		if (p == NULL)
			quit("free2d: attempt to free NULL pointer\n");
	}
	free((void *) p);
}

/** Allocate a 3d array. Note that the storage for the values
 * is one large piece of memory so that it can be read or written
 * easily. The i,j,k value is accessed by saying array[k][j][i].
 * array[k] is a pointer to the kth plane of the array. array[k][j]
 * is a pointer to the jth row in the kth plane of the array.
 *
 * @param n1 size of first dimension of the array.
 * @param n2 size of second dimension of the array.
 * @param n3 size of third dimension of the array.
 * @return a 3d array n1*n2*n3 short values.
 */
short ***
s_alloc3d(long int n1, long int n2, long int n3) {
	short *p;
	short **pp;
	short ***ppp;
	size_t size;
	long i;

	if (n1 < 0 || n2 < 0 || n3 < 0)
		quit("s_alloc3d: you have asked to allocate an array of size %d, %d, %d\n", n1, n2, n3);

	/* first allocate main storage */
	size = (size_t)(n1 * n2 * n3) *sizeof(short);
	if ((p = (short *) malloc(size)) == NULL)
		quit("s_alloc3d: Not enough memory\n");
	for (i = 0; i < n1 * n2 * n3; i++)
		p[i] = 0;

	/* now allocate row pointers */
	size = (size_t)(n2 * n3) *sizeof(short *);
	if ((pp = (short **) malloc(size)) == NULL)
		quit("s_alloc3d: Can't allocate row pointers\n");
	/* point to the rows */
	for (i = 0; i < n2 * n3; i++)
		pp[i] = p + i * n1;

	/* now allocate the plane pointers */
	size = (size_t)(n3) *sizeof(short **);
	if ((ppp = (short ***) malloc(size)) == NULL)
		quit("s_alloc3d: Can't allocate plane pointers\n");
	/* point to the planes */
	for (i = 0; i < n3; i++)
		ppp[i] = pp + i * n2;
	return (ppp);
}

/** Deallocates the memory for a 3d array of shorts.
 * Was allocated using s_alloc3d.
 *
 * @param p pointer to aray of shorts.
 */
void s_free3d(short int ***ppp) {
	short **pp;
	short *p;

	if (ppp == NULL)
		quit("s_free3d: attempt to free NULL pointer\n");
	/* extract pointer to main storage */
	p = ppp[0][0];
	/* extract pointer to row pointer storage */
	pp = ppp[0];
	/* free plane pointer space */
	free((void *) ppp);

	/* free row pointer space */
	if (pp == NULL)
		quit("s_free3d: attempt to free NULL pointer\n");
	free((void *) pp);

	/* free main storage */
	if (p == NULL)
		quit("s_free3d: attempt to free NULL pointer\n");
	free((void *) p);
}

/** Allocate a 4d array. Note that the storage for the values
 * is one large piece of memory so that it can be read or written
 * easily. The i,j,k,l value is accessed by saying array[l][k][j][i].
 *
 * @param n1 size of first dimension of the array.
 * @param n2 size of second dimension of the array.
 * @param n3 size of third dimension of the array.
 * @param n4 size of third dimension of the array.
 * @return a 4d array n1*n2*n3*n4 short values.
 */
short ****
s_alloc4d(long int n1, long int n2, long int n3, long int n4) {
	short *p;
	short **pp;
	short ***ppp;
	short ****pppp;
	size_t size;
	long i;

	if (n1 < 0 || n2 < 0 || n3 < 0 || n4 < 0)
		quit("s_alloc4d: you have asked to allocate an array of size %d, %d, %d, %d\n", n1, n2, n3, n4);

	/* first allocate main storage */
	size = (size_t)(n1 * n2 * n3 * n4) *sizeof(short);
	if ((p = (short *) malloc(size)) == NULL)
		quit("s_alloc4d: Not enough memory\n");
	for (i = 0; i < n1 * n2 * n3 * n4; i++)
		p[i] = 0;

	/* now allocate row pointers */
	size = (size_t)(n2 * n3 * n4) *sizeof(short *);
	if ((pp = (short **) malloc(size)) == NULL)
		quit("s_alloc4d: Can't allocate row pointers\n");
	/* point to the rows */
	for (i = 0; i < n2 * n3 * n4; i++)
		pp[i] = p + i * n1;

	/* now allocate the plane pointers */
	size = (size_t)(n3 * n4) *sizeof(short **);
	if ((ppp = (short ***) malloc(size)) == NULL)
		quit("s_alloc4d: Can't allocate plane pointers\n");
	/* point to the planes */
	for (i = 0; i < n3 * n4; i++)
		ppp[i] = pp + i * n2;

	/* now allocate the cube pointers */
	size = (size_t)(n4) *sizeof(short *);
	if ((pppp = (short ****) malloc(size)) == NULL)
		quit("s_alloc4d: Can't allocate cube pointers\n");
	/* point to the cubes */
	for (i = 0; i < n4; i++)
		pppp[i] = ppp + i * n3;

	return (pppp);
}

/** Deallocates the memory for a 4d array of shorts.
 * Was allocated using s_alloc4d.
 *
 * @param p pointer to aray of shorts.
 */
void s_free4d(short int ****pppp) {
	short ***ppp;
	short **pp;
	short *p;

	if (pppp == NULL)
		quit("free4d: attempt to free NULL pointer\n");
	/* extract pointer to main storage */
	p = pppp[0][0][0];
	/* extract pointer to row pointer storage */
	pp = pppp[0][0];
	/* extract pointer to plane pointer storage */
	ppp = pppp[0];
	/* free cube pointer space */
	free((void *) pppp);

	/* free plane pointer space */
	if (ppp == NULL)
		quit("free4d: attempt to free NULL pointer\n");
	free((void *) ppp);

	/* free row storage */
	if (pp == NULL)
		quit("free4d: attempt to free NULL pointer\n");
	free((void *) pp);

	/* free main storage */
	if (p == NULL)
		quit("free3d: attempt to free NULL pointer\n");
	free((void *) p);
}

/** Allocate and clear a 1d array of char values.
 *
 * @param n1 size of array.
 *
 * @return an array n1 of char values.
 */
char *
c_alloc1d(long int n1) {
	char *p;
	size_t size;
	long i;

	if (n1 < 0)
		quit("c_alloc1d: you have asked to allocate an array of size %d\n", n1);

	/* allocate storage */
	size = (size_t)(n1) *sizeof(char);
	if ((p = (char *) malloc(size)) == NULL)
		quit("c_alloc1d: Not enough memory\n");
	for (i = 0; i < n1; i++)
		p[i] = 0;
	return ((char *) p);
}

/** Deallocates the memory for a 1d array of chars.
 * Was allocated using c_alloc1d.
 *
 * @param p pointer to aray of chars.
 */
void c_free1d(char *p) {
	if (p == NULL)
		quit("free1d: attempt to free NULL pointer\n");
	free((void *) p);
}

/** Allocate and clear a 2d array. Note that the storage for the values
 * is one large piece of memory so that it can be read or written
 * easily. The i,j value is accessed by saying array[j][i], because
 * array[j] is a pointer to the jth row of the matrix.
 * n1 is the number of columns, n2 is the number of rows.
 *
 * @param n1 size of first dimension of the array.
 * @param n2 size of second dimension of the array.
 * @return a 2d array n1*n2 char values.
 */
char **c_alloc2d(long int n1, long int n2) {
	char *p;
	char **pp;
	size_t size;
	long i;

	if (n1 < 0 || n2 < 0)
		quit("c_alloc2d: you have asked to allocate an array of size %d, %d\n", n1, n2);

	/* first allocate main storage */
	size = (size_t)(n1 * n2) *sizeof(char);
	if ((p = (char *) malloc(size)) == NULL)
		quit("c_alloc2d: Not enough memory\n");
	for (i = 0; i < n1 * n2; i++)
		p[i] = 0;

	/* now allocate row pointers */
	size = (size_t)(n2) *sizeof(char *);
	if ((pp = (char **) malloc(size)) == NULL)
		quit("c_alloc2d: Can't allocate row pointers\n");
	/* point to the rows */
	for (i = 0; i < n2; i++)
		pp[i] = p + i * n1;
	return (pp);
}

/** Deallocates the memory for a 2d array of chars.
 * Was allocated using c_alloc2d.
 *
 * @param p pointer to aray of chars.
 */
void c_free2d(char **pp) {
	char *p;
	if (checkMemory) {
		if (pp == NULL)
			quit("free2d: attempt to free NULL pointer\n");
	}
	/* extract pointer to main storage */
	p = pp[0];
	/* free row pointer space */
	free((void *) pp);

	/* free main storage */
	if (checkMemory) {
		if (p == NULL)
			quit("free2d: attempt to free NULL pointer\n");
	}
	free((void *) p);
}

/** Allocate a 3d array. Note that the storage for the values
 * is one large piece of memory so that it can be read or written
 * easily. The i,j,k value is accessed by saying array[k][j][i].
 * array[k] is a pointer to the kth plane of the array. array[k][j]
 * is a pointer to the jth row in the kth plane of the array.
 *
 * @param n1 size of first dimension of the array.
 * @param n2 size of second dimension of the array.
 * @param n3 size of third dimension of the array.
 * @return a 3d array n1*n2*n3 char values.
 */
char ***
c_alloc3d(long int n1, long int n2, long int n3) {
	char *p;
	char **pp;
	char ***ppp;
	size_t size;
	long i;

	if (n1 < 0 || n2 < 0 || n3 < 0)
		quit("c_alloc3d: you have asked to allocate an array of size %d, %d, %d\n", n1, n2, n3);

	/* first allocate main storage */
	size = (size_t)(n1 * n2 * n3) *sizeof(char);
	if ((p = (char *) malloc(size)) == NULL)
		quit("c_alloc3d: Not enough memory\n");
	for (i = 0; i < n1 * n2 * n3; i++)
		p[i] = 0;

	/* now allocate row pointers */
	size = (size_t)(n2 * n3) *sizeof(char *);
	if ((pp = (char **) malloc(size)) == NULL)
		quit("c_alloc3d: Can't allocate row pointers\n");
	/* point to the rows */
	for (i = 0; i < n2 * n3; i++)
		pp[i] = p + i * n1;

	/* now allocate the plane pointers */
	size = (size_t)(n3) *sizeof(char **);
	if ((ppp = (char ***) malloc(size)) == NULL)
		quit("c_alloc3d: Can't allocate plane pointers\n");
	/* point to the planes */
	for (i = 0; i < n3; i++)
		ppp[i] = pp + i * n2;
	return (ppp);
}

/** Deallocates the memory for a 3d array of chars.
 * Was allocated using c_alloc3d.
 *
 * @param p pointer to aray of chars.
 */
void c_free3d(char ***ppp) {
	char **pp;
	char *p;

	if (ppp == NULL)
		quit("free3d: attempt to free NULL pointer\n");
	/* extract pointer to main storage */
	p = ppp[0][0];
	/* extract pointer to row pointer storage */
	pp = ppp[0];
	/* free plane pointer space */
	free((void *) ppp);

	/* free row pointer space */
	if (pp == NULL)
		quit("free3d: attempt to free NULL pointer\n");
	free((void *) pp);

	/* free main storage */
	if (p == NULL)
		quit("free3d: attempt to free NULL pointer\n");
	free((void *) p);
}

/** Allocate a 4d array. Note that the storage for the values
 * is one large piece of memory so that it can be read or written
 * easily. The i,j,k,l value is accessed by saying array[l][k][j][i].
 *
 * @param n1 size of first dimension of the array.
 * @param n2 size of second dimension of the array.
 * @param n3 size of third dimension of the array.
 * @param n4 size of third dimension of the array.
 * @return a 4d array n1*n2*n3*n4 char values.
 */
char ****
c_alloc4d(long int n1, long int n2, long int n3, long int n4) {
	char *p;
	char **pp;
	char ***ppp;
	char ****pppp;
	size_t size;
	long i;

	if (n1 < 0 || n2 < 0 || n3 < 0 || n4 < 0)
		quit("c_alloc4d: you have asked to allocate an array of size %d, %d, %d, %d\n", n1, n2, n3, n4);

	/* first allocate main storage */
	size = (size_t)(n1 * n2 * n3 * n4) *sizeof(char);
	if ((p = (char *) malloc(size)) == NULL)
		quit("c_alloc4d: Not enough memory\n");
	for (i = 0; i < n1 * n2 * n3 * n4; i++)
		p[i] = 0;

	/* now allocate row pointers */
	size = (size_t)(n2 * n3 * n4) *sizeof(char *);
	if ((pp = (char **) malloc(size)) == NULL)
		quit("c_alloc4d: Can't allocate row pointers\n");
	/* point to the rows */
	for (i = 0; i < n2 * n3 * n4; i++)
		pp[i] = p + i * n1;

	/* now allocate the plane pointers */
	size = (size_t)(n3 * n4) *sizeof(char **);
	if ((ppp = (char ***) malloc(size)) == NULL)
		quit("c_alloc4d: Can't allocate plane pointers\n");
	/* point to the planes */
	for (i = 0; i < n3 * n4; i++)
		ppp[i] = pp + i * n2;

	/* now allocate the cube pointers */
	size = (size_t)(n4) *sizeof(char *);
	if ((pppp = (char ****) malloc(size)) == NULL)
		quit("c_alloc4d: Can't allocate cube pointers\n");
	/* point to the cubes */
	for (i = 0; i < n4; i++)
		pppp[i] = ppp + i * n3;

	return (pppp);
}

/** Deallocates the memory for a 4d array of chars.
 * Was allocated using c_alloc4d.
 *
 * @param p pointer to aray of chars.
 */
void c_free4d(char ****pppp) {
	char ***ppp;
	char **pp;
	char *p;

	if (pppp == NULL)
		quit("free4d: attempt to free NULL pointer\n");
	/* extract pointer to main storage */
	p = pppp[0][0][0];
	/* extract pointer to row pointer storage */
	pp = pppp[0][0];
	/* extract pointer to plane pointer storage */
	ppp = pppp[0];
	/* free cube pointer space */
	free((void *) pppp);

	/* free plane pointer space */
	if (ppp == NULL)
		quit("free4d: attempt to free NULL pointer\n");
	free((void *) ppp);

	/* free row storage */
	if (pp == NULL)
		quit("free4d: attempt to free NULL pointer\n");
	free((void *) pp);

	/* free main storage */
	if (p == NULL)
		quit("free4d: attempt to free NULL pointer\n");
	free((void *) p);
}

