/**
 * \file
 * \brief Array utility functions. These functions output, initialise and allocate
 * arrays.
 *
 * \ingroup atUtil
 *
 *
 *	Revisions:
 *	Changed the fprintf statements to use 20 dp precision.
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include <signal.h>
#include <math.h>
#include <string.h>
#include <sjwlib.h>
#include <netcdf.h>
#include <atlantisboxmodel.h>
#include <atlantisMem.h>

/******************************************** Array Allocation **************************************************/

void Util_Init_5D_Double(double *****array, int dim1, int dim2, int dim3, int dim4, int dim5, double value)
{
	int i, j, k, l, m;

	for(i = 0; i < dim1; i++){
		for(j = 0; j < dim2; j++){
			for(k = 0; k < dim3; k++){
				for(l = 0; l < dim4; l++){
					for(m= 0; m < dim5; m++){
						array[i][j][k][l][m] = value;
					}
				}
			}
		}
	}
}

void Util_Init_4D_Double(double ****array, int dim1, int dim2, int dim3, int dim4, double value)
{
	int i, j, k, l;

	for(i = 0; i < dim1; i++){
		for(j = 0; j < dim2; j++){
			for(k = 0; k < dim3; k++){
				for(l = 0; l < dim4; l++){
					array[i][j][k][l] = value;
				}
			}
		}
	}
}

void Util_Init_3D_Double(double ***array, int dim1, int dim2, int dim3,  double value)
{
	int i,j, k;

	for(i = 0; i < dim1; i++){
		for(j = 0; j < dim2; j++){
			for(k = 0; k < dim3; k++){
				array[i][j][k] = value;
			}
		}
	}
}


void Util_Init_3D_Long_Double(long double ***array, int dim1, int dim2, int dim3,  long double value)
{
	int i,j, k;

	for(i = 0; i < dim1; i++){
		for(j = 0; j < dim2; j++){
			for(k = 0; k < dim3; k++){
				array[i][j][k] = value;
			}
		}
	}
}


void Util_Init_4D_Long_Double(long double ****array, int dim1, int dim2, int dim3,  int dim4, long double value)
{
	int i,j, k, l;

	for(i = 0; i < dim1; i++){
		for(j = 0; j < dim2; j++){
			for(k = 0; k < dim3; k++){
				for(l = 0; l < dim4; l++){
					array[i][j][k][l] = value;
				}
			}
		}
	}
}
void Util_Init_2D_Double(double **array, int dim1, int dim2,  double value)
{
	int i, j;

	for(i = 0; i < dim1; i++){
		for(j = 0; j < dim2; j++){
			array[i][j]= value;
		}
	}
}

void Util_Init_2D_Long_Double(long double **array, int dim1, int dim2, long double value)
{
	int i, j;

	for(i = 0; i < dim1; i++){
		for(j = 0; j < dim2; j++){
			array[i][j]= value;
		}
	}
}

void Util_Init_1D_Double(double *array, int dim1, double value)
{
	int i;

	for(i = 0; i < dim1; i++){
		array[i]= value;
	}
}

double ***Util_Alloc_Init_3D_Double(int dim1, int dim2, int dim3, double value)
{
	double ***array;

	//printf("Creating double array with dimensions %d, %d, %d\n", dim1, dim2, dim3);

	array = (double ***)alloc3d(dim1, dim2, dim3);
	Util_Init_3D_Double(array, dim3, dim2, dim1, value);

	return array;
}

long double ***Util_Alloc_Init_3D_Long_Double(int dim1, int dim2, int dim3, double value)
{
	long double ***array;

	//printf("Creating long double array with dimensions %d, %d, %d\n", dim1, dim2, dim3);

	array = (long double ***)d_alloc3longd(dim1, dim2, dim3);
	Util_Init_3D_Long_Double(array, dim3, dim2, dim1, value);

	return array;
}

double **Util_Alloc_Init_2D_Double(int dim1, int dim2,  double value)
{
	double **array;

	//printf("Creating double array with dimensions %d, %d\n", dim1, dim2);

	array = (double **)alloc2d(dim1, dim2);
	Util_Init_2D_Double(array, dim2, dim1, value);

	return array;
}

long double **Util_Alloc_Init_2D_Long_Double(int dim1, int dim2,  double value)
{
	long double **array;

	//printf("Creating long double array with dimensions %d, %d\n", dim1, dim2);

	array = (long double **)d_alloc2longd(dim1, dim2);
	Util_Init_2D_Long_Double(array, dim2, dim1, value);

	return array;
}

double ****Util_Alloc_Init_4D_Double(int dim1, int dim2, int dim3, int dim4, double value)
{
	double ****array;

	//printf("Creating double array with dimensions %d, %d, %d, %d\n", dim1, dim2, dim3, dim4);

	array = (double ****)alloc4d(dim1, dim2, dim3, dim4);
	Util_Init_4D_Double(array, dim4, dim3, dim2, dim1, value);

	return array;
}

long double ****Util_Alloc_Init_4D_Long_Double(int dim1, int dim2, int dim3, int dim4, long double value)
{
	long double ****array;

	array = (long double ****)d_alloc4longd(dim1, dim2, dim3, dim4);
	Util_Init_4D_Long_Double(array, dim4, dim3, dim2, dim1, value);

	return array;
}

double *****Util_Alloc_Init_5D_Double(int dim1, int dim2, int dim3, int dim4, int dim5, double value)
{
	double *****array;

	//printf("Creating double array with dimensions %d, %d, %d, %d, %d\n", dim1, dim2, dim3, dim4, dim5);

	array = (double *****)alloc5d(dim1, dim2, dim3, dim4, dim5);
	Util_Init_5D_Double(array, dim5, dim4, dim3, dim2, dim1, value);

	return array;
}

double *Util_Alloc_Init_1D_Double(int dim1, double value)
{
	double *array;

	//printf("Creating double array with dimensions %d\n", dim1);

	array = (double *)alloc1d(dim1);
	Util_Init_1D_Double(array, dim1, value);

	return array;
}


long double *Util_Alloc_Init_1D_Long_Double(int dim1, long double value)
{
	long double *array;
	int i;

	//printf("Creating long double array with dimensions %d\n", dim1);

	array = (long double *)d_alloc1longd(dim1);

	for(i = 0; i < dim1; i++){
		array[i]= value;
	}

	return array;
}



void Util_Init_5D_Int(int *****array, int dim1, int dim2, int dim3, int dim4, int dim5, int value)
{
	int i, j, k, l, m;

	for(i = 0; i < dim1; i++){
		for(j = 0; j < dim2; j++){
			for(k = 0; k < dim3; k++){
				for(l = 0; l < dim4; l++){
					for(m= 0; m < dim5; m++){
						array[i][j][k][l][m] = value;
					}
				}
			}
		}
	}
}

void Util_Init_4D_Int(int ****array, int dim1, int dim2, int dim3, int dim4, int value)
{
	int i, j, k, l;

	for(i = 0; i < dim1; i++){
		for(j = 0; j < dim2; j++){
			for(k = 0; k < dim3; k++){
				for(l = 0; l < dim4; l++){
					array[i][j][k][l] = value;
				}
			}
		}
	}
}

void Util_Init_3D_Int(int ***array, int dim1, int dim2, int dim3,  int value)
{
	int i, j, k;

	for(i = 0; i < dim1; i++){
		for(j = 0; j < dim2; j++){
			for(k = 0; k < dim3; k++){
				array[i][j][k]= value;
			}
		}
	}
}

void Util_Init_2D_Int(int **array, int dim1, int dim2,  int value)
{
	int i, j;

	for(i = 0; i < dim1; i++){
		for(j = 0; j < dim2; j++){
			array[i][j]= value;
		}
	}
}

void Util_Init_1D_Int(int *array, int dim1, int value)
{
	int i;

	for(i = 0; i < dim1; i++){
		array[i]= value;
	}
}

int ***Util_Alloc_Init_3D_Int(int dim1, int dim2, int dim3, int value)
{
	int ***array;

	//printf("Creating int array with dimensions %d, %d, %d\n", dim1, dim2, dim3);

	array = (int ***)alloc3d(dim1, dim2, dim3);
	if(array == NULL)
		quit("Util_Alloc_Init_3D_Int: Could not allocate array\n");
	Util_Init_3D_Int(array, dim3, dim2, dim1, value);

	return array;
}

int **Util_Alloc_Init_2D_Int(int dim1, int dim2, int value)
{
	int **array;

	//printf("Creating int array with dimensions %d, %d\n", dim1, dim2);

	array = (int **)alloc2d(dim1, dim2);
	Util_Init_2D_Int(array, dim2, dim1, value);

	return array;
}

int ****Util_Alloc_Init_4D_Int(int dim1, int dim2, int dim3, int dim4, int value)
{
	int ****array;

	//printf("Creating int array with dimensions %d, %d, %d, %d\n", dim1, dim2, dim3, dim4);

	array = (int ****)alloc4d(dim1, dim2, dim3, dim4);
	Util_Init_4D_Int(array, dim4, dim3, dim2, dim1, value);

	return array;
}
int *****Util_Alloc_Init_5D_Int(int dim1, int dim2, int dim3, int dim4, int dim5, int value)
{
	int *****array;

	//printf("Creating int array with dimensions %d, %d, %d, %d, %d\n", dim1, dim2, dim3, dim4, dim5);

	array = (int *****)alloc5d(dim1, dim2, dim3, dim4, dim5);
	Util_Init_5D_Int(array, dim5, dim4, dim3, dim2, dim1, value);

	return array;
}
int *Util_Alloc_Init_1D_Int(int dim1, int value)
{
	int *array;

	//printf("Creating int array with dimensions %d\n", dim1);

	array = (int *)alloc1d(dim1);
	Util_Init_1D_Int(array, dim1, value);

	return array;
}


/******************************************** Array Reallocation **************************************************/
double *Util_ReAlloc_1D_Double(int newdim, int olddim, double *oldarray, double value)
{
	int i;
	double *array;
	
	array = Util_Alloc_Init_1D_Double(newdim, value);

	if ( newdim > olddim ) {
		for(i = 0; i < olddim; i++){
			array[i]= oldarray[i];
		}
	} else {
		quit("Util_ReAlloc_1D_Double: Should we realloc if the dimension is smaller?\n");
	}
	free1d(oldarray);

	return array;
}

double ***Util_ReAlloc_3D_Double(int newdim1, int olddim1, int newdim2, int olddim2, int newdim3, int olddim3, double ***oldarray, double value)
{
	int i, j, k;
	double ***array;
	
	array = Util_Alloc_Init_3D_Double(newdim1, newdim2, newdim3, value);

	if (( newdim1 > olddim1 ) || ( newdim2 > olddim2 ) || ( newdim3 > olddim3 )){
		for(i = 0; i < olddim3; i++){
			for(j = 0; j < olddim2; j++){
				for(k = 0; k < olddim1; k++){
					array[i][j][k] = oldarray[i][j][k];
				}
			}
		}

	} else {
		quit("Util_ReAlloc_3D_Double: Should we realloc if the dimension is smaller?\n");
	}
	free3d(oldarray);

	return array;
}

int *Util_ReAlloc_1D_Int(int newdim, int olddim, int *oldarray, int value)
{
	int i;
	int *array;
	
	array = Util_Alloc_Init_1D_Int(newdim, value);

	if ( newdim > olddim ) {
		for(i = 0; i < olddim; i++){
			array[i]= oldarray[i];
		}
	} else {
		quit("Util_ReAlloc_1D_Int: Should we realloc if the dimension is smaller?\n");
	}
	i_free1d(oldarray);

	return array;
}

int **Util_ReAlloc_2D_Int(int newdim1, int olddim1, int newdim2, int olddim2, int **oldarray, int value)
{
	int i, j;
	int **array;
	
	array = Util_Alloc_Init_2D_Int(newdim1, newdim2, value);

	if (( newdim1 > olddim1 ) || ( newdim2 > olddim2 )){
		for(i = 0; i < olddim2; i++){
			for(j = 0; j < olddim1; j++){
				array[i][j] = oldarray[i][j];
			}
		}
	} else {
		quit("Util_ReAlloc_2D_Int: Should we realloc if the dimension is smaller?\n");
	}
	i_free2d(oldarray);

	return array;

}
