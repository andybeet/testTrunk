/*********************************************************************

    File:           sjwalloc.h (memory.h)
    
    Created:        Thu Aug 13 21:24:33 EST 1992
    
    Author:         Stephen Walker
                    CSIRO Division of Oceanography
    
    Purpose:        Include file for memory allocation library
    
    Revisions:      None

    $Id: sjwalloc.h 3302 2012-08-02 02:10:43Z gor171 $


	Changes:
		Bec Gorton
		Moved the 5 and 6 D array allocation functions into this library instead of the
		assessment lib.
*********************************************************************/

#ifdef __cplusplus
extern "C" {
#endif

double *alloc1Input(long int n1);
void free1dInput(double *p);

double ** alloc2dInput(long int n1, long int n2);
void free2dInput(double **pp);


double	*d_alloc1d(long n1);
void	d_free1d(double *p);
double	**d_alloc2d(long n1, long n2);
void	d_free2d(double **pp);
double	***d_alloc3d(long n1, long n2, long n3);
void	d_free3d(double ***ppp);
double	****d_alloc4d(long n1, long n2, long n3, long n4);
void	d_free4d(double ****ppp);


long double  *d_alloc1longd(long n1);
void    d_free1longd(long double *p);
long double  **d_alloc2longd(long n1, long n2);
void    d_free2longd(long double **pp);
long double  ***d_alloc3longd(long n1, long n2, long n3);
void    d_free3longd(long double ***ppp);
long double  ****d_alloc4longd(long n1, long n2, long n3, long n4);
void    d_free4longd(long double ****ppp);
long double *****d_alloc5longd(long n1, long n2, long n3, long n4, long n5);
long double ******d_alloc6longd(long n1, long n2, long n3, long n4, long n5, long n6);
void d_free5longd(long double *****ppppp);
void d_free6longd(long double ******pppppp);


/* Memory handling routines */
// Adding in a defined value so we can work out whick version of Atlantis people are using. Fixes issues in outputSplitter if people are using these libs
// instead of the old ones.
#define _ATLANTIS2_LIBS
double *****d_alloc5d(long n1, long n2, long n3, long n4, long n5);
double ******d_alloc6d(long n1, long n2, long n3, long n4, long n5, long n6);
void d_free5d(double *****ppppp);
void d_free6d(double ******pppppp);



float	*f_alloc1d(long n1);
void	f_free1d(float *p);
float	**f_alloc2d(long n1, long n2);
void	f_free2d(float **pp);
float	***f_alloc3d(long n1, long n2, long n3);
void	f_free3d(float ***ppp);
float	****f_alloc4d(long n1, long n2, long n3, long n4);
void	f_free4d(float ****ppp);
long	*l_alloc1d(long n1);
void	l_free1d(long *p);
long	**l_alloc2d(long n1, long n2);
void	l_free2d(long **pp);
long	***l_alloc3d(long n1, long n2, long n3);
void	l_free3d(long ***ppp);
int	*i_alloc1d(long n1);
void	i_free1d(int *p);
int	**i_alloc2d(long n1, long n2);
void	i_free2d(int **pp);
int	***i_alloc3d(long n1, long n2, long n3);
void	i_free3d(int ***ppp);
int	****i_alloc4d(long n1, long n2, long n3, long n4);
void	i_free4d(int ****ppp);
short	*s_alloc1d(long n1);
void	s_free1d(short *p);
short	**s_alloc2d(long n1, long n2);
void	s_free2d(short **pp);
short	***s_alloc3d(long n1, long n2, long n3);
void	s_free3d(short ***ppp);
short	****s_alloc4d(long n1, long n2, long n3, long n4);
void	s_free4d(short ****ppp);
char	*c_alloc1d(long n1);
void	c_free1d(char *p);
char	**c_alloc2d(long n1, long n2);
void	c_free2d(char **pp);
char	***c_alloc3d(long n1, long n2, long n3);
void	c_free3d(char ***ppp);
char	****c_alloc4d(long n1, long n2, long n3, long n4);
void	c_free4d(char ****ppp);

//#define _ATLANTIS_LONG_DOUBLE

//#ifdef  _ATLANTIS_LONG_DOUBLE
//#define double  long double
//#define double  double
//
//#define alloc1d alloc1longd
//#define free1d free1longd
//
//#define alloc2d alloc2longd
//#define free2d free2longd
//
//#define alloc3d alloc3longd
//#define free3d free3longd
//
//#define alloc4d alloc4longd
//#define free4d free4longd
//
//#else
//
//#define alloc1d alloc1d
//#define free1d free1d
//
//#define alloc2d alloc2d
//#define free2d free2d
//
//#define alloc3d alloc3d
//#define free3d free3d
//
//#define alloc4d alloc4d
//#define free4d free4d
//
//#endif  /* defined(double) */


#define alloc1dInput d_alloc1d
#define free1dInput d_free1d

#define alloc2dInput d_alloc2d
#define free2dInput d_free2d

#define alloc3dInput d_alloc3d
#define free3dInput d_free3d


#define alloc4dInput d_alloc4d
#define free4dInput d_free4d

#ifdef  __cplusplus
}
#endif
