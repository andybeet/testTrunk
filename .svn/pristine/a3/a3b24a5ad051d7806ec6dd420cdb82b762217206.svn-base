/**
 * This header file defines the precision of the model.
 * Be warned that if you want to change the precision from long double to double you will need to do a global find and replace on
 * all printf statements that print out the double values.
 *
 * Its is highly recommended that you turn on the WARN +=-Werror compiler flag. This will treat all warnings as errors.
 *
 *
 * You will need to find and replace the following:
 *
 * %e - change to %e
 * %f - change to %f
 *%.10f%.10f
 *%.20e%.20e
 */

#ifndef ATLANTISMEM_H_
#define ATLANTISMEM_H_

//#define _ATLANTIS_LONG_DOUBLE

#ifdef  _ATLANTIS_LONG_DOUBLE
//#define double  long double
//#define double  double

#define alloc1d d_alloc1longd
#define free1d d_free1longd

#define alloc2d d_alloc2longd
#define free2d d_free2longd

#define alloc3d d_alloc3longd
#define free3d d_free3longd

#define alloc4d d_alloc4longd
#define free4d d_free4longd


#define alloc5d d_alloc5longd
#define free5d d_free5longd

#define alloc6d d_alloc6longd
#define free6d d_free6longd


#else

#define alloc1d d_alloc1d
#define free1d d_free1d

#define alloc2d d_alloc2d
#define free2d d_free2d

#define alloc3d d_alloc3d
#define free3d d_free3d

#define alloc4d d_alloc4d
#define free4d d_free4d

#define alloc5d d_alloc5d
#define free5d d_free5d

#define alloc6d d_alloc6d
#define free6d d_free6d

#endif  /* defined(double) */


#define d_alloc1dInput d_alloc1d
#define d_free1dInput d_free1d

#define d_alloc2dInput d_alloc2d
#define d_free2dInput d_free2d

#define d_alloc3dInput d_alloc3d
#define d_free3dInput d_free3d


#define d_alloc4dInput d_alloc4d
#define d_free4dInput d_free4d



#endif /*ATLANTISMEM_H_*/
