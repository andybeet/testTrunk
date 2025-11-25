/*

    File:           diffusion.c
    
    Created:        Tue Dec 6 11:55:26 EST 1994
    
    Author:         Stephen Walker
                    CSIRO Division of Oceanography

    Returns:        void
    
    Revisions:      none

    $Id: diffusion.c 2761 2011-07-01 04:35:25Z gor171 $

*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sjwlib.h>
#include <atlantisboxmodel.h>
#include <atUtilLib.h>

/* Prototype for tri-diagonal system solver */
void
solvetri(double *Cm1, double *C, double *Cp1, double *rhs,
	 double *x,int imin,int imax);

/* Local temporary storage */
static int max_n = -1;
static double *Cm1;
static double *C;
static double *Cp1;
static double *rhs;
static double *ud;


/** Calculates 1 time step for the diffusion equation.
  *
\begin{verbatim}
dC     d (  dC)
--  = ---(K --)
dt     dx(  dx)
\end{verbatim}
  *		    
  * This routine uses a weighted method using
  * a times the values at the new time step plus
  * (1-a) times the values at the current time
  * step to calculate the spatial derivatives.
  * This leads to a tri-diagonal system of
  * equations to solve for C.
  *
  * Note that consistent units must be used for
  * all these quantities. If SI units are used,
  * then the units are as follows:
  *
  * \begin{itemize}
  * x    - metres (m)
  * dt   - seconds (s)
  * k[0] - units per square metre per second (m-2s-1)
  * k[n] - units per square metre per second (m-2s-1)
  * k[i] - metres squared per second (m2s-1)
  * \end{itemize}
  *    
  * @param n number of concentration values
  * @param c array of concentration values. (n values)
  * @param xc array of x coords at c points (n values)
  * @param k array of diffusion coefficients. This array has n+1 values.
  *	 Values k[0] and k[n] are the boundary fluxes in +ve
  *	 x direction. Other values k[i] are diffusion coefficients at
  *	 xk[i], which is between points xc[i] and xc[i-1]
  * @param xk array of x coords at k points (n+1 values)
  * @param dt time step
  * @param a weighting value, in the range [0,1].
  *	  (a = 1.0) == fully implicit,
  *	  (a = 0.5) == semi implicit
  *	  (a = 0.0) == explicit
  */
                    
void
diffusion1d(int n, double *c, double *xc, double *k,
	    double *xk, double dt, double a)
{
    int i;
    double dx;
    double v;

    /* Sanity checks */
    if( n < 1 )
	quit("diffusion1d: n < 1 (no points!)\n");
    if( a < 0.0 || a > 1.0 )
	quit("diffusion1d: weight value must be in range [0,1]\n");

    /* Only 1 layer - simple calculation */
    if( n == 1 ) {
	c[0] += dt*(k[0]-k[1])/(xk[1]-xk[0] + small_num);
	return;
    }

    /* Multiple layers - allocate temporary storage if necessary */
    if( n > max_n ) {
	if( Cm1 != NULL ) {
	    free1d(Cm1);
	    free1d(C);
	    free1d(Cp1);
	    free1d(rhs);
	    free1d(ud);
	}
	Cm1 = alloc1d(n+1);
	C = alloc1d(n+1);
	Cp1 = alloc1d(n+1);
	rhs = alloc1d(n+1);
	ud = alloc1d(n+1);
	max_n = n;
    }

    /* i=0 boundary */
    i = 0;
    dx = xk[i+1]-xk[i] + small_num;
    v = dt*k[i+1]/((xc[i+1]-xc[i] + small_num)*dx);
    Cm1[i] = 0.0;
    C[i] = 1 + a*v;
    Cp1[i] = -a*v;
    rhs[i] = (1-(1-a)*v)*c[i] + (1-a)*v*c[i+1] + dt*k[i]/dx;

    /* middle points */
    for(i=1; i<n-1; i++) {
	double dxi = xk[i+1]-xk[i] + small_num;
	double vm = dt*k[i]/((xc[i]-xc[i-1] + small_num)*dxi);
	double vp = dt*k[i+1]/((xc[i+1]-xc[i] + small_num)*dxi);

	Cm1[i] = -a*vm;
	C[i] = 1.0 + a*(vm+vp);
	Cp1[i] = -a*vp;
	rhs[i] =  (1-a)*vm*c[i-1]
		+ (1 - (1-a)*(vm+vp))*c[i]
		+ (1-a)*vp*c[i+1];
    }

    /* i=n-1 boundary */
    i = n-1;
    dx = xk[i+1]-xk[i] + small_num;
    v = dt*k[i]/((xc[i]-xc[i-1] + small_num)*dx);
    Cm1[i] = -a*v;
    C[i] = 1 + a*v;
    Cp1[i] = 0.0;
    rhs[i] = (1-a)*v*c[i-1] + (1-(1-a)*v)*c[i] - dt*k[i+1]/dx;

    /* Solve the system */
    solvetri(Cm1,C,Cp1,rhs,c,0,n-1);
}


/*
  Routine to solve tridiagonal system of equations
  Arguments:

  Cm1[i] - coefficient of Xi-1 in ith eqn (lower diagonal)
  C[i]   - coefficients of Xi in ith eqn  (diagonal)
  Cp1[i] - coefficients of Xi+1 in ith eqn (upper diagonal)
  rhs[i] - right hand side of ith eqn
  x[i]   - where to store solved Xi
  imin   - minimum index
  imax   - maximum index

*/
void
solvetri(double *Cm1, double *C, double *Cp1, double *rhs,
	 double *x, int imin, int imax)
{
    int i;
    double div;

    /* check indices - note imin == imax => 1 trivial eqn */
    if( imin >= imax || imin < 0 )
	 quit("solvetri: bad index values\n");
    
    /* DEBUGGING
    for(i=imin; i<=imax; i++)
	fprintf(stderr,"%e %e %e %e\n",Cm1[i],C[i],Cp1[i],rhs[i]);
    */

    div = C[imin];
    if( div == 0.0 )
	quit("solvetri: zero first coefficient\n");
    x[imin] = rhs[imin]/div;
    for(i=imin+1; i<=imax; i++) {
	ud[i] = Cp1[i-1]/div;
	div = C[i]-Cm1[i]*ud[i];
	if( div == 0.0 )
	    quit("solvetri: zero divisor\n");
	x[i] = (rhs[i]-Cm1[i]*x[i-1])/div;
    }
    for(i=imax-1; i>=imin; i--)
	x[i] -= ud[i+1]*x[i+1];
}

/**
 * Free up the diffusion arrays that have been allocated.
 *
 *
 */
void free_diffusion1d() {

	if (Cm1 != NULL) {
		free1d(Cm1);
		free1d(C);
		free1d(Cp1);
		free1d(rhs);
		free1d(ud);
	}
}
