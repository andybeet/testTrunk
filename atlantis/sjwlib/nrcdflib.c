#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <float.h>
#include "sjwlib.h"


/******************************************************************************

This code comes from the 3rd edition of Numerical Recipes - including
 
NORM_CDF_VALUES returns values of the Normal CDF (from 3rd edition numerical recipes).
LOGNORM_CDF_VALUES returns values of the Log normal CDF (from 3rd edition numerical recipes).
CAUCHY_CDF_VALUES returns values of the Cauchy CDF (from 3rd edition numerical recipes).
STUDENTT_CDF_VALUES returns values of the Student t CDF (from 3rd edition numerical recipes).
LOGISTIC_CDF_VALUES returns values of the Logisitic CDF (from 3rd edition numerical recipes).
EXP_CDF_VALUES returns values of the Exponential CDF (from 3rd edition numerical recipes).
WEIBULL_CDF_VALUES returns values of the Weibull CDF (from 3rd edition numerical recipes).
CHISQ_CDF_VALUES returns values of the Chi Square CDF (from 3rd edition numerical recipes).
GAMMA_CDF_VALUES returns values of the Gamma CDF (from 3rd edition numerical recipes).
F_CDF_VALUES returns values of the F CDF (from 3rd edition numerical recipes).
BETA_CDF_VALUES returns values of the Beta CDF (from 3rd edition numerical recipes).
KS_CDF_VALUES returns values of the F CDF (from 3rd edition numerical recipes).
POISSON_NR_CDF_VALUES returns values of the F CDF (from 3rd edition numerical recipes).
BINOMIAL_NR_CDF_VALUES returns values of the F CDF (from 3rd edition numerical recipes).
 
*/



/******************************************************************************/

double norm_cdf_values ( double x, double mu, double sigma )

/******************************************************************************/
//
//  Purpose:
//
//    NORM_CDF_VALUES returns values of the Gaussian Normal CDF.
//
//  Modified:
//
//    14 June 2014
//
//  Author:
//
//    Beth Fulton
//
//  Reference:
//
//    Numerical Recipes 3rd edition.
//
//  Parameters:
//
//    Input: x = point to read cdf
//           mu = mean
//           sigma = standard deviation
//    Output: cumulative probability of a value <= x
//
{
    if (sigma <= 0)
        quit("You can not specify a negative or zero sigma when requesting normal cdf value\n");
    return (0.5 * erfc(-0.707106781186547524 * (x - mu )/ sigma));
    
}

/******************************************************************************/

double lognorm_cdf_values ( double x, double mu, double sigma )

/******************************************************************************/
//
//  Purpose:
//
//    LOGNORM_CDF_VALUES returns values of the Lognormal CDF.
//
//  Modified:
//
//    14 June 2014
//
//  Author:
//
//    Beth Fulton
//
//  Reference:
//
//    Numerical Recipes 3rd edition.
//
//  Parameters:
//
//    Input: x = point to read cdf
//           mu = mean
//           sigma = standard deviation
//    Output: cumulative probability of a value <= x
//
{
    if (sigma <= 0)
        quit("You can not specify a negative or zero sigma when requesting log normal cdf value\n");
    if (x <= 0)
        quit("You can not specify a negative zero x when requesting log normal cdf value\n");
    return (0.5 * erfc(-0.707106781186547524 * (x - mu )/ sigma));
    
}


/******************************************************************************/

double cauchy_cdf_values ( double x, double mu, double sigma )

/******************************************************************************/
//
//  Purpose:
//
//    CAUCHY_CDF_VALUES returns values of the Cauchy Distribution CDF.
//
//  Modified:
//
//    14 June 2014
//
//  Author:
//
//    Beth Fulton
//
//  Reference:
//
//    Numerical Recipes 3rd edition.
//
//  Parameters:
//
//    Input: x = point to read cdf
//           mu = mean
//           sigma = standard deviation
//    Output: cumulative probability of a value <= x
//
{
    if (sigma <= 0)
        quit("You can not specify a negative or zero sigma when requesting cauchy cdf value\n");
    
    if (x == 0.) return 0.0;
    
    return (0.5 * erfc(-0.707106781186547524 * (log(x) - mu) / sigma));
    
}

/******************************************************************************/

double studentt_cdf_values ( int tails, double t, double nu, double mu, double sigma )

/******************************************************************************/
//
//  Purpose:
//
//    STUDENTT_CDF_VALUES returns values of the Student t Distribution CDF.
//
//  Modified:
//
//    14 June 2014
//
//  Author:
//
//    Beth Fulton
//
//  Reference:
//
//    Numerical Recipes 3rd edition.
//
//  Parameters:
//
//    Input: tails = whether want a one tail or two tail cdf
//           t = point to read cdf
//           nu = the power (number of convergent moments)
//           mu = mean
//           sigma = standard deviation
//    Output: cumulative probability of a value <= x
//
// Note replaced Numericalr eceipes SQR() call with explicit squaring
{
    double p, step1;
    
    if (sigma <= 0)
        quit("You can not specify a negative or zero sigma when requesting student t cdf value\n");
    if (nu <= 0)
        quit("You can not specify a negative or zero nu when requesting student t cdf value\n");
    
    if(tails == 1) {
        step1 = (t - mu) / sigma;
        p = 0.5 * betai(0.5 * nu, 0.5, nu/(nu + (step1 * step1)));
        if (t >= mu) return 1.0 - p;
        else return p;
    } else if(tails == 2) {
        if (t < 0.0) quit("Can not have t < 0 in two tailed student t cdf");
        return (1.0 - betai(0.5 * nu, 0.5, nu / (nu + (t * t))));
    } else {
        quit("You can not ask for more than two tails for a student t distribution\n");
    }

    return 0;
}

/******************************************************************************/

double logistic_cdf_values ( double x, double mu, double sigma )

/******************************************************************************/
//
//  Purpose:
//
//    LOGISTIC_CDF_VALUES returns values of the Student t Distribution CDF.
//
//  Modified:
//
//    14 June 2014
//
//  Author:
//
//    Beth Fulton
//
//  Reference:
//
//    Numerical Recipes 3rd edition.
//
//  Parameters:
//
//    Input: x = point to read cdf
//           nu = the power (number of convergent moments)
//           mu = mean
//           sigma = standard deviation
//    Output: cumulative probability of a value <= x
//
{
    double p;
    
    if (sigma <= 0)
        quit("You can not specify a negative or zero sigma when requesting logistic cdf value\n");

    p = exp(-fabs(1.81379936423421785 * (x - mu) / sigma));
    if (x >= mu) return (1.0 / (1.0 + p));    //Because we used abs to control overelse
    return (p / (1.0 + p));
    
    
}

/******************************************************************************/
                         
double exp_cdf_values ( double x, double beta )
                         
/******************************************************************************/
//
//  Purpose:
//
//    EXP_CDF_VALUES returns some values of the Exponential CDF.
//
//  Modified:
//
//    14 June 2014
//
//  Author:
//
//    Beth Fulton
//
//  Reference:
//
//    Numerical Recipes 3rd edition.
//
//  Parameters:
//
//    Input: x = point to read cdf
//           beta = width of the distribution
//    Output: cumulative probability of a value <= x
//
{
    if (x < 0.0) quit("Can not have negative exponential distribution");
    return (1.0 - exp(-beta * x));
}
            
/******************************************************************************/

double weibull_cdf_values ( double x, double beta, double alpha )

/******************************************************************************/
//
//  Purpose:
//
//  WEIBULL_CDF_VALUES returns some values of the Weibull CDF.
//
//  Modified:
//
//    14 June 2014
//
//  Author:
//
//    Beth Fulton
//
//  Reference:
//
//    Numerical Recipes 3rd edition.
//
//  Parameters:
//
//    Input: x = point to read cdf
//           beta = width of the distribution
//           alpha = shape of failure likelihood ( alpha > 0; with alpha > 1 more likely with age, alpha~(0,1) means less likely with age)
//    Output: cumulative probability of a value <= x
//
{
    if (x < 0.0) quit("Can not have negative exponential distribution");
    return (1.0 - exp(pow((- x / beta), alpha)));
}
            
/******************************************************************************/

double chisq_cdf_values ( double x, double nu )

/******************************************************************************/
//
//  Purpose:
//
//  WEIBULL_CDF_VALUES returns some values of the Chi Square CDF.
//
//  Modified:
//
//    14 June 2014
//
//  Author:
//
//    Beth Fulton
//
//  Reference:
//
//    Numerical Recipes 3rd edition.
//
//  Parameters:
//
//    Input: x = point to read cdf
//           nu = degrees of freedom
//    Output: cumulative probability of a value <= x
//
{
    if (x < 0.0) quit("Can not have negative chi square distribution");
    return gammp(0.5 * nu, 0.5 * x);
}

/******************************************************************************/

double gamma_cdf_values ( double x, double alpha, double beta )

/******************************************************************************/
//
//  Purpose:
//
//  GAMMA_CDF_VALUES returns some values of the Gamma CDF.
//
//  Modified:
//
//    14 June 2014
//
//  Author:
//
//    Beth Fulton
//
//  Reference:
//
//    Numerical Recipes 3rd edition.
//
//  Parameters:
//
//    Input: x = point to read cdf
//           alpha, beta = shape defining parameters
//    Output: cumulative probability of a value <= x
//
{
    if (alpha <= 0)
        quit("You can not specify a negative or zero alpha when requesting gamma cdf value\n");
    if (beta <= 0)
        quit("You can not specify a negative or zero beta when requesting gamma cdf value\n");

    if (x < 0.0) quit("Can not have negative gamma distribution");
    return gammp(alpha, beta * x);
}

/******************************************************************************/

double F_cdf_values ( double x, double nu1, double nu2 )

/******************************************************************************/
//
//  Purpose:
//
//  F_CDF_VALUES returns some values of the F CDF.
//
//  Modified:
//
//    14 June 2014
//
//  Author:
//
//    Beth Fulton
//
//  Reference:
//
//    Numerical Recipes 3rd edition.
//
//  Parameters:
//
//    Input: x = point to read cdf
//           nu1, nu2 = shape defining parameters
//    Output: cumulative probability of a value <= x
//
{
    if (nu1 <= 0)
        quit("You can not specify a negative or zero nu1 when requesting F distribution cdf value\n");
    if (nu2 <= 0)
        quit("You can not specify a negative or zero nu2 when requesting F distribution cdf value\n");

    if (x < 0.0) quit("Can not have negative F distribution");
    return (betai(0.5 * nu1, 0.5 * nu2, nu1 * x/ (nu2 + nu1 * x)));
}

/******************************************************************************/

double beta_cdf_values ( double x, double alpha, double beta )

/******************************************************************************/
//
//  Purpose:
//
//  F_CDF_VALUES returns some values of the Gamma CDF.
//
//  Modified:
//
//    14 June 2014
//
//  Author:
//
//    Beth Fulton
//
//  Reference:
//
//    Numerical Recipes 3rd edition.
//
//  Parameters:
//
//    Input: x = point to read cdf
//           alpha, beta = shape defining parameters
//    Output: cumulative probability of a value <= x
//
{
    if (alpha <= 0)
        quit("You can not specify a negative or zero nu1 when requesting F distribution cdf value\n");
    if (beta <= 0)
        quit("You can not specify a negative or zero nu2 when requesting F distribution cdf value\n");
    
    if (x < 0.0) quit("Can not have negative beta distribution");
    if (x > 1.0) quit("Can not x > 1 in beta distribution");
    
    return (betai(alpha, beta, x));
}

/******************************************************************************/

double ks_cdf_values ( double z )

/******************************************************************************/
//
//  Purpose:
//
//  KS_CDF_VALUES returns some values of the Kolmogorov-Smirnov CDF.
//
//  Modified:
//
//    14 June 2014
//
//  Author:
//
//    Beth Fulton
//
//  Reference:
//
//    Numerical Recipes 3rd edition.
//
//  Parameters:
//
//    Input: z = point to read cdf
//    Output: cumulative probability of a value <= z
//
// Note replaced Numericalr eceipes SQR() call with explicit squaring
{
    double x, y;
    if (z < 0.0) quit("Can not have negative Kolmogorov-Smirnov distribution");
    
    if (z == 0.0) return 0.;
    
    if (z < 1.18) {
        y = exp(-1.23370055013616983 / (z * z));
        return (2.25675833419102515 * sqrt(-log(y)) *(y + pow(y,9) + pow(y,25) + pow(y,49)));
    } else {
        x = exp(-2.0 * (z * z));
        return (1.0 - 2.0 * (x - pow(x,4) + pow(x,9)));
    }
}

/******************************************************************************/

double poisson_NR_cdf_values ( int n, double lambda )

/******************************************************************************/
//
//  Purpose:
//
//  POISSON_NR_CDF_VALUES returns some values of the Poisson CDF.
//
//  Modified:
//
//    14 June 2014
//
//  Author:
//
//    Beth Fulton
//
//  Reference:
//
//    Numerical Recipes 3rd edition.
//
//  Parameters:
//
//    Input: n = number of events
//           lambda = mean number of expected events
//    Output: cumulative probability n events
//
{
    double ndub = (double)(n);
    
    if (lambda <= 0)
        quit("You can not specify a negative or zero lambda when requesting poisson distribution cdf value\n");

    if (n < 0.0) quit("Can not have negative Poisson distribution");
    if (n == 0) return 0.;
    
    return gammq(ndub, lambda);
}

/******************************************************************************/

double binomial_NR_cdf_values ( int k, int n, double p )

/******************************************************************************/
//
//  Purpose:
//
//  BIONOMIAL_NR_CDF_VALUES returns some values of the Binomial CDF.
//
//  Modified:
//
//    14 June 2014
//
//  Author:
//
//    Beth Fulton
//
//  Reference:
//
//    Numerical Recipes 3rd edition.
//
//  Parameters:
//
//    Input: k = number of events
//           n = the number of event_x in the sequence of k events
//           p = probability of event_x in any one instance
//    Output: cumulative probability n events
//
{
    double dubk = (double)(k);
    if (n <= 0)
        quit("You can not specify a negative or zero sample size n when requesting binomial distribution cdf value\n");
    if (p <= 0.0)
            quit("You can not specify a negative or zero probability of an event when requesting binomial distribution cdf value\n");
    if (p >= 1.0)
            quit("You can not specify a >= 1.0 probability of an event when requesting binomial distribution cdf value\n");
        
    if (k < 0) quit("Can not have negative Binomial distribution");
    if (k == 0) return 0.0;
    if (k > n) return 1.0;
    return (1.0 - betai(dubk, n - k + 1.0, p));

}

/******************************************************************************/

double negbinomial_cdf_values ( int k, int r, double p )

/******************************************************************************/
//
//  Purpose:
//
//  NEGATIVE_BIONOMIAL_CDF_VALUES returns some values of the Negative Binomial CDF.
//
//  Modified:
//
//    14 June 2014
//
//  Author:
//
//    Beth Fulton
//
//  Reference:
//
//    Numerical Recipes 3rd edition (for inspiration and wikipedia for the cdf formula).
//
//  Parameters:
//
//    Input: k = number of events
//           r = the number of trials
//           p = probability of event in any one instance
//    Output: cumulative probability n events
//
{
    double dubk = (double)(k + 1);
    if (r <= 0)
        quit("You can not specify a negative or zero sample size r when requesting negative binomial distribution cdf value\n");
    if (p <= 0.0)
        quit("You can not specify a negative or zero probability of an event when requesting negative binomial distribution cdf value\n");
    if (p >= 1.0)
        quit("You can not specify a >= 1.0 probability of an event when requesting negative binomial distribution cdf value\n");
                
    if (k < 0) quit("Can not have negative binomial distribution");
    return (1.0 - betai(dubk, r, p));
    
}

/******************************************************************************
 *
 * \brief Returns incomplete beta function Ix(a,b) for positive a and b, and x between 0 and 1.
 *
 *  x is the probabilty of an event
 *
 *  Sourced from Numerical Recipes 3rd ed.
 */
double betai(double a, double b, double x) {
    double bt;
    if ((a <= 0.0) || (b <= 0.0)) quit("You can not have negative parameters for betai\n");
    if ((x < 0.0) || (x > 1.0)) quit("Yu can npt haev a probabilty <0 or > 1 x in routine betai");
    if ((x == 0.0) || (x == 1.0)) return x;
    if (a > 3000 && b > 300) return betaiapprox(a,b,x);   // Switch to quadrature method.
        bt = exp(gammln(a+b)- gammln(a) - gammln(b) + a * log(x) + b * log(1.0 - x));
    if (x < (a + 1.0) / (a + b + 2.0))
        return (bt * betacf(a,b,x) / a);
    else
        return (1.0 - bt * betacf(b,a,1.0 - x) / b);
}

/******************************************************************************
 *
 * \brief Evaluates continued fraction for incomplete beta function by modified Lentz’s method
 *  Should only be called from within betai
 *
 *  Sourced from Numerical Recipes 3rd ed.
 */
double betacf(double a, double b, double x) {
    int m, m2;
    double aa, c, d, del, h, qab, qam, qap;
    double EPS = DBL_EPSILON;                           // numeric_limits<double>::epsilon();
    double FPMIN = DBL_MIN / EPS;                       // numeric_limits<double>::min()/EPS;
    
    qab = a + b;        //These q’s will be used in factors that occur in the coefficients
    qap = a + 1.0;
    qam = a - 1.0;
    c = 1.0;              // First step of Lentz’s method.
    d = 1.0 - qab * x / qap;
    if (fabs(d) < FPMIN)
        d = FPMIN;
    d = 1.0 / d;
    h = d;
    for (m=1; m < 10000; m++) {
        m2 = 2 * m;
        aa = m * (b - m) * x / ((qam + m2) * (a + m2));
        d = 1.0 + aa * d;   // One step (the even one) of the recurrence.
        if (fabs(d) < FPMIN)
            d = FPMIN;
        c = 1.0 + aa / c;
        if (fabs(c) < FPMIN)
            c = FPMIN;
        d = 1.0 / d;
        h *= (d * c);
        aa = -(a + m) * (qab + m) * x/ ((a + m2) * (qap + m2));
        d = 1.0 + aa * d;   // Next step of the recurrence (the odd one).
        if (fabs(d) < FPMIN)
            d = FPMIN;
        c = 1.0 + aa / c;
        if (fabs(c) < FPMIN)
            c = FPMIN;
        d = 1.0 / d;
        del = d * c;
        h *= del;
        if (fabs(del - 1.0) <= EPS)
            break;   //Are we done?
    }
    return h;
}
/******************************************************************************
 *
 * \brief Evaluates continued fraction for incomplete beta function by quadrature. Returns Ix(a,b)
 *
 *  Should only be called from within betai
 *
 *  Sourced from Numerical Recipes 3rd ed.
 *
 * Note replaced Numericalr eceipes SQR() call with explicit squaring
 */
double betaiapprox(double a, double b, double x) {
    int j;
    double xu, t, sum, ans;
    double a1 = a - 1.0;
    double b1 = b - 1.0;
    double mu = a / ( a + b);
    double lnmu = log(mu);
    double lnmuc = log(1.0 - mu);
    
    //Abscissas and weights for Gauss-Legendre quadrature.
    double y[18] = {0.0021695375159141994,
        0.011413521097787704,0.027972308950302116,0.051727015600492421,
        0.082502225484340941, 0.12007019910960293,0.16415283300752470,
        0.21442376986779355, 0.27051082840644336, 0.33199876341447887,
        0.39843234186401943, 0.46931971407375483, 0.54413605556657973,
        0.62232745288031077, 0.70331500465597174, 0.78649910768313447,
        0.87126389619061517, 0.95698180152629142};
    
   double w[18] = {0.0055657196642445571,
        0.012915947284065419,0.020181515297735382,0.027298621498568734,
        0.034213810770299537,0.040875750923643261,0.047235083490265582,
        0.053244713977759692,0.058860144245324798,0.064039797355015485,
        0.068745323835736408,0.072941885005653087,0.076598410645870640,
        0.079687828912071670,0.082187266704339706,0.084078218979661945,
        0.085346685739338721,0.085983275670394821};
    
    
    t = sqrt(a * b / (( a + b ) * ( a + b ) * ( a + b + 1.0)));
    if (x > (a / (a + b))) {                  //Set how far to integrate into the tail:
        if (x >= 1.0) return 1.0;
        xu = min(1.0, max(mu + 10.0 * t, x + 5.0 * t));
    } else {
       if (x <= 0.0) return 0.0;
        xu = max(0.0, min(mu - 10.0 * t, x - 5.0 * t));
    }
    sum = 0;
    for (j=0; j < 18; j++) {                    // Gauss-Legendre.
        t = x + (xu - x) * y[j];
        sum += w[j] * exp(a1 * (log(t) - lnmu) + b1 * (log(1-t) - lnmuc));
    }
    ans = sum * (xu - x) * exp(a1 * lnmu - gammln(a) + b1 * lnmuc - gammln(b) + gammln(a+b));
    return (ans>0.0? 1.0-ans : -ans);
}

/******************************************************************************
 *
 * \brief Returns the incomplete gamma function P(a, x).
 *
 *  Sourced from Numerical Recipes 3rd ed.
 */
double gammp(double a, double x) {
    int inta = (int)(a);
    
    if ((x < 0.0) || (a <= 0.0)) quit("Bad args in gammp\n");
    if (x == 0.0) return 0.0;
    else if (inta >= 100) return gammpapprox(a,x,1);    // Quadrature.
    else if (x < (a + 1.0)) return gser(a,x);           // Use the series representation.
    else return (1.0 - gcf(a,x));                       // Use the continued fraction representation.
}

/******************************************************************************
 *
 * \brief Returns the incomplete gamma function Q(a,x) ~ 1 - P(a,x)
 *
 *  Sourced from Numerical Recipes 3rd ed.
 */
double gammq(double a, double x) {
    int inta = (int)(a);
    
    if ((x < 0.0) || (a <= 0.0)) quit("Bad args in gammq");
    if (x == 0.0) return 1.0;
    else if (inta >= 100) return gammpapprox(a,x,0);        // Quadrature.
    else if (x < (a + 1.0)) return 1.0 - gser(a,x);         // Use the series representation.
    else return gcf(a,x);                                   // Use the continued fraction representation.
}


/******************************************************************************
 *
 * \brief Returns the incomplete gamma function P(a,x) evaluated by its series representation.
 * Also sets lnG(a) as gln. 
 *
 * Called from within other gamma functions
 *
 *  Sourced from Numerical Recipes 3rd ed.
 */
double gser(double a, double x) {
    double sum, del, ap;
    double EPS = DBL_EPSILON;                           // numeric_limits<double>::epsilon();
    //double FPMIN = DBL_MIN / EPS;                       // numeric_limits<double>::min()/EPS;
    double gln = gammln(a);
    
    ap = a;
    del = sum = 1.0 / a;
    for (;;) {
        ++ap;
        del *= x/ap;
        sum += del;
        if (fabs(del) < fabs(sum) * EPS) {
            return sum * exp(-x + a * log(x) - gln);
        }
    }
}
/******************************************************************************
 *
 * \brief Returns the incomplete gamma function Q(a, x) evaluated by its continued fraction representation
 * Also sets lnG(a) as gln.
 *
 * Called from within other gamma functions
 *
 *  Sourced from Numerical Recipes 3rd ed.
 */
double gcf(double a, double x) {
    
    double EPS = DBL_EPSILON;                           // numeric_limits<double>::epsilon();
    double FPMIN = DBL_MIN / EPS;                       // numeric_limits<double>::min()/EPS;
    int i;
    double an, b, c, d, del, h;
    double gln = gammln(a);
    
    b = x + 1.0 - a;            // Set up for evaluating continued fraction by modified Lentz’s method with b0 = 0.
    c = 1.0 / FPMIN;
    d = 1.0 / b;
    h = d;
    for (i = 1; ; i++) {        // Iterate to convergence.
        an = -i*(i-a);
        b += 2.0;
        d = an * d + b;
        if (fabs(d) < FPMIN)
            d = FPMIN;
        c = b + an / c;
        if (fabs(c) < FPMIN)
            c = FPMIN;
        d = 1.0 / d;
        del = d * c;
        h *= del;
        if (fabs(del-1.0) <= EPS)
            break;
    }
    return (exp(-x + a * log(x) - gln) * h);    // Put factors in front.
}

/******************************************************************************
 *
 * \brief Incomplete gamma by quadrature. Returns P(a,x) or Q(a,x) when psig is 1 or 0, respectively.
 *
 * Called from within other gamma functions
 *
 *  Sourced from Numerical Recipes 3rd ed.
 */

double gammpapprox(double a, double x, int psig) {
    int j;
    double xu, t, sum, ans;
    double a1 = a - 1.0;
    double lna1 = log(a1);
    double sqrta1 = sqrt(a1);
    double gln = gammln(a);
    
    //Abscissas and weights for Gauss-Legendre quadrature.
    int ngau = 18;
    double y[18] = {0.0021695375159141994, 0.011413521097787704, 0.027972308950302116, 0.051727015600492421, 0.082502225484340941, 0.12007019910960293, 0.16415283300752470, 0.21442376986779355, 0.27051082840644336, 0.33199876341447887, 0.39843234186401943, 0.46931971407375483, 0.54413605556657973, 0.62232745288031077, 0.70331500465597174, 0.78649910768313447, 0.87126389619061517, 0.95698180152629142};
    
    double w[18] = {0.0055657196642445571, 0.012915947284065419,0.020181515297735382,0.027298621498568734,0.034213810770299537,0.040875750923643261,0.047235083490265582,0.053244713977759692,0.058860144245324798,0.064039797355015485,0.068745323835736408,0.072941885005653087,0.076598410645870640,0.079687828912071670,0.082187266704339706,0.084078218979661945, 0.085346685739338721,0.085983275670394821};

    // Set how far to integrate into the tail:
    if (x > a1) xu = max(a1 + 11.5 * sqrta1, x + 6.0 * sqrta1);
    else xu = max(0.0, min(a1 - 7.5 * sqrta1, x - 5.0 * sqrta1));
    sum = 0;
    for (j = 0; j < ngau; j++) {        // Gauss-Legendre.
        t = x + (xu - x) * y[j];
        sum += w[j] * exp(-(t - a1) + a1 * (log(t) - lna1));
    }
    ans = sum * (xu - x) * exp(a1 * (lna1 - 1.0) - gln);
    return (psig?(ans>0.0? 1.0-ans:-ans):(ans>=0.0? ans:1.0+ans));
}

/******************************************************************************
 *
 * \brief Returns the value ln(G(xx) for xx > 0.
 *
 * Called from within other gamma functions
 *
 *  Sourced from Numerical Recipes 3rd ed.
 */

double gammln(double xx) {
    int j;
    double x, tmp, y, ser;
    double cof[14] = {57.1562356658629235,-59.5979603554754912,
        14.1360979747417471,-0.491913816097620199,.339946499848118887e-4,
        .465236289270485756e-4,-.983744753048795646e-4,.158088703224912494e-3,
        -.210264441724104883e-3,.217439618115212643e-3,-.164318106536763890e-3,
        .844182239838527433e-4,-.261908384015814087e-4,.368991826595316234e-5};
    
    if (xx <= 0) quit("Bad arg in gammln\n");
    y = x = xx;
    tmp = x + 5.24218750000000000;              // Rational 671/128.
    tmp = (x + 0.5) * log(tmp) - tmp;
    ser = 0.999999999999997092;
    for (j=0; j < 14; j++)
        ser += cof[j] / ++y;
    return (tmp + log(2.5066282746310005 * ser / x));
}
