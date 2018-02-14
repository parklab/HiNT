/*A library for calculating gamma function and incomplete gamma function
*/

#ifndef GAMMA_HX_XI
#define GAMMA_HX_XI

#include <math.h>
#include <float.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <sys/time.h>


double gammln(const double x); /*log gamma function*/

double gammp(const double a, const double x); /*lower incomplete function*/
/* 1/gamma(x) int_0^x exp(-t)t^{a-1}dt */
double gammq(const double a, const double x);
/*upper incomplete function*/

double pchisq(int k, const double x, int lower_tail); 
/* calculate Pr(X<=x) or Pr(X>x) for chisquare distribution with k degree of freedom
   if lower_tail==1, calculate Pr(X<=x) otherwise Pr(X>x)
*/

double qchisq_bisection(int k, const double p_in, int lower_tail);
double qchisq_newton(int k, const double p, int lower_tail);
/*use Newton-Raphson to determin x such that Pr(X<=x)=p or Pr(X>x)=p for chi-sq distribution with k degrees of freedom*/
double qchisq(int k, const double p, int lower_tail);

double betai(const double a, const double b, const double x);
/*Calculate the incomplete beta function I_x(a,b) = 1/B(a,b) \int_0^xt^{a-1}(1-t)^{b-1}dt, a,b>0*/
/*Based on this function, one can easily calculate cumulative distribution function of a binormial distribution, i.e
  for binomial distribution Binom(n,p), we have
  Pr(X>=k) = I_p(k,n-k+1)
*/

double pbinom(const double k, const double n, const double p, const int lower_tail);
/* calculate Pr(X<=k) or Pr(X>=k) for Binom(n,k) depending on lower_tial ==1 or lower_tail==0; 
   In principle, k and n should be integers, but this function can also calculate the similar value for non-integers
*/




/*The following functions is for generating random numbers*/


void seed_set(int seed); /*reset the seed of the random number generator rand_lp as seed.*/
double rand_lp();
/*Long period (> 2 * 10^18) random number generator of L'Ecuyer with Bays-Durham shuffle and added safeguards. Returns a uniform random deviate between 0.0 and 1.0 (exclusive of the end point values). Call with idum a negative integer to initialize; thereafter, do not alter idum between successive deviates in a sequence. RNMX should approximate the largest floating value that is less than 1
See Numerical Recipes in C++, page 286.
*/

int genSeed(); /*generate a seed for the random number generator*/

void db_shuffle(double *array, int n);
void db_shuffle_int(int *array, int n);

double rgamma1(double alpha, double beta);
/* Generate Gamma distribution
 * See Monte Carlo Statistical Methods (C. Robert and G. Casella 2004) Page 66
 * alpha: the shape parameter
 * beta: the scale parameter
 * the density of gamma distribution is f(x;alpha,beta) = x^{alpha-1}\frac{e^{-x/beta}}{beta^{alpha}Gamma(alpha)}
 * return -1, if error and print an error message
 * */

double rbeta(double alpha, double beta);
/*generate random number from beta distribution*/
//double rbeta_mv(double mean, double var);
/*Given the mean the the variation of a beta distribution, generate random variable from it*/

int rDirichlet(double *alpha, double *x,int n);
/* Generate random variates from  Dirichlet distribution with parameter alpha,
 * both alpha and x shoulb be of length n
 * if error, return -1;
 * The gererated value will be stored in x;
 * */

double rpois(const double lambda); /*generate random variable with mean lambda*/

double rnegbinom(const double n, const double p); /*generate random number from NegBinom(n,p), n is the pre-specified number of failures*/

double rnegbinom_mv(const double m, const double v); /*use mean and variance as parameters to generate random variable from negative binomial distribution*/

double rbinom(const double pp, const int n); /*generate random variable from Binom(n, pp)*/

double rnorm(const double mu, const double sd); /*generate random variable from normal distribution with mean mu and standard deviation sd*/

#endif
