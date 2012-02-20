/*********************************
 *
 *  Single Molecule Diffusion Simulator
 *    Random Number Generation (Uniform and Poisson destributions)
 *    Code adapted from the GPL routines available at:
 *      http://www.agner.org/random.html (Agner Fog)
 *    The uniform distribution algorithm is 'Mersenne Twister'.
 *
 *    $Id: random.c,v 1.2 2008/09/26 20:57:50 mculbert Exp $
 *
 *********************************/

#include "random.h"
#include <stdlib.h>
#include <math.h>
#include <time.h>

#ifdef __WIN32__
//int errno;
#endif

#define M 175
#define R 19
#define TEMU 11
#define TEMS 7
#define TEMT 15
#define TEML 17
#define MATRIX_A 0xE4BD75F5
#define TEMB     0x655E5280
#define TEMC     0xFFD58000

static int PoissonInver(double L, RandomStream *rs);	   // poisson by inversion
static int PoissonRatioUniforms(double L, RandomStream *rs);  // pois by ratio of uniforms
static int PoissonLow(double L, RandomStream *rs);		   // poisson for extremely low L
static double LnFac(int n);		/* ln(n!) */

#define FAK_LEN 1024        /* length of factorial table used by LnFac */
static const double SHAT1 = 2.943035529371538573;    /* 8/e */
static const double SHAT2 = 0.8989161620588987408;   /* 3-sqrt(12/e) */

static RandomStream internalStream;

void seedRnd(int seed, RandomStream *rs)
{
  unsigned int s = (unsigned int)seed;
  if (rs == NULL) rs = &internalStream;
  for (rs->mti=0; rs->mti < N; rs->mti++)
  {
    s = s * 29943829 - 1;
    rs->mt[rs->mti] = s;
  }
  rs->haveNormal = 0;
//  srandom(seed+Random());
}

void initRnd() { seedRnd(time(0), NULL); }

unsigned int Random(RandomStream *rs)
{
  /* generate 32 random bits */
  unsigned int y;
    const unsigned int LOWER_MASK = (1u << R) - 1; /* lower R bits */
    const unsigned int UPPER_MASK = -1 << R;       /* upper 32-R bits */
    int kk, km;
  if (rs == NULL) rs = &internalStream;

  if (rs->mti >= N) {
    /* generate N words at one time */
    for (kk=0, km=M; kk < N-1; kk++) {
      y = (rs->mt[kk] & UPPER_MASK) | (rs->mt[kk+1] & LOWER_MASK);
      rs->mt[kk] = rs->mt[km] ^ (y >> 1) ^ (-(signed int)(y & 1) & MATRIX_A);
      if (++km >= N) km = 0;}

    y = (rs->mt[N-1] & UPPER_MASK) | (rs->mt[0] & LOWER_MASK);
    rs->mt[N-1] = rs->mt[M-1] ^ (y >> 1) ^ (-(signed int)(y & 1) & MATRIX_A);
    rs->mti = 0;}

  y = rs->mt[rs->mti++];

  /* Tempering (May be omitted): */
  y ^=  y >> TEMU;
  y ^= (y << TEMS) & TEMB;
  y ^= (y << TEMT) & TEMC;
  y ^=  y >> TEML;

  return y;                                                          
}

double DRandom(RandomStream *rs)
{
  /* output random float number in the interval 0 <= x < 1 */
  union {
    double f;
    unsigned int i[2];}
  convert;
  /* get 32 random bits and convert to float */
  unsigned int r = Random(rs);
  convert.i[0] =  r << 20;
  convert.i[1] = (r >> 12) | 0x3FF00000;
  return convert.f - 1.0;
}

double Normal(double mean, double sd, RandomStream *rs)
{
  double x1, x2, w;
  if (rs == NULL) rs = &internalStream;
  if (rs->haveNormal)
  {
    rs->haveNormal = 0;
    return rs->normal*sd + mean;
  }
  do
  {
    x1 = 2 * DRandom(rs) - 1;
    x2 = 2 * DRandom(rs) - 1;
    w = x1*x1 + x2*x2;
  } while (w >= 1);
  w = sqrt( -2*log(w)/w );
  rs->haveNormal = 1;
  rs->normal = x1*w;
  return x2*w*sd + mean;
}

unsigned int Poisson(double L, RandomStream *rs)
{
/*
   This function generates a random variate with the poisson distribution.

   Uses inversion by chop-down method for L < 17, and ratio-of-uniforms
   method for L >= 17.

   For L < 1.E-6 numerical inaccuracy is avoided by direct calculation.
*/
 
  /*------------------------------------------------------------------
  //                 choose method
  //------------------------------------------------------------------*/
  if (L < 17) {
    if (L < 1.E-6) {
      if (L == 0) return 0;
      if (L < 0) return 0; 
	/* FatalError("Parameter negative in poisson function"); */
    
      /*--------------------------------------------------------------
      // calculate probabilities
      //--------------------------------------------------------------
      // For extremely small L we calculate the probabilities of x = 1
      // and x = 2 (ignoring higher x). The reason for using this 
      // method is to prevent numerical inaccuracies in other methods.
      //--------------------------------------------------------------*/
      return PoissonLow(L, rs);}
    
    else {
    
      /*--------------------------------------------------------------
      // inversion method
      //--------------------------------------------------------------
      // The computation time for this method grows with L.
      // Gives overflow for L > 80
      //--------------------------------------------------------------*/
      return PoissonInver(L, rs);}}
      
  else {
    if (L > 2.E9) return 0;
	/* FatalError("Parameter too big in poisson function"); */

    /*----------------------------------------------------------------
    // ratio-of-uniforms method
    //----------------------------------------------------------------
    // The computation time for this method does not depend on L.
    // Use where other methods would be slower.
    //----------------------------------------------------------------*/
    return PoissonRatioUniforms(L, rs);}
}

int PoissonLow(double L, RandomStream *rs) {
/*
   This subfunction generates a random variate with the poisson 
   distribution for extremely low values of L.

   The method is a simple calculation of the probabilities of x = 1
   and x = 2. Higher values are ignored.

   The reason for using this method is to avoid the numerical inaccuracies 
   in other methods.
*/   
  double d, r;
  d = sqrt(L);
  if (DRandom(rs) >= d) return 0;
  r = DRandom(rs) * d;
  if (r > L * (1.-L)) return 0;
  if (r > 0.5 * L*L * (1.-L)) return 1;
  return 2;}

int PoissonInver(double L, RandomStream *rs) {
/*
   This subfunction generates a random variate with the poisson 
   distribution using inversion by the chop down method (PIN).

   Execution time grows with L. Gives overflow for L > 80.

   The value of bound must be adjusted to the maximal value of L.
*/   
  const int bound = 130;          /* safety bound. Must be > L + 8*sqrt(L). */
  double p_f0 = exp(-L);             /* value at x=0 */
  double r;                          /* uniform random number */
  double f;                          /* function value */
  int x;                        /* return value */

  while (1) {  
    r = DRandom(rs);  x = 0;  f = p_f0;
    do {                    /* recursive calculation: f(x) = f(x-1) * L / x */
      r -= f;
      if (r <= 0) return x;
      x++;
      f *= L;
      r *= x;}                       /* instead of f /= x */
    while (x <= bound);}}

int PoissonRatioUniforms(double L, RandomStream *rs) {
/*
   This subfunction generates a random variate with the poisson 
   distribution using the ratio-of-uniforms rejection method (PRUAt).

   Execution time does not depend on L, except that it matters whether L
   is within the range where ln(n!) is tabulated.

   Reference: E. Stadlober: "The ratio of uniforms approach for generating
   discrete random variates". Journal of Computational and Applied Mathematics,
   vol. 31, no. 1, 1990, pp. 181-189.
*/
  double p_a;                              /* hat center */
  double p_h;                              /* hat width */
  double p_g;                              /* ln(L) */
  double p_q;                              /* value at mode */
  int p_bound;                        /* upper bound */
  int mode;                           /* mode */
  double u;                                /* uniform random */
  double lf;                               /* ln(f(x)) */
  double x;                                /* real sample */
  int k;                              /* integer sample */

    p_a = L + 0.5;                          /* hat center */
    mode = (int)L;                     /* mode */
    p_g  = log(L);
    p_q = mode * p_g - LnFac(mode);         /* value at mode */
    p_h = sqrt(SHAT1 * (L+0.5)) + SHAT2;    /* hat width */
    p_bound = (int)(p_a + 6.0 * p_h);  /* safety-bound */

  while(1) {
    u = DRandom(rs);
    if (u == 0) continue;                   /* avoid division by 0 */
    x = p_a + p_h * (DRandom(rs) - 0.5) / u;
    if (x < 0 || x >= p_bound) continue;  /* reject if outside valid range */
    k = (int)(x);
    lf = k * p_g - LnFac(k) - p_q;
    if (lf >= u * (4.0 - u) - 3.0) break;   /* quick acceptance */
    if (u * (u - lf) > 1.0) continue;       /* quick rejection */
    if (2.0 * log(u) <= lf) break;}         /* final acceptance */
  return(k);}

double LnFac(int n) {
  /* log factorial function. gives natural logarithm of n! */

  double sum, n1, r;
  int i;

  static double fac_table[FAK_LEN];   /* table of ln(n!) */
  static int fac_table_initialized = 0;

  /* define constants */
  static const double        /* coefficients in Stirling approximation   */
    C0 =  0.918938533204672722,   /* ln(sqrt(2*pi)) */
    C1 =  1./12., 
    C3 = -1./360.;
    /* C5 =  1./1260.,  // use r^5 term if FAK_LEN < 50 */
    /* C7 = -1./1680.;  // use r^7 term if FAK_LEN < 20 */

  if (n < FAK_LEN) {
    if (n <= 1) {
      /* if (n < 0) FatalError("Parameter negative in LnFac function"); */
      return 0;}
    if (!fac_table_initialized) { /* first time. Must initialize table */
      /* make table of ln(n!) */
      sum = fac_table[0] = 0.;
      for (i=1; i<FAK_LEN; i++) {
        sum += log(i);
        fac_table[i] = sum;}
      fac_table_initialized = 1;}
    return fac_table[n];}
    
  /* not found in table. use Stirling approximation */
  n1 = n;  r  = 1. / n1;
  return (n1 + 0.5)*log(n1) - n1 + C0 + r*(C1 + r*r*C3);}
