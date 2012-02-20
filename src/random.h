/********************************************
 *
 *  Single Molecule Diffusion Simulator
 *    Random Number Distributions (uniform, Poissonian)
 *    $Id: random.h,v 1.2 2008/09/26 20:57:50 mculbert Exp $
 *
 *******************************************/
      

#ifndef __SMDS_RANDOM_H__
#define __SMDS_RANDOM_H__

#define N 351
typedef struct {
  /* private */
  unsigned int mt[N];
  int mti, haveNormal;
  double normal;
} RandomStream;

void initRnd(void);
void seedRnd(int seed, RandomStream *);
unsigned int Random(RandomStream *);
double DRandom(RandomStream *);
double Normal(double mean, double sd, RandomStream *);
unsigned int Poisson(double mean, RandomStream *);

#ifdef __WIN32__
#define srandom srand
#define random rand
#endif

#endif
