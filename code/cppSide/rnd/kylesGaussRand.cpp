//  Made by Kyle Mede, code copied from files in /rnd folder
//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

#include "kylesGaussRand.h"

uint32_t CRandomMersenne::BRandom() {
   // Generate 32 random bits
   uint32_t y;

   if (mti >= MERS_N) {
      // Generate MERS_N words at one time
      const uint32_t LOWER_MASK = (1LU << MERS_R) - 1;       // Lower MERS_R bits
      const uint32_t UPPER_MASK = 0xFFFFFFFF << MERS_R;      // Upper (32 - MERS_R) bits
      static const uint32_t mag01[2] = {0, MERS_A};

      int kk;
      for (kk=0; kk < MERS_N-MERS_M; kk++) {
         y = (mt[kk] & UPPER_MASK) | (mt[kk+1] & LOWER_MASK);
         mt[kk] = mt[kk+MERS_M] ^ (y >> 1) ^ mag01[y & 1];}

      for (; kk < MERS_N-1; kk++) {
         y = (mt[kk] & UPPER_MASK) | (mt[kk+1] & LOWER_MASK);
         mt[kk] = mt[kk+(MERS_M-MERS_N)] ^ (y >> 1) ^ mag01[y & 1];}

      y = (mt[MERS_N-1] & UPPER_MASK) | (mt[0] & LOWER_MASK);
      mt[MERS_N-1] = mt[MERS_M-1] ^ (y >> 1) ^ mag01[y & 1];
      mti = 0;
   }
   y = mt[mti++];

   // Tempering (May be omitted):
   y ^=  y >> MERS_U;
   y ^= (y << MERS_S) & MERS_B;
   y ^= (y << MERS_T) & MERS_C;
   y ^=  y >> MERS_L;

   return y;
}


double CRandomMersenne::Random() {
   // Output random float number in the interval 0 <= x < 1
   // Multiply by 2^(-32)
   return (double)BRandom() * (1./(65536.*65536.));
}

void CRandomMersenne::Init0(int seed) {
   // Seed generator
   const uint32_t factor = 1812433253UL;
   mt[0]= seed;
   for (mti=1; mti < MERS_N; mti++) {
      mt[mti] = (factor * (mt[mti-1] ^ (mt[mti-1] >> 30)) + mti);
   }
}

void CRandomMersenne::RandomInit(int seed) {
   // Initialize and seed
   Init0(seed);

   // Randomize some more
   for (int i = 0; i < 37; i++) BRandom();
}

//  END OF CODE FROM MARSENNE.CPP

/***********************************************************************
Constructor
***********************************************************************/
StochasticLib1::StochasticLib1 (int seed)
: STOC_BASE(seed) {
   // Initialize variables for various distributions
   normal_x2_valid = 0;
}

/***********************************************************************
Normal distribution
***********************************************************************/

double StochasticLib1::Normal(double m, double s) {
   // normal distribution with mean m and standard deviation s
   double normal_x1;                   // first random coordinate (normal_x2 is member of class)
   double w;                           // radius
   if (normal_x2_valid) {              // we have a valid result from last call
      normal_x2_valid = 0;
      return normal_x2 * s + m;
   }
   // make two normally distributed variates by Box-Muller transformation
   do {
      normal_x1 = 2. * Random() - 1.;
      normal_x2 = 2. * Random() - 1.;
      w = normal_x1*normal_x1 + normal_x2*normal_x2;
   }
   while (w >= 1. || w < 1E-30);
   w = sqrt(log(w)*(-2./w));
   normal_x1 *= w;  normal_x2 *= w;    // normal_x1 and normal_x2 are independent normally distributed variates
   normal_x2_valid = 1;                // save normal_x2 for next call
   return normal_x1 * s + m;
}

double StochasticLib1::NormalTrunc(double m, double s, double limit) {
   // Truncated normal distribution
   // The tails are cut off so that the output
   // is in the interval from (m-limit) to (m+limit)
   //if (limit < s) FatalError("limit out of range in NormalTrunc function");
   double x;
   do {
      x = Normal(0., s);
   } while (fabs(x) > limit); // reject if beyond limit
   return x + m;
}





