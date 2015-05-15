//  Made by Kyle Mede, code copied from files in /rnd folder
//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

#include <math.h>
#include "randomc.h"

// STOC_BASE defines which base class to use for the non-uniform
// random number generator classes StochasticLib1, 2, and 3.
#define STOC_BASE CRandomMersenne     // define random number generator base class

/***********************************************************************
         Class StochasticLib1
***********************************************************************/

class StochasticLib1 : public STOC_BASE {
   // This class encapsulates the random variate generating functions.
   // May be derived from any of the random number generators.
public:
   StochasticLib1 (int seed);          // Constructor
   double Normal(double m, double s);  // Normal distribution
   double NormalTrunc(double m, double s, double limit); // Truncated normal distribution

   // functions used internally
protected:
   // Variables specific to each distribution:
   // Variables used by Normal distribution
   double normal_x2;  int normal_x2_valid;
};

