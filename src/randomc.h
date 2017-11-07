/*****************************   randomc.h   **********************************
 * Author:        Agner Fog
 * Date created:  1997
 * Last modified: 2008-11-16
 * Project:       randomc.h
 * Source URL:    www.agner.org/random
 *
 * Description:
 * This header file contains class declarations and other definitions for the
 * randomc class library of uniform random number generators in C++ language.
 *
 * Overview of classes:
 * ====================
 *
 * class CRandomMersenne:
 * Random number generator of type Mersenne twister.
 * Source file mersenne.cpp
 *
 * class CRandomMother:
 * Random number generator of type Mother-of-All (Multiply with carry).
 * Source file mother.cpp
 *
 * class CRandomSFMT:
 * Random number generator of type SIMD-oriented Fast Mersenne Twister.
 * The class definition is not included here because it is not
 * portable to all platforms. See sfmt.h and sfmt.cpp for details.
 *
 * Member functions (methods):
 * ===========================
 *
 * All these classes have identical member functions:
 *
 * Constructor(int seed):
 * The seed can be any integer. The time may be used as seed.
 * Executing a program twice with the same seed will give the same sequence
 * of random numbers. A different seed will give a different sequence.
 *
 * void RandomInit(int seed);
 * Re-initializes the random number generator with a new seed.
 *
 * void RandomInitByArray(int const seeds[], int NumSeeds);
 * In CRandomMersenne and CRandomSFMT only: Use this function if you want
 * to initialize with a seed with more than 32 bits. All bits in the seeds[]
 * array will influence the sequence of random numbers generated. NumSeeds
 * is the number of entries in the seeds[] array.
 *
 * double Random();
 * Gives a floating point random number in the interval 0 <= x < 1.
 * The resolution is 32 bits in CRandomMother and CRandomMersenne, and
 * 52 bits in CRandomSFMT.
 *
 * int IRandom(int min, int max);
 * Gives an integer random number in the interval min <= x <= max.
 * (max-min < MAXINT).
 * The precision is 2^-32 (defined as the difference in frequency between
 * possible output values). The frequencies are exact if max-min+1 is a
 * power of 2.
 *
 * int IRandomX(int min, int max);
 * Same as IRandom, but exact. In CRandomMersenne and CRandomSFMT only.
 * The frequencies of all output values are exactly the same for an
 * infinitely long sequence. (Only relevant for extremely long sequences).
 *
 * uint32_t BRandom();
 * Gives 32 random bits.
 *
 *
 * Example:
 * ========
 * The file EX-RAN.CPP contains an example of how to generate random numbers.
 *
 *
 * Library version:
 * ================
 * Optimized versions of these random number generators are provided as function
 * libraries in randoma.zip. These function libraries are coded in assembly
 * language and support only x86 platforms, including 32-bit and 64-bit
 * Windows, Linux, BSD, Mac OS-X (Intel based). Use randoma.h from randoma.zip
 *
 *
 * Non-uniform random number generators:
 * =====================================
 * Random number generators with various non-uniform distributions are
 * available in stocc.zip (www.agner.org/random).
 *
 *
 * Further documentation:
 * ======================
 * The file ran-instructions.pdf contains further documentation and
 * instructions for these random number generators.
 *
 * Copyright 1997-2008 by Agner Fog.
 * GNU General Public License http://www.gnu.org/licenses/gpl.html
 *******************************************************************************/

#ifndef FINITE_CHROMOSOME_ALWAYS_RECOM_RANDOMC_H_
#define FINITE_CHROMOSOME_ALWAYS_RECOM_RANDOMC_H_

// Define integer types with known size: int32_t, uint32_t, int64_t, uint64_t.
// If this doesn't work then insert compiler-specific definitions here:
#if defined(__GNUC__) || (defined(_MSC_VER) && _MSC_VER >= 1600)
// Compilers supporting C99 or C++0x have stdint.h defining these integer types
#include <stdint.h>
#define INT64_SUPPORTED  // Remove this if the compiler doesn't support 64-bit integers
#elif defined(_WIN16) || defined(__MSDOS__) || defined(_MSDOS)
// 16 bit systems use long int for 32 bit integer.
typedef   signed long int int32_t;
typedef unsigned long int uint32_t;
#elif defined(_MSC_VER)
// Older Microsoft compilers have their own definition
typedef   signed __int32  int32_t;
typedef unsigned __int32 uint32_t;
typedef   signed __int64  int64_t;
typedef unsigned __int64 uint64_t;
#define INT64_SUPPORTED  // Remove this if the compiler doesn't support 64-bit integers
#else
// This works with most compilers
typedef signed int          int32_t;
typedef unsigned int       uint32_t;
typedef long long           int64_t;
typedef unsigned long long uint64_t;
#define INT64_SUPPORTED  // Remove this if the compiler doesn't support 64-bit integers
#endif


/***********************************************************************
 System-specific user interface functions
 ***********************************************************************/

//void EndOfProgram(void);                 // System-specific exit code (userintf.cpp)

//void FatalError(const char *ErrorText);  // System-specific error reporting (userintf.cpp)

#include <iostream>

/***********************************************************************
 Define random number generator classes
 ***********************************************************************/

class CRandomMersenne {                // Encapsulate random number generator
    // Choose which version of Mersenne Twister you want:
#if 0
    // Define constants for type MT11213A:
#define MERS_N   351
#define MERS_M   175
#define MERS_R   19
#define MERS_U   11
#define MERS_S   7
#define MERS_T   15
#define MERS_L   17
#define MERS_A   0xE4BD75F5
#define MERS_B   0x655E5280
#define MERS_C   0xFFD58000
#else
    // or constants for type MT19937:
#define MERS_N   624
#define MERS_M   397
#define MERS_R   31
#define MERS_U   11
#define MERS_S   7
#define MERS_T   15
#define MERS_L   18
#define MERS_A   0x9908B0DF
#define MERS_B   0x9D2C5680
#define MERS_C   0xEFC60000
#endif

public:
    CRandomMersenne(int seed) {         // Constructor
        RandomInit(seed); LastInterval = 0; normal_x2_valid = 0;
    }

    void RandomInit(int seed);          // Re-seed
    void RandomInitByArray(int const seeds[], int NumSeeds);  // Seed by more than 32 bits
    int IRandom(int min, int max);      // Output random integer
    int IRandomX(int min, int max);     // Output random integer, exact
    double Random();                    // Output random float
    uint32_t BRandom();                 // Output random bits

    double normal(double m, double s);  // included from stocc package
    int Poisson (double L);             // included from stocc package

    double LnFac(int n);
    // Variables used by Poisson distribution
    double pois_L_last;                                     // previous value of L
    double pois_f0;                                         // value at x=0 or at mode
    double pois_a;                                          // hat center
    double pois_h;                                          // hat width
    double pois_g;                                          // ln(L)
    int  pois_bound;                                        // upper bound
    // constant for LnFac function:
    static const int FAK_LEN = 1024;       // length of factorial table
    const double SHAT1 = 2.943035529371538573;    // 8/e
    const double SHAT2 = 0.8989161620588987408;   // 3-sqrt(12/e)

private:
    void Init0(int seed);               // Basic initialization procedure
    uint32_t mt[MERS_N];                // State vector
    int mti;                            // Index into mt
    uint32_t LastInterval;              // Last interval length for IRandomX
    uint32_t RLimit;                    // Rejection limit used by IRandomX
    double normal_x2;
    bool normal_x2_valid;

    int PoissonLow(double L);           // Poisson helper function
    int PoissonInver(double L);         // Poisson helper function
    int PoissonRatioUniforms(double L); // Poisson helper function


};

void set_seed(int seed);
double uniform();
int random_number(int n);
double normal(double m, double s);
double poisson(double lambda);

#endif  // FINITE_CHROMOSOME_ALWAYS_RECOM_RANDOMC_H_
