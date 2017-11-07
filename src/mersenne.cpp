/**************************   mersenne.cpp   **********************************
 * Author:        Agner Fog
 * Date created:  2001
 * Last modified: 2008-11-16
 * Project:       randomc.h
 * Platform:      Any C++
 * Description:
 * Random Number generator of type 'Mersenne Twister'
 *
 * This random number generator is described in the article by
 * M. Matsumoto & T. Nishimura, in:
 * ACM Transactions on Modeling and Computer Simulation,
 * vol. 8, no. 1, 1998, pp. 3-30.
 * Details on the initialization scheme can be found at
 * http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/emt.html
 *
 * Further documentation:
 * The file ran-instructions.pdf contains further documentation and
 * instructions.
 *
 * Copyright 2001-2008 by Agner Fog.
 * GNU General Public License http://www.gnu.org/licenses/gpl.html
 *******************************************************************************/

#include "randomc.h"
#include <cmath>
#include <iostream>



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


int CRandomMersenne::IRandom(int min, int max) {
    // Output random integer in the interval min <= x <= max
    // Relative error on frequencies < 2^-32
    if (max <= min) {
        if (max == min) return min; else return 0x80000000;
    }
    // Multiply interval with random and truncate
    int r = int((double)(uint32_t)(max - min + 1) * Random() + min);
    if (r > max) r = max;
    return r;
}


/***********************************************************************
 Poisson distribution
 ***********************************************************************/
int CRandomMersenne::Poisson (double L) {
    /*
     This function generates a random variate with the poisson distribution.

     Uses inversion by chop-down method for L < 17, and ratio-of-uniforms
     method for L >= 17.

     For L < 1.E-6 numerical inaccuracy is avoided by direct calculation.
     */

    //------------------------------------------------------------------
    //                 choose method
    //------------------------------------------------------------------
    if (L < 17) {
        if (L < 1.E-6) {
            if (L == 0) return 0;

            //--------------------------------------------------------------
            // calculate probabilities
            //--------------------------------------------------------------
            // For extremely small L we calculate the probabilities of x = 1
            // and x = 2 (ignoring higher x). The reason for using this
            // method is to prevent numerical inaccuracies in other methods.
            //--------------------------------------------------------------
            return PoissonLow(L);
        }
        else {
            //--------------------------------------------------------------
            // inversion method
            //--------------------------------------------------------------
            // The computation time for this method grows with L.
            // Gives overflow for L > 80
            //--------------------------------------------------------------
            return PoissonInver(L);
        }
    } else {

        //----------------------------------------------------------------
        // ratio-of-uniforms method
        //----------------------------------------------------------------
        // The computation time for this method does not depend on L.
        // Use where other methods would be slower.
        //----------------------------------------------------------------
        return PoissonRatioUniforms(L);
    }
}


/***********************************************************************
 Subfunctions used by poisson
 ***********************************************************************/
int CRandomMersenne::PoissonLow(double L) {
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
    if (Random() >= d) return 0;
    r = Random() * d;
    if (r > L * (1.-L)) return 0;
    if (r > 0.5 * L*L * (1.-L)) return 1;
    return 2;
}


int CRandomMersenne::PoissonInver(double L) {
    /*
     This subfunction generates a random variate with the poisson
     distribution using inversion by the chop down method (PIN).

     Execution time grows with L. Gives overflow for L > 80.

     The value of bound must be adjusted to the maximal value of L.
     */
    const int bound = 130;              // safety bound. Must be > L + 8*sqrt(L).
    double r;                           // uniform random number
    double f;                           // function value
    int x;                          // return value

    if (L != pois_L_last) {             // set up
        pois_L_last = L;
        pois_f0 = exp(-L);               // f(0) = probability of x=0
    }
    while (1) {
        r = Random();  x = 0;  f = pois_f0;
        do {                             // recursive calculation: f(x) = f(x-1) * L / x
            r -= f;
            if (r <= 0) return x;
            x++;
            f *= L;
            r *= x;                       // instead of f /= x
        }
        while (x <= bound);
    }
}

/***********************************************************************
 Log factorial function
 ***********************************************************************/
double CRandomMersenne::LnFac(int n) {
    // log factorial function. gives natural logarithm of n!

    // define constants
    static const double                 // coefficients in Stirling approximation
    C0 =  0.918938533204672722,      // ln(sqrt(2*pi))
    C1 =  1./12.,
    C3 = -1./360.;
    // C5 =  1./1260.,                  // use r^5 term if FAK_LEN < 50
    // C7 = -1./1680.;                  // use r^7 term if FAK_LEN < 20
    // static variables
    static double fac_table[FAK_LEN];   // table of ln(n!):
    static int initialized = 0;         // remember if fac_table has been initialized

    if (n < FAK_LEN) {
        if (n <= 1) {
            return 0;
        }
        if (!initialized) {              // first time. Must initialize table
            // make table of ln(n!)
            double sum = fac_table[0] = 0.;
            for (int i=1; i<FAK_LEN; i++) {
                sum += log(double(i));
                fac_table[i] = sum;
            }
            initialized = 1;
        }
        return fac_table[n];
    }
    // not found in table. use Stirling approximation
    double  n1, r;
    n1 = n;  r  = 1. / n1;
    return (n1 + 0.5)*log(n1) - n1 + C0 + r*(C1 + r*r*C3);
}


int CRandomMersenne::PoissonRatioUniforms(double L) {
    /*
     This subfunction generates a random variate with the poisson
     distribution using the ratio-of-uniforms rejection method (PRUAt).

     Execution time does not depend on L, except that it matters whether L
     is within the range where ln(n!) is tabulated.

     Reference: E. Stadlober: "The ratio of uniforms approach for generating
     discrete random variates". Journal of Computational and Applied Mathematics,
     vol. 31, no. 1, 1990, pp. 181-189.
     */
    double u;                                          // uniform random
    double lf;                                         // ln(f(x))
    double x;                                          // real sample
    int k;                                             // integer sample

    if (pois_L_last != L) {
        pois_L_last = L;                                // Set-up
        pois_a = L + 0.5;                               // hat center
        int mode = (int)L;                              // mode
        pois_g  = log(L);
        pois_f0 = mode * pois_g - LnFac(mode);          // value at mode
        pois_h = sqrt(SHAT1 * (L+0.5)) + SHAT2;         // hat width
        pois_bound = (int)(pois_a + 6.0 * pois_h);      // safety-bound
    }
    while(1) {
        u = Random();
        if (u == 0) continue;                           // avoid division by 0
        x = pois_a + pois_h * (Random() - 0.5) / u;
        if (x < 0 || x >= pois_bound) continue;         // reject if outside valid range
        k = (int)(x);
        lf = k * pois_g - LnFac(k) - pois_f0;
        if (lf >= u * (4.0 - u) - 3.0) break;           // quick acceptance
        if (u * (u - lf) > 1.0) continue;               // quick rejection
        if (2.0 * log(u) <= lf) break;                  // final acceptance
    }
    return k;
}


CRandomMersenne rndgen(5); //the one random number generator

void set_seed(int seed)
{
    rndgen.RandomInit(seed);
}

double uniform()
{
    return rndgen.Random();
}

int random_number(int n)
{
    return rndgen.IRandom(0,n-1);
}

double poisson(double lambda)
{
    return rndgen.Poisson(lambda);
}






