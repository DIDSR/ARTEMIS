/*
 * Definition file: Mersenne twister and box muller
 * Created: Yuan Fang
 * Date: 2013-04-18
 *
 * Additional notes:
 * Modification history
 *	1. changed genrand_real1 & genrand_real3 return types from float to double
 *	2. changed genrand_real1 to return [min,max] instead of [0,1]
 *  3. additional definition of box muller method by Everett F. Carter Jr.
 *
 *
 */

#include<math.h>
#include<stdio.h>
#include<stdlib.h>

////////////////////////////////////////////////////////////////////////////
//			  MERSENNE TWISTER				  							  //
////////////////////////////////////////////////////////////////////////////
//
// Random numbers (mersenne twister).
// A C-program for MT19937, with initialization improved 2002/2/10. Coded by 
// Takuji Nishimura and Makoto Matsumoto. This is a faster version by taking  
// Shawn Cokus's optimization, Matthe Bellew's simplification, Isaku Wada's 
// real version.
// Copyright (C) 1997 - 2002, Makoto Matsumoto and Takuji Nishimura, All 
// rights reserved. Redistribution and use in source and binary forms, with 
// or without modification, are permitted provided that the following 
// conditions are met:
//    Redistributions of source code must retain the above copyright notice, 
//    this list of conditions and the following disclaimer. 
//    Redistributions in binary form must reproduce the above copyright 
//    notice, this list of conditions and the following disclaimer in the 
//    documentation and/or other materials provided with the distribution.
//    The names of its contributors may not be used to endorse or promote 
//    products derived from this software without specific prior written 
//    permission.
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE 
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF 
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS 
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN 
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) 
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE 
// POSSIBILITY OF SUCH DAMAGE. Any feedback is very welcome. 
// http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/emt.html 
// email: m-mat @ math.sci.hiroshima-u.ac.jp (remove space)

// Modification history
//	1. changed genrand_real1 & genrand_real3 return types from float to double
//	2. changed genrand_real1 to return [min,max] instead of [0,1]
//
////////////////////////////////////////////////////////////////////////////

 
/* Period parameters */
#define N 624
#define M 397
#define MATRIX_A 0x9908b0dfUL   /* constant vector a */
#define UPPER_MASK 0x80000000UL /* most significant w-r bits */
#define LOWER_MASK 0x7fffffffUL /* least significant r bits */
#define pi 3.1415926

static unsigned long mt[N]; /* the array for the state vector  */
static int mti=N+1; /* mti==N+1 means mt[N] is not initialized */

/* initializes mt[N] with a seed */
void init_genrand(unsigned long s)
{
        mt[0]= s & 0xffffffffUL;
        for (mti=1; mti<N; mti++) {
                mt[mti] = (1812433253UL * (mt[mti-1] ^ (mt[mti-1] >> 30)) + mti);
                /* See Knuth TAOCP Vol2. 3rd Ed. P.106 for multiplier. */
                /* In the previous versions, MSBs of the seed affect   */
                /* only MSBs of the array mt[].                        */
                /* 2002/01/09 modified by Makoto Matsumoto             */
                mt[mti] &= 0xffffffffUL;
                /* for >32 bit machines */
        }
}

/* initialize by an array with array-length */
/* init_key is the array for initializing keys */
/* key_length is its length */
/* slight change for C++, 2004/2/26 */
void init_by_array(unsigned long init_key[], int key_length)
{
    int i, j, k;
    init_genrand(19650218UL);
    i=1; j=0;
    k = (N>key_length ? N : key_length);
    for (; k; k--) {
            mt[i] = (mt[i] ^ ((mt[i-1] ^ (mt[i-1] >> 30)) * 1664525UL))
                    + init_key[j] + j; /* non linear */
            mt[i] &= 0xffffffffUL; /* for WORDSIZE > 32 machines */
            i++; j++;
            if (i>=N) { mt[0] = mt[N-1]; i=1; }
            if (j>=key_length) j=0;
    }
	for (k=N-1; k; k--) {
            mt[i] = (mt[i] ^ ((mt[i-1] ^ (mt[i-1] >> 30)) * 1566083941UL))
                    - i; /* non linear */
            mt[i] &= 0xffffffffUL; /* for WORDSIZE > 32 machines */
            i++;
            if (i>=N) { mt[0] = mt[N-1]; i=1; }
    }

	mt[0] = 0x80000000UL; /* MSB is 1; assuring non-zero initial array */
}

/* generates a random number on [0,0xffffffff]-interval */
unsigned long genrand_int32(void)
{
    unsigned long y;
    static unsigned long mag01[2]={0x0UL, MATRIX_A};
    /* mag01[x] = x * MATRIX_A  for x=0,1 */

    if (mti >= N) { /* generate N words at one time */
            int kk;

            if (mti == N+1)   /* if init_genrand() has not been called, */
                    init_genrand(5489UL); /* a default initial seed is used */

            for (kk=0;kk<N-M;kk++) {
                    y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
                    mt[kk] = mt[kk+M] ^ (y >> 1) ^ mag01[y & 0x1UL];
            }
            for (;kk<N-1;kk++) {
                    y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
                    mt[kk] = mt[kk+(M-N)] ^ (y >> 1) ^ mag01[y & 0x1UL];
            }
            y = (mt[N-1]&UPPER_MASK)|(mt[0]&LOWER_MASK);
            mt[N-1] = mt[M-1] ^ (y >> 1) ^ mag01[y & 0x1UL];

            mti = 0;
    }

	y = mt[mti++];

    /* Tempering */
    y ^= (y >> 11);
    y ^= (y << 7) & 0x9d2c5680UL;
    y ^= (y << 15) & 0xefc60000UL;
    y ^= (y >> 18);

    return y;
}

/* generates a random number on [min,max]-real-interval */
double randMT(double min, double max)
{
	return (double)genrand_int32()*(1.0/4294967295.0)*(max-min)+min;
	/* divided by 2^32-1 */
}

/* generates a random number on (0,1)-real-interval */
double genrand_real3(void)
{
        return (((double)genrand_int32()) + 0.5)*(1.0/4294967296.0);
        /* divided by 2^32 */
}

/*
                      (c) Copyright 1994, Everett F. Carter Jr.
                          Permission is granted by the author to use
			  this software for any application provided this
			  copyright notice is preserved.
*/
double box_muller(double m, double s)	/* normal random variate generator */
{				        /* mean m, standard deviation s */
	double MIN = 1e-20;
	double x1, x2, w, y1;
	static double y2;
	static int use_last = 0;
	if (use_last)		        /* use value from previous call */
	{
		y1 = y2;
		use_last = 0;
	} else {
		do {
			x1 = 2.0 * randMT(MIN,1) - 1.0;
			x2 = 2.0 * randMT(MIN,1) - 1.0;
			w = x1 * x1 + x2 * x2;
		} while ( w >= 1.0 );
		w = sqrt( (-2.0 * log( w ) ) / w );
		y1 = x1 * w;
		y2 = x2 * w;
		use_last = 1;
	}
	return( m + y1 * s );
}

//Returns as a floating-point number an integer value that is a random deviate drawn //from a Poisson distribution of mean xm, using ran1(idum) as a source of uniform
//random deviates.
float poidev(float xm){
	float gammln(float xx);
	//float ran1(long *idum);
	static float sq,alxm,g,oldm=(-1.0); //oldm is a flag for whether xm has changed
	float em,t,y; //since last call.
	if (xm < 12.0) { //Use direct method.
		if (xm != oldm) {
			oldm=xm;
			g=exp(-xm); //If xm is new, compute the exponential.
		}
		em = -1;
		t=1.0;

		//Instead of adding exponential deviates it is equivalent to multiply uniform deviates.
		//We never actually have to take the log, merely compare to the pre-computed exponential.
		do {
			++em;
			t *= randMT(0,1);
		} while (t > g);
	} else { //Use rejection method.
		if (xm != oldm) {
			//If xm has changed since the last call, then pre-compute some functions that occur below.
			oldm=xm;
			sq=sqrt(2.0*xm);
			alxm=log(xm);
			g=xm*alxm-gammln(xm+1.0);
			//The function gammln is the natural log of the gamma function, as given in x6:1.
		}
		do {
			do { //y is a deviate from a Lorentzian comparison function.
				y=tan(pi*randMT(0,1));
				em=sq*y+xm; //em is y, shifted and scaled.
			} while (em < 0.0); //Reject if in regime of zero probability.

			em=floor(em); //The trick for integer-valued distributions.
			t=0.9*(1.0+y*y)*exp(em*alxm-gammln(em+1.0)-g);

			/*	The ratio of the desired distribution to the comparison function; we accept or
				reject by comparing it to another uniform deviate. The factor 0.9 is chosen so
				that t never exceeds 1.  */
		} while (randMT(0,1) > t);
	}
	return em;
}

//Returns the value ln[..(xx)] for xx > 0.
//Internal arithmetic will be done in double precision, a nicety that you can omit if ^Lve-^Lgure
//accuracy is good enough.
float gammln(float xx) {
	double x,y,tmp,ser;
	static double cof[6]={76.18009172947146,-86.50532032941677,24.01409824083091,-1.231739572450155,0.1208650973866179e-2,-0.5395239384953e-5};
	int j;
	y=x=xx;
	tmp=x+5.5;
	tmp -= (x+0.5)*log(tmp);
	ser=1.000000000190015;
	for (j=0;j<=5;j++) ser += cof[j]/++y;
	return -tmp+log(2.5066282746310005*ser/x);
}

