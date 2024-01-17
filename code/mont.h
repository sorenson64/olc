#ifndef _MONTBIG
#define _MONTBIG

#include "bigint.h"
#include "uint.h"

using namespace std;

// montgomery multiplication

class MontgomeryBig
{
  public:
  static int k;  // power of 2 for R  R=2^k or 1<<k
  static uint128 N;
  static uint128 Nprime; // N*Nprime = -1 mod R
  static uint128 R, Rminus1, RmodN;  // R=2^k
  static uint128 B, Bminus1, BmodN;  // B=2^64

  static int mink(uint128 x)
  {
    int k=1;
    uint128 one=1;
    while( (one<<k) <= x ) k++;
    return k;
  }

  inline static void setup(uint128 modulus)
  {
    N=modulus;
    k=mink(modulus);
    R=1; R<<=k;
    Rminus1=R-1;
    RmodN=R%N;
    B=1; B<<=64;
    Bminus1=B-1;
    BmodN=B%N;
    Nprime=inv128(N,R);
    Nprime=R-Nprime; // switch sign
  }

  inline static void topair(uint64_t y[2], uint128 x)
  {
    y[0]=x & ((((int128)1)<<64)-1);
    y[1]=x >> 64;
  }
  inline static uint128 frompair(uint64_t x[])
  {
    uint128 y=x[1];
    y <<=64;
    y = y | x[0];
    return y;
  }
  static void multiply(uint64_t ans[4], uint128 X, uint128 Y)
  {
    uint64_t x[2],y[2];  uint128 z[3];
    topair(x,X); topair(y,Y);

    z[0]=((uint128)x[0])*y[0];
    z[1]= ((uint128)x[1])*y[0] + ((uint128)x[0])*y[1]; // Karatsuba?
    z[2]=((uint128)x[1])*y[1];

    z[1]+= z[0]>>64;
    ans[0] = z[0] & Bminus1;

    z[2]+= z[1]>>64;
    ans[1] = z[1] & Bminus1;

    ans[3] = z[2]>>64;
    ans[2] = z[2] & Bminus1;
  }

  static uint128 multmodR(uint128 X, uint128 Y)
  {
    uint64_t x[2],y[2];  uint128 z[2];
    topair(x,X&Rminus1); topair(y,Y&Rminus1);

    z[0]=( ((uint128)x[0])*y[0] )&Rminus1;
    z[1]=(  ((uint128)x[1])*y[0] + ((uint128)x[0])*y[1] )&Bminus1;
    // z[2] not needed - shifted by 128 bits so divisible by R

    return (z[0]+(z[1]<<64))&Rminus1;
  }

  static inline uint128 toMont(uint128 x)
  {
    uint128 y=x; return multmodN(y,RmodN);
  }
  static inline uint128 fromMont(uint128 T)
  {
    // return T Rinv mod N
    uint128 m = multmodR(T,Nprime); // ( (T & Rminus1) * Nprime ) & Rminus1;
    uint128 t = ( T + m*N) >> k;
    return (t>=N)?(t-N):t;
  }
  static inline uint64_t mult(uint64_t a, uint64_t b)
  {
    uint128 biga=a;
    return fromMont(biga*b);
  }

  static uint128 powmod(int128 xin, int128 e)
  // uses Montgomery form -- assumes modulus=N is odd
  // assumes setup has been called, but xin is Not in montgomery form
  {
    uint128 y=toMont(1); 
    uint128 x=toMont(xin);
    while(e>0)
    {
      if(e&1==1) y=mult(y,x); // if e is odd
      e >>= 1;  // e = e/2
      x = mult(x,x);
    }
    return fromMont(y);
  }

  static inline uint128 powmod(int128 xin, int128 e, uint128 m)
  {
    setup(m);
    return powmod(xin,e);
  }

};

// these puppies must be defined - C++ quirk, IMO

int Montgomery::k;  
uint64_t Montgomery::N;
uint64_t Montgomery::Nprime; 
uint128 Montgomery::R; 
uint128 Montgomery::Rminus1;
uint64_t Montgomery::RmodN;

#endif
