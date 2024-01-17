#ifndef _MONTBIG
#define _MONTBIG

#include "bigint.h"
#include "int.h"

using namespace std;

// montgomery multiplication for multiprecision
// https://en.wikipedia.org/wiki/Montgomery_modular_multiplication

class MontgomeryBig
{
  public:
  static const int k=128;  // power of 2 for R  R=2^k or 1<<k
  static const uint128 B=((uint128)1)<<64; 
  static const uint128 Bminus1=B-1; 
  static uint128 N;
  static uint128 Nprime; // N*Nprime = -1 mod B
  static uint128 BmodN;  // B=2^64
  static uint128 R2modN;  // R^2 mod N
  static uint64_t r2modn[2], n[2];

  inline static void setup(uint128 modulus)
  {
    N=modulus;
    BmodN=B%N; // likely just B if N>B (common)
    Nprime=inv128(N,B);
    Nprime=B-Nprime; // switch sign

    R2modN=BmodN;  // start with 64 bits
    for(int i=65; i<=256; i++)
    {
      R2modN <<=1;
      if(R2modN>N) R2modN-=N;
    }

    topair(r2modn,R2modN);
    topair(n,N);
/*
    cout << "N="<<N<<" B="<<B<<" B%N="<<BmodN<<" B-1="<<Bminus1
	    <<" Nprime="<<Nprime
	    <<" R2modN="<<R2modN << endl;
*/
  }

  inline static void topair(uint64_t y[], uint128 x)
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

  /*
  static void multiply(uint64_t ans[], uint128 X, uint128 Y)
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
  */

  static inline void 
	  multiply(uint64_t xy[/*3*/], uint64_t x, uint64_t y[/*2*/])
//  multiplies a scalar (int64) by a vector length 2
  {
    uint64_t carry=0;
    uint128 temp;
    temp=((uint128)x)*y[0];
    xy[0]=temp&Bminus1;
    carry=temp>>64;

    temp=carry+((uint128)x)*y[1];
    xy[1]=temp&Bminus1;
    xy[2]=temp>>64;
  }
  static inline void addto(uint64_t y[/*3*/], uint64_t x[/*3*/])
// does y=y+x
  {
    uint128 temp;
    uint64_t carry=0;

    temp=((uint128)x[0])+y[0];
    y[0]=temp&Bminus1;
    carry=temp>>64;

    temp=((uint128)x[1])+y[1]+carry;
    y[1]=temp&Bminus1;
    carry=temp>>64;

    temp=((uint128)x[2])+y[2]+carry;
    y[2]=temp&Bminus1;
    //carry=temp>>64;  not needed
  }

  static inline void REDC( uint64_t ans[], uint64_t a[], uint64_t b[] )
  {
    uint64_t q, c[3],ab[3],nq[3];
    //for(int i=0; i<3; i++) c[i]=0; 
    //c[0]=c[1]=c[2]=0;
    //for(int i=0; i<2; i++)
    {
      int i=0;
      //multiply(ab,a[i],b);
      //addto(c,ab);
      multiply(c,a[i],b);
      q =( ((int128)Nprime)*c[0]) & Bminus1;
      multiply(nq,q,n);
      addto(c,nq);
      // C >>= 64;
      c[0]=c[1]; c[1]=c[2]; c[2]=0;

      i=1;
      multiply(ab,a[i],b);
      addto(c,ab);
      q =( ((int128)Nprime)*c[0]) & Bminus1;
      multiply(nq,q,n);
      addto(c,nq);
      // C >>= 64;
      c[0]=c[1]; c[1]=c[2]; c[2]=0;
    }
    uint128 C=frompair(c);
    if(C>=N) C=C-N;
    topair(ans,C);
  }

  static inline void toMont(uint64_t ans[/*2*/], uint128 X)
  {
    uint64_t x[2];
    topair(x,X);
    REDC(ans,x,r2modn);
  }
  static inline uint128 fromMont(uint64_t x[/*2*/])
  {
    uint64_t ans[2], one[2]={1,0};
    REDC(ans,x,one);
    return frompair(ans);
  }

  static uint128 powmod(uint128 xin, uint128 e)
  // uses Montgomery form -- assumes modulus=N is odd
  // assumes setup has been called, but xin is Not in montgomery form
  {
    uint64_t x[2], y[2];
    toMont(x,xin);
    toMont(y,1);
    while(e>0)
    {
      if(e&1==1) REDC(y,y,x); // if e is odd
      e >>= 1;  // e = e/2
      REDC(x,x,x);
    }
    return fromMont(y);
  }

  static inline uint128 powmod(uint128 xin, uint128 e, uint128 m)
  {
    setup(m);
    return powmod(xin,e);
  }

};

// these puppies must be defined - C++ quirk, IMO

  uint128 MontgomeryBig::N;
  uint128 MontgomeryBig::Nprime; // N*Nprime = -1 mod B
  uint128 MontgomeryBig::BmodN;  // B=2^64
  uint128 MontgomeryBig::R2modN;  // R^2 mod N
  uint64_t MontgomeryBig::r2modn[2];
  uint64_t MontgomeryBig::n[2];

#endif
