#ifndef _MONT2
#define _MONT2

#include "bigint.h"
#include "uint.h"

using namespace std;

// montgomery multiplication for modulus=p^2

class Montgomery2
{
  public:
  static int k;  // power of 2 for R  R=2^k or 1<<k
  static uint64_t p;  // square root of modulus
  static uint64_t q; // pq=1 mod R
  static uint128 R; // we use this larger size
  //static uint128 Rminus1;
  static uint64_t Rminus1;
  static uint128 Rmodp;

  static int mink(int64_t x)
  {
    int k=1;
    uint128 one=1;
    while( (one<<k) <= x ) k++;
    return k;
  }

  inline static void setup(int kin, uint64_t pin)
  // pin is the square root of the modulus
  {
    p=pin;
    k=kin;
    R=1; R<<=k;
    Rminus1=R-1;
    Rmodp=R%p;
    q=inv(p,R);

/*
cout << p << "=p  k=" << k << " q=" << q <<
   " R="<<R<<" Rminus1="<<Rminus1<<
   " Rmodp=" << Rmodp<< endl;
*/
    // we are assuming here that k<=64.
  }

  static inline void toMont(uint64_t y[2], uint128 xin)
  {
    uint64_t x[2];
    x[0]=xin%p; x[1]=(xin/p)%p;
    y[0]=(Rmodp*x[0])%p;
    uint64_t carry=((R*x[0])/p)%p;
    y[1]=(carry+Rmodp*x[1])%p;

//    cout <<  "should be:" << (xin*R)%p << " " << ((xin*R)/p)%p << endl;
//    cout <<  " but is:  "<< y[0] << " " << y[1] << endl;
  }
  static inline uint128 fromMont(const uint64_t x[])
  {
    const uint64_t one[2]={1,0};
    uint64_t z[2];
    mult(z,x,one);
    return ( ((uint128)(z[1]%p))*p+z[0] );
  }
  static inline void mult(uint64_t zout[], const uint64_t x[], const uint64_t y[])
  // z=xy mod p^2
  {
    uint128 t,u,v,tprime,uprime,vprime;
    int64_t z[2];

/*
cout << "inside mult: x="<< x[0] <<"+p*"<< x[1] << "\t y="
     << y[0] << "+p*"<<y[1] << endl;
*/

    t=x[0]; t=t*y[0];
    //u=q; u=u*(t & Rminus1);
    u=q; u=u*((uint64_t)t & (uint64_t)Rminus1);
    v=p; v=v*((uint64_t)u & Rminus1);
    tprime=((uint128)x[0])*y[1]+((uint128)x[1])*y[0]+(u & Rminus1);
    uprime=q; uprime*=((uint64_t)tprime & Rminus1);
    vprime=p; vprime*=((uint64_t)uprime & Rminus1);

    z[0]=(uint64_t)(t >> k)-(uint64_t)(v >> k);
    z[1]=(uint64_t)(tprime >> k)-(uint64_t)(vprime >> k);
    if(z[0]<0) { z[0]+=p; z[1]--; }
    if(z[1]<0) z[1]+=p;
    else if(z[1]>p) z[1]-=p;
    zout[0]=z[0]; zout[1]=z[1];
//cout << "answer is " << z[0] << "+p*"<<z[1] << endl;
  }

  static uint128 powmod(uint128 xin, uint128 e)
  // uses Montgomery form -- assumes modulus=N is odd
  // assumes setup has been called, but xin is Not in montgomery form
  {
    uint64_t y[2]; toMont(y,1); 
    uint64_t x[2]; toMont(x,xin);
    while(e>0)
    {
      if(e&1==1) mult(y,y,x); // if e is odd
      e >>= 1;  // e = e/2
      mult(x,x,x);
    }
    return fromMont(y);
  }

  static inline uint128 powmod(uint128 xin, uint128 e, uint64_t p)
  {
    setup( mink(p), p);
    return powmod(xin,e);
  }

};

// these puppies must be defined - C++ quirk, IMO

int Montgomery2::k;  
uint64_t Montgomery2::p;
uint64_t Montgomery2::q; 
uint128 Montgomery2::R; 
//uint128 Montgomery2::Rminus1;
uint64_t Montgomery2::Rminus1;
uint128 Montgomery2::Rmodp;

#endif
