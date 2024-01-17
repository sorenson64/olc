#ifndef _GMPPRINT
#define _GMPPRINT

#include "gmp.h"
#include "int.h"

using namespace std;

void gmpprint(ostream &out, const mpz_t &x)
{
  int size=mpz_sizeinbase(x,10);
  char s[size+10];
  mpz_get_str(s,10,x);
  out << s;
}

inline ostream& operator<< (ostream& out, const mpz_t &x)
{ gmpprint(out,x); return out; }


void bigint2mpz(mpz_t & y, uint128 x)
{
  mpz_set_ui(y, x>>64);
  mpz_mul_2exp(y,y,64);
  mpz_add_ui(y,y, x & 0xFFFFFFFFFFFFFFFF);
}

#endif
