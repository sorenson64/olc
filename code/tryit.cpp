#include<iostream>
#include<gmp.h>
#include"gmputil.h"
#include "int.h"

using namespace std;

// convert an int128 to a mpz_t
//
void bi2mpz(mpz_t & y, uint128 x)
{
  mpz_set_ui(y,(uint64_t)( x>>64));
  cout << (x>>64) << endl;
  cout << y << endl;
  mpz_mul_2exp(y,y,64);
  cout << y << endl;
  mpz_add_ui(y,y,(uint64_t)(x & 0xFFFFFFFFFFFFFFFF));
}

int main()
{
  char s[]= "400000000148309773081";
  uint128 x=atobi(s);
  //x<<=64;
  cout << "int 128: " << x << endl;

 // cout << sizeof(unsigned long) << endl; // 8 chars, or 64 bits

  mpz_t y;
  mpz_init(y);

  bi2mpz(y,x);

  cout << "mpz: " << y << endl;
  return 0;
}
