#include<iostream>
#include<gmp.h>

#include"montbig.h"
#include"gmputil.h"

using namespace std;

int main()
{
  uint128 base, e, n;
  cout << "Enter base: "; cin>> base;
  cout << "Enter exp: "; cin>> e;
  cout << "Enter modulus: "; cin>> n;

  uint64_t z[2];
  MontgomeryBig::topair(z,n);
  n=MontgomeryBig::frompair(z);
  cout << " pairs; n="<<n<<endl;


  cout << "Mont=" << MontgomeryBig::powmod(base,e,n) << endl;

  mpz_t B,E,M,ans;
  mpz_init(B); bigint2mpz(B,base);
  mpz_init(E); bigint2mpz(E,e);
  mpz_init(M); bigint2mpz(M,n);
  mpz_init(ans);
  mpz_powm(ans,B,E,M);
  cout << "GMP=" << ans << endl;


  return 0;
}
