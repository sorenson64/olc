#ifndef _RPRIMETEST
#define _RPRIMETEST

#include "bigint.h"
#include "gmputil.h"

mpz_t two,Q,N,answer,save;
mpz_t y; 

bool issquare(int128 n)
{
  bigint2mpz(y,n);
  bool answer=mpz_perfect_square_p(y);
  return answer;
}

int64_t countprimetests=0;

void primetestinits() {
  mpz_init(y);
  mpz_init_set_ui(two,2);
  mpz_init(Q); 
  mpz_init(N); 
  mpz_init(answer);
  mpz_init(save);
}
// Theorem 4.1.5 from Crandall and Pomerance
// Brillhart-Lehmer-Selfridge
// assumes R is prime and n^1/3 < R < n^1/2, R | n-1
bool primetest(int128 n, int64_t R)
{
  countprimetests++;
  int128 q=(n-1)/R;

  //cout << "In primetest: n="<<n<<endl;

  bigint2mpz(Q,q);
  bigint2mpz(N,n);

  mpz_mul_ui(N,Q,R);
  mpz_add_ui(N,N,1);

  //cout << "R="<<R<<" Q="<<Q<<" N="<<N<<endl;

  //check gcd(a^(n-1)/R -1, n)=1
  mpz_powm(save,two,Q,N);
  mpz_sub_ui(answer,save,1);
  mpz_gcd(answer,answer,N);
  int64_t ans=mpz_get_si(answer);
  //cout << "GCD check, ans="<<ans<<endl;
  if(ans!=1) return false;

  //check a^(n-1) mod n ==1
  mpz_powm_ui(answer,save,R,N);
  ans=mpz_get_si(answer);
  //cout << "Fermat check, ans="<<ans<<endl;
  if(ans!=1) return false;

  //write n = c2*R^2 + c1*R + a0
  //check c1^2-4*c2 not a square

  int64_t c1= q % R;
  int64_t c2= (q-c1)/R;
  //cout << "c1="<<c1<<" c2="<<c2<<endl;
  bool result=!issquare(((int128)c1)*c1 - 4*c2);
  //cout << "Square check, ans="<<result<<endl;
  return result;
}

#endif
