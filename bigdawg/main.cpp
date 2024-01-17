#include<iostream>
#include<fstream>
#include<bitset>
#include<cmath>
#include<gmp.h>

using namespace std;

const int64_t minintervalsize=128;
const int64_t numintervals=256; // 256 for desktop
const int64_t maskwidth=minintervalsize*numintervals*4;

#include "int.h"
#include "smallprimes.h"
#include "rlist.h"
#include "masks.h"
#include "primetest.h"

// for data reporting

uint64_t countprimesfound;
//uint64_t countprimetests; // this is declared in primetest.h
uint64_t countfails;
uint64_t countR[5];
uint64_t countfilter, countfix, countfix2;
uint64_t countT;
long double totalT;
uint64_t countk, totalk;
uint64_t nbegin, nend;

void countreset()
{
  countprimesfound=0;
  countprimetests=0;
  countfails=0;
  countR[0]= countR[1]= countR[2]= countR[3]= countR[4]=0;
  countfilter=countfix=countfix2=0;
  countT=0;
  totalT=0;
  countk=totalk=0;
}

void countprint(ostream & os)
{
  os << nbegin <<" "<< nend 
     << " \t"<<countprimesfound<<" "<<countprimetests<<" "<<countfails
     << " \t"<<countR[0]<<" "<<countR[1]<<" "<<countR[2]
     <<" "<<countR[3]<<" "<<countR[4]
     << " \t"<<countfilter<<" "<<countfix<<" "<<countfix2
     << " \t"<<(((long double)totalT)/countT)
     << " \t"<<(((long double)totalk)/countk) << endl;
}

int64_t findR(int128 n, int col)
{
  long double bound=powl(n,1.0L/3.0L);
  int i=0;
  while(i<rlistlen && rlist[i++][col]<=bound) ;
  return rlist[i][col];
}

int64_t findR(int128 n) { return findR(n,0); }

int64_t findstart(int p, int64_t M, int128 qleft)
{
  int Mmodp=M%p;
  int qmodp=qleft%p;
  int start=p-(inv(Mmodp,p)+qmodp)%p;
  if(start==p) start=0;
  return start;
}

bool filter(int128 n)
{
  countfilter++;
  // trial division by small primes vi GCD
  const int64_t m=3*5*7*11*13*17;
  int64_t rem=n%m;
  return (gcd(m,rem)==1);
}

mpz_t N1,N2;

bool fixit2(int128 nstart, int128 nstop)
{
  countfix2++;
  bigint2mpz(N1,nstart);
  bigint2mpz(N2,nstop);

  while(mpz_cmp(N1,N2)<0) // while N1<N2
  {
    mpz_nextprime(N1,N1);
    if(mpz_cmp(N1,N2)>0)
    {
      return false;
    }
    int val=mpz_probab_prime_p(N1,64);
    if(val==2) 
    {
      return true;
    }
    if(val==1) // need to prove primality later
    {
      ofstream f("primes2check.txt",ios::app);
      f << N1 << endl;
      f.close();
      return true;
    }
  }
  return false;
}

bool fixit(int128 nleft, int128 nright)
{
  countfix++;
  for(int r=1; r<5; r++)
  {
    int64_t R=findR(nright,r);
    int128 qstart=nleft/(2*R)+1;
    int128 qstop=nright/(2*R);
    int64_t qdelta=qstop-qstart;
    int128 n=qstart*2*R+1;
    for(int i=0; i<qdelta; i++)
    {
      if(filter(n))
      {
        if(primetest(n,R)) { countR[r]++; return true; }
      }
      n=n+2*R;
    }
  }
  return fixit2(nleft,nright);
}

int main()
{
  makemasks();
  primetestinits();
  countreset();
  mpz_init(N1); 
  mpz_init(N2); 
  Smallprimes p(maskwidth);
  bitset<maskwidth> x; // bit "vector" for sieving
  //printmasks();
  int64_t n=(int64_t)1e10; // start point
  nbegin=n;

  int64_t loopcount=0;
  bitset<maskwidth> y;

  //while(true) // BIG LOOP on n
  //while(n<2e9+1e7) // BIG LOOP on n
  while(n<1e10+1e7) // BIG LOOP on n
{
  loopcount++;

  //const int t=numintervals;
  int64_t T=numintervals;

  //int64_t R=findR(((int128)(n+t))*(n+t)); // need R>(n+t)^{2/3}
  int64_t R=findR(((int128)(n+T+T))*(n+T+T)); // need R>(n+2t)^{2/3}

  // compute the modulus M=R*2*small primes
  int64_t M=2*R;
  int k=1;

  while(n/(M*primes[k])>minintervalsize) { M *= primes[k++]; }

  while(n/(2*M)>minintervalsize) { M *= 2; }

  totalk+=k; countk++;

  // compute T
  int64_t intervalsize=n/M;

  T=(maskwidth/2)/(intervalsize+1);

  totalT+=T; countT++;

  //T=numintervals;

  //cout << "  T="<<T<<" numintervals="<<numintervals<<endl;

  // numbers are of the form M*q+1, from q=left to right
  // so qleft is bit position 0 of the bitset x
  int128 nsquared=((int128)n)*((int128)n);
  int128 nplustsquared=((int128)n+T)*((int128)n+T);
  int128 qleft=(nsquared/M)+1; // n squared is not prime
  int128 qright=(nplustsquared/M);
  if((n+T)%M==0) qright--;

  /* * /
  cout << "n=" << n << " R="<<R<<" k="<<k<<" M="<<M<<" intv size="<< n/M
	  <<endl;
  cout << "target n="<<n+T
	  << " each interval averages "<< (2.0*n*T-T*T)/(2*T*M) 
	  << " total size=" << (2.0*n*T-T*T)/M 
	  << endl;
  cout << "qleft=" << qleft
	  << " qright=" << qright 
	  << " qdelta=" << qright-qleft+1
	  << endl;
  cout << "maskwidth="<<maskwidth<<endl;
  / * */
  int64_t qdelta=qright-qleft+1;

  x.reset(); // all bits are zero
  // do masks first
  for(int i=1; i<masklen; i++) // not using k here on purpose
    if(M%primes[i]!=0) // don't do primes in the modulus!
      {
	int p=primes[i];
        int start=findstart(p,M,qleft);
	// check that (M*(qleft+start)+1 is divisible by p
	//if( ((M%p)*(qleft+start)%p +1 )%p != 0) 
	//{ cout << "Error! "<< p << endl; break; }
	y=mask[i];
	y <<= start;
	x |= y;
	/*
	cout << primes[i] << endl;
	cout << "---" << endl;
	cout << mask[i] << endl;
	cout << "---" << endl;
	cout << (mask[i]<<start) << endl;
	cout << "---" << endl;
	cout << x << endl;
	cout << endl;
	*/
      }
  for(int i=masklen; i<p.length(); i++)
  {
    int start= findstart(p[i],M,qleft);
    //if( ((M%p[i])*(qleft+start)%p[i] +1 )%p[i] != 0) 
	//{ cout << "Error! "<< p[i] << endl; break; }
    for(int j=start; j<qdelta; j+=p[i]) x.set(j);
  }
  //cout << endl << "---" << endl;
  //cout << x << endl;


  for(int t=0; t<T; t++)
  {
    // left half (A)  [(n+t)^2,(n+t)(n+t+1)]
    int128 nstart=((int128)(n+t))*(n+t);
    int128 qstart=nstart/M+1;
    int128 nstop=((int128)(n+t))*(n+t+1);
    int128 qstop=nstop/M;

    bool done=false;
    for(int j=qstart-qleft;!done && j<qstop-qleft; j++)
      if(!x.test(j))
      {
        int128 n=((int128)M)*(qleft+j)+1;
        if(primetest( n, R )) done=true;
      }

    if(done) countR[0]++;
    if(!done) countfails++;
    if(!done) done=fixit(nstart,nstop);

    if(done) countprimesfound++;

    // right half (B)  [(n+t)(n+t+1),(n+t+1)^2]
    nstart=((int128)(n+t))*(n+t+1);
    qstart=nstart/M+1;
    nstop=((int128)(n+t+1))*(n+t+1);
    qstop=nstop/M;

    done=false;
    for(int j=qstart-qleft;!done && j<qstop-qleft; j++)
      if(!x.test(j))
      {
        int128 n=((int128)M)*(qleft+j)+1;
        if(primetest( n, R )) done=true;
      }

    if(done) countR[0]++;
    if(!done) countfails++;
    if(!done) done=fixit(nstart,nstop);

    if(done) countprimesfound++;
  }

  /*
  cout << "Total Failure count="<<failcount<<endl;
  cout << "Total Primes found="<<primesfound<<endl;
  cout << "Number of prime tests performed="<<numtests
	  << " success rate="<< (100.0*primesfound)/numtests<<endl;
  cout << "Processed from "<< n <<"^2 to "<< n+numintervals <<"^2"<<endl;
  */

  n+=T;

  if(loopcount%10000==0)
  {
  cout << "n="<<n<<endl;
  cout << "Failure count="<<countfails<<endl;
  cout << "Primes found="<<countprimesfound<<endl;
  cout << "Number of prime tests performed="<<countprimetests
	  << " success rate="<< (100.0*countprimesfound)/countprimetests<<endl;
  }
} // BIG loop on n
  nend=n;
  
/*
  cout << "\nFinal n="<<n<<endl;
  cout << "Failure count="<<failcount<<endl;
  cout << "Primes found="<<primesfound<<endl;
  cout << "Number of prime tests performed="<<numtests
	  << " success rate="<< (100.0*primesfound)/numtests<<endl;
*/
  countprint(cout);

  return 0;
}
