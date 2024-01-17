#ifndef _smallprimes
#define _smallprimes

#include<vector>
#include<cstdint>
#include<iostream>
#include<iomanip>

using namespace std;

class Smallprimes
{

  private:
    vector<int64_t> p;
    uint64_t pos;

  public:
    Smallprimes():pos(0) {}
    Smallprimes(int64_t x):pos(0) { find(x); }

    inline int64_t operator[](uint64_t i) const { return p[i]; }
    inline int64_t max() const                  { return p.back(); }
    inline int64_t length() const               { return p.size(); }
    inline void reset(uint64_t start=0)         { pos=start; }
    inline int64_t next()                       { return p[pos++]; }
    inline bool hasNext() const                 { return pos<p.size(); }
    inline int operator!() const                { return p.size()==0; }

  //private:
    void find(int64_t x) // find primes up to x
    {
      if(x<2) return;
      if(p.size()>0 && max()>=x) return;
      pos=0;
      vector<bool> s;
      s.resize(x+1,true);

      s[0]=s[1]=false;
      for(int64_t i=4; i<=x; i+=2) s[i]=false;
      for(int64_t r=3; r*r<=x; r+=2)
        for(int64_t q=r*r; q<=x; q+=2*r) s[q]=false;

      int64_t count=0;
      for(int64_t i=2; i<=x; i++)
        if(s[i]) count++;

      p.clear();
      p.reserve(count);
      for(int64_t i=2; i<=x; i++)
        if(s[i]) p.push_back(i);
    }

  public:
    static void testme()
    {
       int primesupto100[]={
         2, 3, 5, 7, 11, 
         13, 17, 19, 23, 29, 
         31, 37, 41, 43, 47, 
         53, 59, 61, 67, 71, 
         73, 79, 83, 89, 97};
       int length=25;
       Smallprimes p(100);
       if(p.length()!=length)
         { cout << "Fail!\n"; return; }
       if(p.max()!=97)
         { cout << "Fail!!\n"; return; }
       for(int i=0; i<length; i++)
         if(p[i]!=primesupto100[i]) 
           { cout << "Fail!!!\n"; return; }
       cout << "Pass.\n";
    }

};

ostream & operator<<(ostream & os, const Smallprimes & P)
{
  int64_t i,p;
  int width=1;
  for(p=P.max(); p>0; p/=10) width++;
  for(i=0; i<P.length(); i++)
    { os << setw(width) << P[i] ; if(i%10==9) os << endl; }
  os << endl;
  return os;
}

#endif
