#ifndef _MASKS
#define _MASKS

#include<iostream>
#include "int.h"

using namespace std;

//const int maskwidth=64; // this is defined in the main
bitset<maskwidth> mask[15];
const int masklen=15;

void makemasks()
{
  for(int i=1; i<masklen; i++)
  {
    int p=primes[i];
    mask[i].reset();
    for(int q=0; q<maskwidth; q+=p) mask[i].set(q);
  }
}

void printmasks()
{
  for(int i=1; i<masklen; i++)
  {
    cout << primes[i] <<": "<< mask[i] << endl;
  }
}

#endif
