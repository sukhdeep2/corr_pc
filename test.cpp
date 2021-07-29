#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <fstream>
#include <string>
#include <sstream>
#include <gsl/gsl_math.h>
using namespace std;

int f1(double a)
{
  return a;
}

int f2(int a,int (*f)(double))
{
  return f((double)a)*2;
}

double drand()
{
  srand(123123123);
  return rand() / (RAND_MAX + 1.);
}

/*int main ()
{
  int (*f)(double);
  f=&f1;
  int a=10;
  int b=f2(a,f1);
  string ad="aaaa";
  string d="bbbb";
  string c="aaaa";
  if (ad.compare(c)==0) cout<< ad;
  //  if (a==d) cout<< d;
  //  if(d==c) cout<< d;
  cout<<b<<endl;

  string file="coor.cpp";
  string comm="wc -l ";
  comm.append(file);
  
  long al=system(comm.c_str());
  cout<<al<<endl;
  return 0;
  };*/

int main()
{
  int j=0;
  cout<<RAND_MAX<<endl;
  while (j<10)
    {
      cout<<drand()<<endl;
      j+=1;
    }
}
