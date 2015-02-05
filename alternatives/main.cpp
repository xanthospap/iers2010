#include "shclass.hpp"
#include <ctime>

int main ()
{
  double d, grn, gre;
  std::clock_t start;
  double dmjd;
  
  printf ("\nOne call to each function costs:");
  dmjd = 55055e0;
  
  start = std::clock();
  iers2010::gpt (dmjd,0.6708665767e0,-1.393397187e0,812.546e0,d,grn,gre);
  printf ("\n\tTime for gpt v1 was %15.5f ms.", (std::clock() - start) / (double)(CLOCKS_PER_SEC / 1000));
  
  start = std::clock();
  iers2010::gpt2 (dmjd,0.6708665767e0,-1.393397187e0,812.546e0,d,grn,gre);
  printf ("\n\tTime for gpt v2 was %15.5f ms.", (std::clock() - start) / (double)(CLOCKS_PER_SEC / 1000));
  
  start = std::clock();
  iers2010::gpt3 (dmjd,0.6708665767e0,-1.393397187e0,812.546e0,d,grn,gre);
  printf ("\n\tTime for gpt v3 was %15.5f ms.", (std::clock() - start) / (double)(CLOCKS_PER_SEC / 1000));
  
  int calls = 10000;
  printf ("\nMultiple (%5i) to each function cost:", calls);
  dmjd = 55055e0;
  start = std::clock();
  for (int i = 0;i<calls;i++) {
    iers2010::gpt (dmjd,0.6708665767e0,-1.393397187e0,812.546e0,d,grn,gre);
    dmjd += (i/1e03);
  }
  printf ("\n\tTime for gpt v1 was %15.5f ms.", (std::clock() - start) / (double)(CLOCKS_PER_SEC / 1000));
  
  dmjd = 55055e0;
  start = std::clock();
  for (int i = 0;i<calls;i++) {
    iers2010::gpt2 (dmjd,0.6708665767e0,-1.393397187e0,812.546e0,d,grn,gre);
    dmjd += (i/1e03);
   }
  printf ("\n\tTime for gpt v2 was %15.5f ms.", (std::clock() - start) / (double)(CLOCKS_PER_SEC / 1000));
  
  dmjd = 55055e0;
  start = std::clock();
  for (int i = 0;i<calls;i++) {
      iers2010::gpt3 (dmjd,0.6708665767e0,-1.393397187e0,812.546e0,d,grn,gre);
      dmjd += (i/1e03);
  }
  printf ("\n\tTime for gpt v3 was %15.5f ms.", (std::clock() - start) / (double)(CLOCKS_PER_SEC / 1000));

  printf ("\n\n");
  return 0;
}
