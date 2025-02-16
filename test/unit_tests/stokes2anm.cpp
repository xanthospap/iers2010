#include "stokes_coefficients.hpp"
#include "eigen3/Eigen/Eigen"
#ifdef NDEBUG
#undef NDEBUG
#endif
#include <cassert>
#include <iostream>
#include <random>

// @snippet example_usage

/** Tests to validate the StokesCoeffs.to_anm() function.
 * 
 * To get the element corresponding to C(l,m) or S(l,m) from the (resulting) 
 * anm matrix, use the 'getters' anmC and anmS.
 */

std::uniform_real_distribution<double> unif(-1e0, 1e0);
std::default_random_engine re;

using namespace dso;

double anmC(int row, int col, const Eigen::MatrixXd *anm) {
  assert(row >= col );
  return anm->operator()(row, col);
}
double anmS(int row, int col, const Eigen::MatrixXd *anm) {
  assert( row >= col );
  return (col==0)?0e0:anm->operator()(col-1, row);
}

int main() {

  // Test 1:
  {
  printf("Test Case 1\n");
    const int deg=10;
    const int ord=10;
  StokesCoeffs s(deg);

  for (int m = 0; m <= ord; m++) {
    for (int l = m; l <= deg; l++) {
        s.C(l, m) = 1e0;
        if (m!=0) s.S(l,m) = 2e0;
    }
}

  auto anm = s.to_anm();
  //std::cout<<"Anm\n"<<anm<<"\n";

  for (int row=0;row<=deg;row++) {
    for (int col=0;col<=row;col++) {
      assert( s.C(row,col) == anmC(row,col,&anm) );
      if(col!=0)
        assert( s.S(row,col) == anmS(row,col,&anm) );
    }
  }
  }
  
  // Test 2:
  {
  printf("Test Case 2\n");
    const int deg=20;
    const int ord = 20;
  StokesCoeffs s(deg);

  for (int m = 0; m <= ord; m++) {
    for (int l = m; l <= deg; l++) {
        s.C(l, m) = unif(re);
        if (m!=0) s.S(l, m) = unif(re);
    }
}

  auto anm = s.to_anm();
  //std::cout<<"\nAnm\n"<<anm;

  for (int row=0;row<=deg;row++) {
    for (int col=0;col<=row;col++) {
      assert( s.C(row,col) == anmC(row,col,&anm) );
      if(col!=0)
        assert( s.S(row,col) == anmS(row,col,&anm) );
    }
  }
  }
  
  // Test 3:
  {
  printf("Test Case 3\n");
    const int deg=26;
    const int ord =19;
  StokesCoeffs s(deg, ord);

  for (int m = 0; m <= ord; m++) {
    for (int l = m; l <= deg; l++) {
        s.C(l, m) = unif(re);
        if (m!=0) s.S(l, m) = unif(re);
    }
}

  auto anm = s.to_anm();

  for (int row=0;row<=deg;row++) {
    for (int col=0;col<=std::min(row,ord);col++) {
      assert( s.C(row,col) == anmC(row,col,&anm) );
      if(col!=0)
        assert( s.S(row,col) == anmS(row,col,&anm) );
    }
  }
  }
  
  // Test 4:
  {
  printf("Test Case 4\n");
    const int deg=26;
    const int ord =25;
  StokesCoeffs s(deg, ord);

  for (int m = 0; m <= ord; m++) {
    for (int l = m; l <= deg; l++) {
        s.C(l, m) = unif(re);
        if (m!=0) s.S(l, m) = unif(re);
    }
}

  auto anm = s.to_anm();

  for (int row=0;row<=deg;row++) {
    for (int col=0;col<=std::min(row,ord);col++) {
      assert( s.C(row,col) == anmC(row,col,&anm) );
      if(col!=0)
        assert( s.S(row,col) == anmS(row,col,&anm) );
    }
  }
  }
  
  // Test 5:
  {
  printf("Test Case 5\n");
    const int deg=126;
    const int ord =126;
  StokesCoeffs s(deg, ord);

  for (int m = 0; m <= ord; m++) {
    for (int l = m; l <= deg; l++) {
        s.C(l, m) = unif(re);
        if( m!=0) s.S(l, m) = unif(re);
    }
}
  auto anm = s.to_anm();

  for (int row=0;row<=deg;row++) {
    for (int col=0;col<=std::min(row,ord);col++) {
      assert( s.C(row,col) == anmC(row,col,&anm) );
      if(col!=0)
        assert( s.S(row,col) == anmS(row,col,&anm) );
    }
  }
  }

  return 0;
}

// @endsnippet
