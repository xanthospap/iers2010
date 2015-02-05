#include <stdio.h>
#include <cmath>

namespace iers2010 {

int gpt (const double& dmjd,const double& dlat,const double& dlon,
  const double& dhgt,double& pres,double& temp,double& undu);
int gpt2 (const double& dmjd,const double& dlat,const double& dlon,
  const double& dhgt,double& pres,double& temp,double& undu);

class SphericalHarmonics {
public:
  SphericalHarmonics ();
  SphericalHarmonics (const double& x, const double& y, const double& z);
  ~SphericalHarmonics ();
  SphericalHarmonics (const SphericalHarmonics&);
  SphericalHarmonics& operator = (const SphericalHarmonics&);
  
  inline double x () const {return mx;}
  inline double y () const {return my;}
  inline double z () const {return mz;}
  inline double v (const int& i, const int& j) const {return mv[i][j];}
  inline double w (const int& i, const int& j) const {return mw[i][j];}
  
  inline void reset (const double& x, const double& y, const double& z)
  {
    mx = x;
    my = y;
    mz = z;
    this->populate ();
    return;
  }
  
private:
  void populate ();
  
  int mdegree, morder;
  double mx,my, mz;
  double** mv;
  double** mw;
};

}; //  end of namespace