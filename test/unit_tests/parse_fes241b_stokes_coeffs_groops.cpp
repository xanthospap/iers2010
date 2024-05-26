#include "ocean_tide.hpp"
#include <array>

constexpr const int DEGREE = 180;
constexpr const int ORDER = 180;

using namespace dso;

int main(int argc, char *argv[]) {
  if (argc != 2) {
    fprintf(stderr,
            "Usage: %s [FES2014 GROOPS STOKES COEFFS DIR] "
            "\n",
            argv[0]);
    return 1;
  }
  
  /* create an OceanTide instance using the fes2014b gfc files */
  auto fes14b =
      dso::initFes2014bFromIcgem(argv[1], "fes2014b_n180_version20170520");

  auto w = fes14b.find_wave(DoodsonConstituent({2,2,-3,0,0,1}));
  /* 272.556 (t2) sin */
  assert(w->sin_coeffs().C(3, 1) == 8.715122060108e-12);
  assert(w->sin_coeffs().C(99, 77) == -3.921804503917e-15);
  assert(w->sin_coeffs().C(118, 34) == 3.118436203255e-15);
  assert(w->sin_coeffs().C(136, 120) == -3.566973052004e-15);
  assert(w->sin_coeffs().C(154, 30) == -7.594205054205e-16);
  assert(w->sin_coeffs().C(169, 0) == -1.862898266095e-15);
  assert(w->sin_coeffs().C(176, 124) == 9.650263945381e-16);
  assert(w->sin_coeffs().C(180, 180) == -3.951604599436e-17);
  assert(w->sin_coeffs().S(180, 180) == -6.467020899650e-16);
  assert(w->sin_coeffs().S(180, 112) == -1.297283014988e-15);
  assert(w->sin_coeffs().S(180, 0) == 0.000000000000e+00);
  assert(w->sin_coeffs().S(115, 78) == -5.941469186423e-15);
  assert(w->sin_coeffs().S(83, 45) == -9.946589886221e-15);
  assert(w->sin_coeffs().S(65, 22) == -7.411832484347e-15);
  assert(w->sin_coeffs().S(36, 17) == -9.632379117539e-14);
  assert(w->sin_coeffs().S(3, 1) == 2.120773048563e-12);

  /* 272.556 (t2) cos */
  assert(w->cos_coeffs().C(1, 0) == 1.871758347879e-12);
  assert(w->cos_coeffs().C(4, 3) == 2.915306134701e-13);
  assert(w->cos_coeffs().C(58, 17) == -3.900490900286e-15);
  assert(w->cos_coeffs().C(98, 23) == 3.678628530507e-15);
  assert(w->cos_coeffs().C(133, 24) == 2.527394775140e-15);
  assert(w->cos_coeffs().C(161, 85) == -2.150436730519e-15);
  assert(w->cos_coeffs().C(162, 81) == 1.302551598344e-15);
  assert(w->cos_coeffs().C(168, 91) == -3.884585962159e-16);
  assert(w->cos_coeffs().S(175, 106) == -1.354076194331e-16);
  assert(w->cos_coeffs().S(180, 179) == -4.205281211999e-16);
  assert(w->cos_coeffs().S(180, 180) == 4.992869238187e-16);
  assert(w->cos_coeffs().S(122, 84) == 7.120204985437e-16);
  assert(w->cos_coeffs().S(122, 83) == -1.547292127783e-15);
  assert(w->cos_coeffs().S(109, 34) == 3.451831921860e-15);
  assert(w->cos_coeffs().S(46, 24) == -6.694585279134e-14);
  assert(w->cos_coeffs().S(19, 19) == 1.738665271027e-14);
  
  /* 057.555 (ssa) sin */
  w = fes14b.find_wave(DoodsonConstituent({0,0,2,0,0,0}));
  assert(w->sin_coeffs().C(1  ,  0  ) ==   7.214255048141e-14 ); 
  assert(w->sin_coeffs().C(1  ,  1  ) ==  -1.879708896727e-12 ); 
  assert(w->sin_coeffs().C(43 ,  17 ) ==  -2.875346779425e-16 );
  assert(w->sin_coeffs().C(43 ,  43 ) ==  -2.805333544103e-16 );
  assert(w->sin_coeffs().C(161,   79) ==  -5.668770010875e-17 );
  assert(w->sin_coeffs().C(161,  161) ==   1.789432249234e-16 );
  assert(w->sin_coeffs().C(180,  179) ==   4.042790223195e-17 );
  assert(w->sin_coeffs().C(180,  180) ==  -1.491075395915e-16 );
  assert(w->sin_coeffs().S(180,  180) ==   -5.957725359298e-17);
  assert(w->sin_coeffs().S(180,  179) ==    1.660812826773e-18);
  assert(w->sin_coeffs().S(180,  144) ==   -6.037516789349e-17);
  assert(w->sin_coeffs().S(82 ,  28 ) ==    1.753796204003e-17);
  assert(w->sin_coeffs().S(63 ,   0 ) ==    0.000000000000e+00);
  assert(w->sin_coeffs().S(62 ,  62 ) ==    1.257003557915e-15);
  assert(w->sin_coeffs().S(38 ,  38 ) ==   -3.205744262855e-15);
  assert(w->sin_coeffs().S(4  ,  3  ) ==   -1.458355775156e-13);

  /* 057.555 (ssa) cos */
  assert(w->cos_coeffs().C(0,    0) ==   0.000000000000e+00);
  assert(w->cos_coeffs().C(1,    0) ==   2.008319822233e-11); 
  assert(w->cos_coeffs().C(1,    1) ==  -1.671287595828e-11); 
  assert(w->cos_coeffs().C(2,    0) ==  -5.635396618747e-11); 
  assert(w->cos_coeffs().C(2,    1) ==  -1.151456060440e-12); 
  assert(w->cos_coeffs().C(2,    2) ==  -1.585229398111e-12); 
  assert(w->cos_coeffs().C(3,    0) ==  -6.255322306044e-13); 
  assert(w->cos_coeffs().C(3,    1) ==   4.392068567700e-12); 
  assert(w->cos_coeffs().S(3,    2) ==  -6.676908308010e-13);
  assert(w->cos_coeffs().S(3,    3) ==  -3.629935284785e-12);
  assert(w->cos_coeffs().S(4,    0) ==   0.000000000000e+00);
  assert(w->cos_coeffs().S(4,    1) ==   7.675250119918e-13);
  assert(w->cos_coeffs().S(4,    2) ==   3.221882652340e-13);
  assert(w->cos_coeffs().S(4,    3) ==   1.144177157346e-12);
  assert(w->cos_coeffs().S(4,    4) ==  -5.209013657490e-12);
  assert(w->cos_coeffs().S(5,    1) ==   2.759455708717e-12);
  
  /* 056.554 (sa) sin */
  w = fes14b.find_wave(DoodsonConstituent({0,0,1,0,0,-1}));
  for (int d=0; d<180; d++) {
    for (int c=0; c<=d; c++) {
      assert(w->sin_coeffs().S(c,d) == 0e0);
      assert(w->sin_coeffs().C(c,d) == 0e0);
    }
  }
  
  /* 056.554 (sa) cos */
  assert(w->cos_coeffs().C(180,  170) ==  -1.575152348962e-17); 
  assert(w->cos_coeffs().C(180,  171) ==  -6.095537636953e-17); 
  assert(w->cos_coeffs().C(180,  172) ==  -5.164873665599e-17); 
  assert(w->cos_coeffs().C(180,  173) ==   3.679700343749e-17); 
  assert(w->cos_coeffs().C(180,  174) ==   2.477579259748e-17); 
  assert(w->cos_coeffs().C(180,  175) ==  -5.106305346296e-17); 
  assert(w->cos_coeffs().C(180,  176) ==   2.408677426627e-17); 
  assert(w->cos_coeffs().C(180,  177) ==   3.514117503618e-17); 
  assert(w->cos_coeffs().S(180,  178) ==   1.110745915902e-16);
  assert(w->cos_coeffs().S(180,  179) ==  -4.375835721966e-17);
  assert(w->cos_coeffs().S(180,  180) ==   7.638358564079e-17);
  assert(w->cos_coeffs().S(180,  170) ==   1.017557501383e-16);
  assert(w->cos_coeffs().S(180,  171) ==   5.991122958607e-17);
  assert(w->cos_coeffs().S(180,  172) ==   5.509014486131e-17);
  assert(w->cos_coeffs().S(180,  173) ==  -6.649193664111e-17);
  assert(w->cos_coeffs().S(180,  174) ==  -4.091627432427e-17);

  /* 491.555 (s4) sin */
  w = fes14b.find_wave(DoodsonConstituent({4,4,-4,0,0,0}));
  assert(w->sin_coeffs().C( 1 ,   0 ) ==  -5.861960130196e-14); 
  assert(w->sin_coeffs().C( 0 ,   0 ) ==   0.000000000000e+00); 
  assert(w->sin_coeffs().C(123,   12) ==  -4.689970409562e-16); 
  assert(w->sin_coeffs().C(123,   36) ==   3.760804829115e-16); 
  assert(w->sin_coeffs().C(174,   52) ==  -1.524583727687e-16); 
  assert(w->sin_coeffs().C(174,   40) ==  -1.120707451998e-16); 
  assert(w->sin_coeffs().C(180,  179) ==  -2.669644076679e-17); 
  assert(w->sin_coeffs().C(180,  180) ==   3.733082106223e-16); 
  assert(w->sin_coeffs().S(180,  100) ==  -3.188302717552e-16);
  assert(w->sin_coeffs().S(180,   10) ==   4.151804496720e-16);
  assert(w->sin_coeffs().S(180,    0) ==   0.000000000000e+00);
  assert(w->sin_coeffs().S(180,    1) ==   2.969370901708e-16);
  assert(w->sin_coeffs().S(180,    2) ==  -1.563327860999e-16);
  assert(w->sin_coeffs().S(179,   38) ==   6.104870433984e-17);
  assert(w->sin_coeffs().S(174,    1) ==  -5.591085225555e-16);
  assert(w->sin_coeffs().S(173,  173) ==   1.692920999945e-16);

  assert(w->cos_coeffs().C(2   , 1  ) ==  -7.456795667932e-14); 
  assert(w->cos_coeffs().C(2   , 2  ) ==   1.188216538568e-13);  
  assert(w->cos_coeffs().C(3   , 0  ) ==  -1.047547078628e-13); 
  assert(w->cos_coeffs().C(154 ,  44) ==   3.480394634832e-17); 
  assert(w->cos_coeffs().C(154 ,  45) ==   3.592307506442e-16); 
  assert(w->cos_coeffs().C(154 ,  46) ==  -3.056791744095e-16); 
  assert(w->cos_coeffs().C(99  , 99 ) ==   7.078449921728e-16); 
  assert(w->cos_coeffs().C(100 ,   0) ==  -4.776111721987e-16); 
  assert(w->cos_coeffs().S(99  , 99 ) ==   5.722764917164e-16);
  assert(w->cos_coeffs().S(100 ,   1) ==   7.836057816335e-16);
  assert(w->cos_coeffs().S(160 , 159) ==   2.634530442708e-17);
  assert(w->cos_coeffs().S(160 , 160) ==  -2.137127658690e-17);
  assert(w->cos_coeffs().S(161 ,   1) ==   2.555597341328e-16);
  assert(w->cos_coeffs().S(161 ,   2) ==   6.322543009537e-16);
  assert(w->cos_coeffs().S(35  ,  1 ) ==   1.059983763834e-14);
  assert(w->cos_coeffs().S(35  ,  5 ) ==   6.633617360653e-15);

  return 0;
}
