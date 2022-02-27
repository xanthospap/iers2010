#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest.h"
#include "iers2010.hpp"

using iers2010::pmsdnut2;
#define PRECISION 1e-6

TEST_CASE(
    "PMSDNUT2: Comparisson of results based on Fortran implementation.") {
      double mjd, dx, dxr, dy, dyr;
      SUBCASE("Example   1: Checking for descripancies  > PRECISION") {
        mjd = 38777.273437500000000e0;
        dxr = -13.332430339347264e0;
        dyr = 36.989813911749813e0;
        pmsdnut2(mjd, dx, dy);
        REQUIRE(std::abs(dx - dxr) < PRECISION);
        REQUIRE(std::abs(dy - dyr) < PRECISION);
      }
      SUBCASE("Example   2: Checking for descripancies  > PRECISION") {
        mjd = 43280.726562500000000e0;
        dxr = -12.692326579475079e0;
        dyr = 18.309054292704225e0;
        pmsdnut2(mjd, dx, dy);
        REQUIRE(std::abs(dx - dxr) < PRECISION);
        REQUIRE(std::abs(dy - dyr) < PRECISION);
      }
      SUBCASE("Example   3: Checking for descripancies  > PRECISION") {
        mjd = 64722.160156250000000e0;
        dxr = -18.585504882515195e0;
        dyr = 4.816093506779952e0;
        pmsdnut2(mjd, dx, dy);
        REQUIRE(std::abs(dx - dxr) < PRECISION);
        REQUIRE(std::abs(dy - dyr) < PRECISION);
      }
      SUBCASE("Example   4: Checking for descripancies  > PRECISION") {
        mjd = 54519.675781250000000e0;
        dxr = -0.003947469537895e0;
        dyr = -0.681882443563385e0;
        pmsdnut2(mjd, dx, dy);
        REQUIRE(std::abs(dx - dxr) < PRECISION);
        REQUIRE(std::abs(dy - dyr) < PRECISION);
      }
      SUBCASE("Example   5: Checking for descripancies  > PRECISION") {
        mjd = 57066.449218750000000e0;
        dxr = 22.620261956612371e0;
        dyr = -3.806792963710297e0;
        pmsdnut2(mjd, dx, dy);
        REQUIRE(std::abs(dx - dxr) < PRECISION);
        REQUIRE(std::abs(dy - dyr) < PRECISION);
      }
      SUBCASE("Example   6: Checking for descripancies  > PRECISION") {
        mjd = 46284.660156250000000e0;
        dxr = 0.798039003477055e0;
        dyr = 10.352197096098882e0;
        pmsdnut2(mjd, dx, dy);
        REQUIRE(std::abs(dx - dxr) < PRECISION);
        REQUIRE(std::abs(dy - dyr) < PRECISION);
      }
      SUBCASE("Example   7: Checking for descripancies  > PRECISION") {
        mjd = 40377.191406250000000e0;
        dxr = 27.259339653100419e0;
        dyr = 4.949551339799813e0;
        pmsdnut2(mjd, dx, dy);
        REQUIRE(std::abs(dx - dxr) < PRECISION);
        REQUIRE(std::abs(dy - dyr) < PRECISION);
      }
      SUBCASE("Example   8: Checking for descripancies  > PRECISION") {
        mjd = 62085.742187500000000e0;
        dxr = 9.893871693187558e0;
        dyr = 1.613518661360981e0;
        pmsdnut2(mjd, dx, dy);
        REQUIRE(std::abs(dx - dxr) < PRECISION);
        REQUIRE(std::abs(dy - dyr) < PRECISION);
      }
      SUBCASE("Example   9: Checking for descripancies  > PRECISION") {
        mjd = 62100.199218750000000e0;
        dxr = -6.375195665589236e0;
        dyr = 4.915673115121399e0;
        pmsdnut2(mjd, dx, dy);
        REQUIRE(std::abs(dx - dxr) < PRECISION);
        REQUIRE(std::abs(dy - dyr) < PRECISION);
      }
      SUBCASE("Example  10: Checking for descripancies  > PRECISION") {
        mjd = 70875.992187500000000e0;
        dxr = -28.493615221399235e0;
        dyr = -8.975756212385377e0;
        pmsdnut2(mjd, dx, dy);
        REQUIRE(std::abs(dx - dxr) < PRECISION);
        REQUIRE(std::abs(dy - dyr) < PRECISION);
      }
      SUBCASE("Example  11: Checking for descripancies  > PRECISION") {
        mjd = 51937.242187500000000e0;
        dxr = -13.233523720191146e0;
        dyr = 8.565011499426708e0;
        pmsdnut2(mjd, dx, dy);
        REQUIRE(std::abs(dx - dxr) < PRECISION);
        REQUIRE(std::abs(dy - dyr) < PRECISION);
      }
      SUBCASE("Example  12: Checking for descripancies  > PRECISION") {
        mjd = 56607.449218750000000e0;
        dxr = -1.927087563038225e0;
        dyr = 8.368961537325319e0;
        pmsdnut2(mjd, dx, dy);
        REQUIRE(std::abs(dx - dxr) < PRECISION);
        REQUIRE(std::abs(dy - dyr) < PRECISION);
      }
      SUBCASE("Example  13: Checking for descripancies  > PRECISION") {
        mjd = 67311.914062500000000e0;
        dxr = -22.222715400927886e0;
        dyr = 9.062172233022295e0;
        pmsdnut2(mjd, dx, dy);
        REQUIRE(std::abs(dx - dxr) < PRECISION);
        REQUIRE(std::abs(dy - dyr) < PRECISION);
      }
      SUBCASE("Example  14: Checking for descripancies  > PRECISION") {
        mjd = 39948.613281250000000e0;
        dxr = -17.930588555995982e0;
        dyr = -12.881553329378832e0;
        pmsdnut2(mjd, dx, dy);
        REQUIRE(std::abs(dx - dxr) < PRECISION);
        REQUIRE(std::abs(dy - dyr) < PRECISION);
      }
      SUBCASE("Example  15: Checking for descripancies  > PRECISION") {
        mjd = 40597.175781250000000e0;
        dxr = -20.502031422883949e0;
        dyr = -0.274067640506393e0;
        pmsdnut2(mjd, dx, dy);
        REQUIRE(std::abs(dx - dxr) < PRECISION);
        REQUIRE(std::abs(dy - dyr) < PRECISION);
      }
      SUBCASE("Example  16: Checking for descripancies  > PRECISION") {
        mjd = 56960.722656250000000e0;
        dxr = 19.428101079715891e0;
        dyr = 5.656061494783518e0;
        pmsdnut2(mjd, dx, dy);
        REQUIRE(std::abs(dx - dxr) < PRECISION);
        REQUIRE(std::abs(dy - dyr) < PRECISION);
      }
      SUBCASE("Example  17: Checking for descripancies  > PRECISION") {
        mjd = 61820.636718750000000e0;
        dxr = 9.203449873235408e0;
        dyr = -26.513664689878667e0;
        pmsdnut2(mjd, dx, dy);
        REQUIRE(std::abs(dx - dxr) < PRECISION);
        REQUIRE(std::abs(dy - dyr) < PRECISION);
      }
      SUBCASE("Example  18: Checking for descripancies  > PRECISION") {
        mjd = 39025.503906250000000e0;
        dxr = -21.924010865646462e0;
        dyr = 1.330976620638623e0;
        pmsdnut2(mjd, dx, dy);
        REQUIRE(std::abs(dx - dxr) < PRECISION);
        REQUIRE(std::abs(dy - dyr) < PRECISION);
      }
      SUBCASE("Example  19: Checking for descripancies  > PRECISION") {
        mjd = 51934.675781250000000e0;
        dxr = 22.146911454655982e0;
        dyr = -9.787070409489433e0;
        pmsdnut2(mjd, dx, dy);
        REQUIRE(std::abs(dx - dxr) < PRECISION);
        REQUIRE(std::abs(dy - dyr) < PRECISION);
      }
      SUBCASE("Example  20: Checking for descripancies  > PRECISION") {
        mjd = 41057.714843750000000e0;
        dxr = -29.400915120438608e0;
        dyr = -8.079515589587679e0;
        pmsdnut2(mjd, dx, dy);
        REQUIRE(std::abs(dx - dxr) < PRECISION);
        REQUIRE(std::abs(dy - dyr) < PRECISION);
      }
      SUBCASE("Example  21: Checking for descripancies  > PRECISION") {
        mjd = 53105.757812500000000e0;
        dxr = -33.212295016724866e0;
        dyr = -1.652680670211601e0;
        pmsdnut2(mjd, dx, dy);
        REQUIRE(std::abs(dx - dxr) < PRECISION);
        REQUIRE(std::abs(dy - dyr) < PRECISION);
      }
      SUBCASE("Example  22: Checking for descripancies  > PRECISION") {
        mjd = 62357.019531250000000e0;
        dxr = 24.941949664846614e0;
        dyr = 7.471561625305756e0;
        pmsdnut2(mjd, dx, dy);
        REQUIRE(std::abs(dx - dxr) < PRECISION);
        REQUIRE(std::abs(dy - dyr) < PRECISION);
      }
      SUBCASE("Example  23: Checking for descripancies  > PRECISION") {
        mjd = 58997.296875000000000e0;
        dxr = 24.555824746206149e0;
        dyr = -12.053532958043480e0;
        pmsdnut2(mjd, dx, dy);
        REQUIRE(std::abs(dx - dxr) < PRECISION);
        REQUIRE(std::abs(dy - dyr) < PRECISION);
      }
      SUBCASE("Example  24: Checking for descripancies  > PRECISION") {
        mjd = 70729.179687500000000e0;
        dxr = 19.982583649351966e0;
        dyr = -17.119095127846730e0;
        pmsdnut2(mjd, dx, dy);
        REQUIRE(std::abs(dx - dxr) < PRECISION);
        REQUIRE(std::abs(dy - dyr) < PRECISION);
      }
      SUBCASE("Example  25: Checking for descripancies  > PRECISION") {
        mjd = 67834.671875000000000e0;
        dxr = -20.028288116131737e0;
        dyr = 18.187162378371710e0;
        pmsdnut2(mjd, dx, dy);
        REQUIRE(std::abs(dx - dxr) < PRECISION);
        REQUIRE(std::abs(dy - dyr) < PRECISION);
      }
      SUBCASE("Example  26: Checking for descripancies  > PRECISION") {
        mjd = 56865.863281250000000e0;
        dxr = 1.378763300540413e0;
        dyr = 18.752518362701863e0;
        pmsdnut2(mjd, dx, dy);
        REQUIRE(std::abs(dx - dxr) < PRECISION);
        REQUIRE(std::abs(dy - dyr) < PRECISION);
      }
      SUBCASE("Example  27: Checking for descripancies  > PRECISION") {
        mjd = 41920.781250000000000e0;
        dxr = -3.452444941864645e0;
        dyr = 20.869543370447850e0;
        pmsdnut2(mjd, dx, dy);
        REQUIRE(std::abs(dx - dxr) < PRECISION);
        REQUIRE(std::abs(dy - dyr) < PRECISION);
      }
      SUBCASE("Example  28: Checking for descripancies  > PRECISION") {
        mjd = 61228.531250000000000e0;
        dxr = -9.980183315098490e0;
        dyr = 6.463256120233743e0;
        pmsdnut2(mjd, dx, dy);
        REQUIRE(std::abs(dx - dxr) < PRECISION);
        REQUIRE(std::abs(dy - dyr) < PRECISION);
      }
      SUBCASE("Example  29: Checking for descripancies  > PRECISION") {
        mjd = 53054.750000000000000e0;
        dxr = 13.029101871781496e0;
        dyr = -23.351792525395410e0;
        pmsdnut2(mjd, dx, dy);
        REQUIRE(std::abs(dx - dxr) < PRECISION);
        REQUIRE(std::abs(dy - dyr) < PRECISION);
      }
      SUBCASE("Example  30: Checking for descripancies  > PRECISION") {
        mjd = 62852.558593750000000e0;
        dxr = 8.633977719149692e0;
        dyr = 0.704074975289882e0;
        pmsdnut2(mjd, dx, dy);
        REQUIRE(std::abs(dx - dxr) < PRECISION);
        REQUIRE(std::abs(dy - dyr) < PRECISION);
      }
      SUBCASE("Example  31: Checking for descripancies  > PRECISION") {
        mjd = 70038.554687500000000e0;
        dxr = -2.588141557958082e0;
        dyr = -0.061471355287246e0;
        pmsdnut2(mjd, dx, dy);
        REQUIRE(std::abs(dx - dxr) < PRECISION);
        REQUIRE(std::abs(dy - dyr) < PRECISION);
      }
      SUBCASE("Example  32: Checking for descripancies  > PRECISION") {
        mjd = 64949.582031250000000e0;
        dxr = -1.114786541938696e0;
        dyr = 23.236755259783717e0;
        pmsdnut2(mjd, dx, dy);
        REQUIRE(std::abs(dx - dxr) < PRECISION);
        REQUIRE(std::abs(dy - dyr) < PRECISION);
      }
      SUBCASE("Example  33: Checking for descripancies  > PRECISION") {
        mjd = 47778.378906250000000e0;
        dxr = -27.649611257485336e0;
        dyr = 1.039298003896462e0;
        pmsdnut2(mjd, dx, dy);
        REQUIRE(std::abs(dx - dxr) < PRECISION);
        REQUIRE(std::abs(dy - dyr) < PRECISION);
      }
      SUBCASE("Example  34: Checking for descripancies  > PRECISION") {
        mjd = 40391.851562500000000e0;
        dxr = -10.978274111322630e0;
        dyr = 21.078348503385211e0;
        pmsdnut2(mjd, dx, dy);
        REQUIRE(std::abs(dx - dxr) < PRECISION);
        REQUIRE(std::abs(dy - dyr) < PRECISION);
      }
      SUBCASE("Example  35: Checking for descripancies  > PRECISION") {
        mjd = 64051.136718750000000e0;
        dxr = -4.376702790334100e0;
        dyr = 11.066864941149387e0;
        pmsdnut2(mjd, dx, dy);
        REQUIRE(std::abs(dx - dxr) < PRECISION);
        REQUIRE(std::abs(dy - dyr) < PRECISION);
      }
      SUBCASE("Example  36: Checking for descripancies  > PRECISION") {
        mjd = 50038.210937500000000e0;
        dxr = -5.649118191023394e0;
        dyr = 2.650996574451751e0;
        pmsdnut2(mjd, dx, dy);
        REQUIRE(std::abs(dx - dxr) < PRECISION);
        REQUIRE(std::abs(dy - dyr) < PRECISION);
      }
      SUBCASE("Example  37: Checking for descripancies  > PRECISION") {
        mjd = 60497.414062500000000e0;
        dxr = -15.000866021807997e0;
        dyr = -33.417044158735571e0;
        pmsdnut2(mjd, dx, dy);
        REQUIRE(std::abs(dx - dxr) < PRECISION);
        REQUIRE(std::abs(dy - dyr) < PRECISION);
      }
      SUBCASE("Example  38: Checking for descripancies  > PRECISION") {
        mjd = 64750.121093750000000e0;
        dxr = -13.299554297600395e0;
        dyr = 2.134642497935274e0;
        pmsdnut2(mjd, dx, dy);
        REQUIRE(std::abs(dx - dxr) < PRECISION);
        REQUIRE(std::abs(dy - dyr) < PRECISION);
      }
      SUBCASE("Example  39: Checking for descripancies  > PRECISION") {
        mjd = 72811.484375000000000e0;
        dxr = -0.154836852510011e0;
        dyr = 3.845261296604908e0;
        pmsdnut2(mjd, dx, dy);
        REQUIRE(std::abs(dx - dxr) < PRECISION);
        REQUIRE(std::abs(dy - dyr) < PRECISION);
      }
      SUBCASE("Example  40: Checking for descripancies  > PRECISION") {
        mjd = 51313.492187500000000e0;
        dxr = -23.934058192940046e0;
        dyr = -13.330860307587482e0;
        pmsdnut2(mjd, dx, dy);
        REQUIRE(std::abs(dx - dxr) < PRECISION);
        REQUIRE(std::abs(dy - dyr) < PRECISION);
      }
      SUBCASE("Example  41: Checking for descripancies  > PRECISION") {
        mjd = 47248.519531250000000e0;
        dxr = 13.824856414906312e0;
        dyr = -7.453931753219721e0;
        pmsdnut2(mjd, dx, dy);
        REQUIRE(std::abs(dx - dxr) < PRECISION);
        REQUIRE(std::abs(dy - dyr) < PRECISION);
      }
      SUBCASE("Example  42: Checking for descripancies  > PRECISION") {
        mjd = 72520.046875000000000e0;
        dxr = 22.098018177947786e0;
        dyr = 18.916071411549062e0;
        pmsdnut2(mjd, dx, dy);
        REQUIRE(std::abs(dx - dxr) < PRECISION);
        REQUIRE(std::abs(dy - dyr) < PRECISION);
      }
      SUBCASE("Example  43: Checking for descripancies  > PRECISION") {
        mjd = 63590.925781250000000e0;
        dxr = -7.582614070993224e0;
        dyr = -5.651613344913026e0;
        pmsdnut2(mjd, dx, dy);
        REQUIRE(std::abs(dx - dxr) < PRECISION);
        REQUIRE(std::abs(dy - dyr) < PRECISION);
      }
      SUBCASE("Example  44: Checking for descripancies  > PRECISION") {
        mjd = 64645.445312500000000e0;
        dxr = 16.826167195125198e0;
        dyr = 11.570949654020081e0;
        pmsdnut2(mjd, dx, dy);
        REQUIRE(std::abs(dx - dxr) < PRECISION);
        REQUIRE(std::abs(dy - dyr) < PRECISION);
      }
      SUBCASE("Example  45: Checking for descripancies  > PRECISION") {
        mjd = 61146.011718750000000e0;
        dxr = 10.678617994316941e0;
        dyr = -0.458717990362933e0;
        pmsdnut2(mjd, dx, dy);
        REQUIRE(std::abs(dx - dxr) < PRECISION);
        REQUIRE(std::abs(dy - dyr) < PRECISION);
      }
      SUBCASE("Example  46: Checking for descripancies  > PRECISION") {
        mjd = 41258.757812500000000e0;
        dxr = 17.838406078366820e0;
        dyr = -29.604850696539462e0;
        pmsdnut2(mjd, dx, dy);
        REQUIRE(std::abs(dx - dxr) < PRECISION);
        REQUIRE(std::abs(dy - dyr) < PRECISION);
      }
      SUBCASE("Example  47: Checking for descripancies  > PRECISION") {
        mjd = 60463.003906250000000e0;
        dxr = 12.882113050681902e0;
        dyr = -1.033099895637409e0;
        pmsdnut2(mjd, dx, dy);
        REQUIRE(std::abs(dx - dxr) < PRECISION);
        REQUIRE(std::abs(dy - dyr) < PRECISION);
      }
      SUBCASE("Example  48: Checking for descripancies  > PRECISION") {
        mjd = 69158.804687500000000e0;
        dxr = -17.702404629184372e0;
        dyr = -8.515368874662554e0;
        pmsdnut2(mjd, dx, dy);
        REQUIRE(std::abs(dx - dxr) < PRECISION);
        REQUIRE(std::abs(dy - dyr) < PRECISION);
      }
      SUBCASE("Example  49: Checking for descripancies  > PRECISION") {
        mjd = 48130.093750000000000e0;
        dxr = -2.850101517863636e0;
        dyr = -18.400397382560808e0;
        pmsdnut2(mjd, dx, dy);
        REQUIRE(std::abs(dx - dxr) < PRECISION);
        REQUIRE(std::abs(dy - dyr) < PRECISION);
      }
      SUBCASE("Example  50: Checking for descripancies  > PRECISION") {
        mjd = 53755.183593750000000e0;
        dxr = -6.779783181243143e0;
        dyr = 4.321800693584171e0;
        pmsdnut2(mjd, dx, dy);
        REQUIRE(std::abs(dx - dxr) < PRECISION);
        REQUIRE(std::abs(dy - dyr) < PRECISION);
      }
      SUBCASE("Example  51: Checking for descripancies  > PRECISION") {
        mjd = 65096.144531250000000e0;
        dxr = 0.435240164984248e0;
        dyr = 14.766719145037863e0;
        pmsdnut2(mjd, dx, dy);
        REQUIRE(std::abs(dx - dxr) < PRECISION);
        REQUIRE(std::abs(dy - dyr) < PRECISION);
      }
      SUBCASE("Example  52: Checking for descripancies  > PRECISION") {
        mjd = 55175.988281250000000e0;
        dxr = -16.689978102145073e0;
        dyr = 1.688091067215500e0;
        pmsdnut2(mjd, dx, dy);
        REQUIRE(std::abs(dx - dxr) < PRECISION);
        REQUIRE(std::abs(dy - dyr) < PRECISION);
      }
      SUBCASE("Example  53: Checking for descripancies  > PRECISION") {
        mjd = 46930.585937500000000e0;
        dxr = -37.850651852764315e0;
        dyr = -14.797991790071515e0;
        pmsdnut2(mjd, dx, dy);
        REQUIRE(std::abs(dx - dxr) < PRECISION);
        REQUIRE(std::abs(dy - dyr) < PRECISION);
      }
      SUBCASE("Example  54: Checking for descripancies  > PRECISION") {
        mjd = 48206.882812500000000e0;
        dxr = -0.801384974968732e0;
        dyr = -5.413753167033020e0;
        pmsdnut2(mjd, dx, dy);
        REQUIRE(std::abs(dx - dxr) < PRECISION);
        REQUIRE(std::abs(dy - dyr) < PRECISION);
      }
      SUBCASE("Example  55: Checking for descripancies  > PRECISION") {
        mjd = 51104.894531250000000e0;
        dxr = -3.162905162832105e0;
        dyr = -4.262953191169372e0;
        pmsdnut2(mjd, dx, dy);
        REQUIRE(std::abs(dx - dxr) < PRECISION);
        REQUIRE(std::abs(dy - dyr) < PRECISION);
      }
      SUBCASE("Example  56: Checking for descripancies  > PRECISION") {
        mjd = 44481.339843750000000e0;
        dxr = -7.297482955754317e0;
        dyr = 4.309331512697030e0;
        pmsdnut2(mjd, dx, dy);
        REQUIRE(std::abs(dx - dxr) < PRECISION);
        REQUIRE(std::abs(dy - dyr) < PRECISION);
      }
      SUBCASE("Example  57: Checking for descripancies  > PRECISION") {
        mjd = 55477.058593750000000e0;
        dxr = -16.715666861946541e0;
        dyr = 0.242607833065277e0;
        pmsdnut2(mjd, dx, dy);
        REQUIRE(std::abs(dx - dxr) < PRECISION);
        REQUIRE(std::abs(dy - dyr) < PRECISION);
      }
      SUBCASE("Example  58: Checking for descripancies  > PRECISION") {
        mjd = 69603.921875000000000e0;
        dxr = 22.253811589307890e0;
        dyr = 9.376549231049733e0;
        pmsdnut2(mjd, dx, dy);
        REQUIRE(std::abs(dx - dxr) < PRECISION);
        REQUIRE(std::abs(dy - dyr) < PRECISION);
      }
      SUBCASE("Example  59: Checking for descripancies  > PRECISION") {
        mjd = 70000.132812500000000e0;
        dxr = 9.303948486265627e0;
        dyr = -21.168343032637281e0;
        pmsdnut2(mjd, dx, dy);
        REQUIRE(std::abs(dx - dxr) < PRECISION);
        REQUIRE(std::abs(dy - dyr) < PRECISION);
      }
      SUBCASE("Example  60: Checking for descripancies  > PRECISION") {
        mjd = 40841.382812500000000e0;
        dxr = -21.133599160127453e0;
        dyr = -15.758244963868266e0;
        pmsdnut2(mjd, dx, dy);
        REQUIRE(std::abs(dx - dxr) < PRECISION);
        REQUIRE(std::abs(dy - dyr) < PRECISION);
      }
      SUBCASE("Example  61: Checking for descripancies  > PRECISION") {
        mjd = 69843.429687500000000e0;
        dxr = 10.073690085401307e0;
        dyr = 14.402798592556453e0;
        pmsdnut2(mjd, dx, dy);
        REQUIRE(std::abs(dx - dxr) < PRECISION);
        REQUIRE(std::abs(dy - dyr) < PRECISION);
      }
      SUBCASE("Example  62: Checking for descripancies  > PRECISION") {
        mjd = 56095.546875000000000e0;
        dxr = -28.003759576189911e0;
        dyr = 1.112066812184740e0;
        pmsdnut2(mjd, dx, dy);
        REQUIRE(std::abs(dx - dxr) < PRECISION);
        REQUIRE(std::abs(dy - dyr) < PRECISION);
      }
      SUBCASE("Example  63: Checking for descripancies  > PRECISION") {
        mjd = 56500.027343750000000e0;
        dxr = 4.946193092754573e0;
        dyr = 0.817000446014318e0;
        pmsdnut2(mjd, dx, dy);
        REQUIRE(std::abs(dx - dxr) < PRECISION);
        REQUIRE(std::abs(dy - dyr) < PRECISION);
      }
      SUBCASE("Example  64: Checking for descripancies  > PRECISION") {
        mjd = 49722.156250000000000e0;
        dxr = -15.661377867550097e0;
        dyr = -0.340670462519145e0;
        pmsdnut2(mjd, dx, dy);
        REQUIRE(std::abs(dx - dxr) < PRECISION);
        REQUIRE(std::abs(dy - dyr) < PRECISION);
      }
      SUBCASE("Example  65: Checking for descripancies  > PRECISION") {
        mjd = 72660.031250000000000e0;
        dxr = 5.736141307635520e0;
        dyr = -7.089533251391602e0;
        pmsdnut2(mjd, dx, dy);
        REQUIRE(std::abs(dx - dxr) < PRECISION);
        REQUIRE(std::abs(dy - dyr) < PRECISION);
      }
      SUBCASE("Example  66: Checking for descripancies  > PRECISION") {
        mjd = 55733.132812500000000e0;
        dxr = 12.950069292899482e0;
        dyr = 1.952038217605005e0;
        pmsdnut2(mjd, dx, dy);
        REQUIRE(std::abs(dx - dxr) < PRECISION);
        REQUIRE(std::abs(dy - dyr) < PRECISION);
      }
      SUBCASE("Example  67: Checking for descripancies  > PRECISION") {
        mjd = 47905.195312500000000e0;
        dxr = -18.361600286145869e0;
        dyr = 4.952180857761604e0;
        pmsdnut2(mjd, dx, dy);
        REQUIRE(std::abs(dx - dxr) < PRECISION);
        REQUIRE(std::abs(dy - dyr) < PRECISION);
      }
      SUBCASE("Example  68: Checking for descripancies  > PRECISION") {
        mjd = 41878.531250000000000e0;
        dxr = -25.966770058894092e0;
        dyr = -8.691607171643376e0;
        pmsdnut2(mjd, dx, dy);
        REQUIRE(std::abs(dx - dxr) < PRECISION);
        REQUIRE(std::abs(dy - dyr) < PRECISION);
      }
      SUBCASE("Example  69: Checking for descripancies  > PRECISION") {
        mjd = 71324.281250000000000e0;
        dxr = 6.509319061234912e0;
        dyr = -0.930842548356267e0;
        pmsdnut2(mjd, dx, dy);
        REQUIRE(std::abs(dx - dxr) < PRECISION);
        REQUIRE(std::abs(dy - dyr) < PRECISION);
      }
      SUBCASE("Example  70: Checking for descripancies  > PRECISION") {
        mjd = 41294.582031250000000e0;
        dxr = 8.783278218248695e0;
        dyr = 2.613112826527714e0;
        pmsdnut2(mjd, dx, dy);
        REQUIRE(std::abs(dx - dxr) < PRECISION);
        REQUIRE(std::abs(dy - dyr) < PRECISION);
      }
      SUBCASE("Example  71: Checking for descripancies  > PRECISION") {
        mjd = 55964.296875000000000e0;
        dxr = -5.251264910642107e0;
        dyr = 21.928814824772275e0;
        pmsdnut2(mjd, dx, dy);
        REQUIRE(std::abs(dx - dxr) < PRECISION);
        REQUIRE(std::abs(dy - dyr) < PRECISION);
      }
      SUBCASE("Example  72: Checking for descripancies  > PRECISION") {
        mjd = 51959.039062500000000e0;
        dxr = -22.855760395927462e0;
        dyr = 7.092530022622603e0;
        pmsdnut2(mjd, dx, dy);
        REQUIRE(std::abs(dx - dxr) < PRECISION);
        REQUIRE(std::abs(dy - dyr) < PRECISION);
      }
      SUBCASE("Example  73: Checking for descripancies  > PRECISION") {
        mjd = 48281.820312500000000e0;
        dxr = -23.971303314348010e0;
        dyr = -26.105198299506547e0;
        pmsdnut2(mjd, dx, dy);
        REQUIRE(std::abs(dx - dxr) < PRECISION);
        REQUIRE(std::abs(dy - dyr) < PRECISION);
      }
      SUBCASE("Example  74: Checking for descripancies  > PRECISION") {
        mjd = 70158.109375000000000e0;
        dxr = -20.630597323088693e0;
        dyr = -8.202158517106499e0;
        pmsdnut2(mjd, dx, dy);
        REQUIRE(std::abs(dx - dxr) < PRECISION);
        REQUIRE(std::abs(dy - dyr) < PRECISION);
      }
      SUBCASE("Example  75: Checking for descripancies  > PRECISION") {
        mjd = 56962.589843750000000e0;
        dxr = 2.702019351417760e0;
        dyr = 9.672240253901050e0;
        pmsdnut2(mjd, dx, dy);
        REQUIRE(std::abs(dx - dxr) < PRECISION);
        REQUIRE(std::abs(dy - dyr) < PRECISION);
      }
      SUBCASE("Example  76: Checking for descripancies  > PRECISION") {
        mjd = 54718.835937500000000e0;
        dxr = 21.056732026852060e0;
        dyr = 10.679598249858994e0;
        pmsdnut2(mjd, dx, dy);
        REQUIRE(std::abs(dx - dxr) < PRECISION);
        REQUIRE(std::abs(dy - dyr) < PRECISION);
      }
      SUBCASE("Example  77: Checking for descripancies  > PRECISION") {
        mjd = 71091.835937500000000e0;
        dxr = 4.731916779332427e0;
        dyr = 8.099831049938961e0;
        pmsdnut2(mjd, dx, dy);
        REQUIRE(std::abs(dx - dxr) < PRECISION);
        REQUIRE(std::abs(dy - dyr) < PRECISION);
      }
      SUBCASE("Example  78: Checking for descripancies  > PRECISION") {
        mjd = 40481.101562500000000e0;
        dxr = -15.745423309021476e0;
        dyr = -22.098121857600258e0;
        pmsdnut2(mjd, dx, dy);
        REQUIRE(std::abs(dx - dxr) < PRECISION);
        REQUIRE(std::abs(dy - dyr) < PRECISION);
      }
      SUBCASE("Example  79: Checking for descripancies  > PRECISION") {
        mjd = 64925.863281250000000e0;
        dxr = 7.840138311274127e0;
        dyr = 24.535689852633936e0;
        pmsdnut2(mjd, dx, dy);
        REQUIRE(std::abs(dx - dxr) < PRECISION);
        REQUIRE(std::abs(dy - dyr) < PRECISION);
      }
      SUBCASE("Example  80: Checking for descripancies  > PRECISION") {
        mjd = 65224.281250000000000e0;
        dxr = 18.806110856089990e0;
        dyr = -28.955611420380741e0;
        pmsdnut2(mjd, dx, dy);
        REQUIRE(std::abs(dx - dxr) < PRECISION);
        REQUIRE(std::abs(dy - dyr) < PRECISION);
      }
      SUBCASE("Example  81: Checking for descripancies  > PRECISION") {
        mjd = 67203.734375000000000e0;
        dxr = 28.640332401621116e0;
        dyr = -18.533995447993849e0;
        pmsdnut2(mjd, dx, dy);
        REQUIRE(std::abs(dx - dxr) < PRECISION);
        REQUIRE(std::abs(dy - dyr) < PRECISION);
      }
      SUBCASE("Example  82: Checking for descripancies  > PRECISION") {
        mjd = 43068.066406250000000e0;
        dxr = 0.273374050199649e0;
        dyr = -15.813701151568903e0;
        pmsdnut2(mjd, dx, dy);
        REQUIRE(std::abs(dx - dxr) < PRECISION);
        REQUIRE(std::abs(dy - dyr) < PRECISION);
      }
      SUBCASE("Example  83: Checking for descripancies  > PRECISION") {
        mjd = 39306.859375000000000e0;
        dxr = 15.780702361589208e0;
        dyr = 31.331860840446666e0;
        pmsdnut2(mjd, dx, dy);
        REQUIRE(std::abs(dx - dxr) < PRECISION);
        REQUIRE(std::abs(dy - dyr) < PRECISION);
      }
      SUBCASE("Example  84: Checking for descripancies  > PRECISION") {
        mjd = 62415.699218750000000e0;
        dxr = 1.801900085507055e0;
        dyr = -7.636268845384032e0;
        pmsdnut2(mjd, dx, dy);
        REQUIRE(std::abs(dx - dxr) < PRECISION);
        REQUIRE(std::abs(dy - dyr) < PRECISION);
      }
      SUBCASE("Example  85: Checking for descripancies  > PRECISION") {
        mjd = 68592.125000000000000e0;
        dxr = -14.141901725537833e0;
        dyr = -6.912269990587727e0;
        pmsdnut2(mjd, dx, dy);
        REQUIRE(std::abs(dx - dxr) < PRECISION);
        REQUIRE(std::abs(dy - dyr) < PRECISION);
      }
      SUBCASE("Example  86: Checking for descripancies  > PRECISION") {
        mjd = 60391.707031250000000e0;
        dxr = 4.307563696254488e0;
        dyr = -9.945641533252617e0;
        pmsdnut2(mjd, dx, dy);
        REQUIRE(std::abs(dx - dxr) < PRECISION);
        REQUIRE(std::abs(dy - dyr) < PRECISION);
      }
      SUBCASE("Example  87: Checking for descripancies  > PRECISION") {
        mjd = 64056.351562500000000e0;
        dxr = -1.366276705250699e0;
        dyr = -7.975768824746169e0;
        pmsdnut2(mjd, dx, dy);
        REQUIRE(std::abs(dx - dxr) < PRECISION);
        REQUIRE(std::abs(dy - dyr) < PRECISION);
      }
      SUBCASE("Example  88: Checking for descripancies  > PRECISION") {
        mjd = 63685.226562500000000e0;
        dxr = -2.756077875199102e0;
        dyr = 0.953578400843045e0;
        pmsdnut2(mjd, dx, dy);
        REQUIRE(std::abs(dx - dxr) < PRECISION);
        REQUIRE(std::abs(dy - dyr) < PRECISION);
      }
      SUBCASE("Example  89: Checking for descripancies  > PRECISION") {
        mjd = 73101.156250000000000e0;
        dxr = -4.889012309749827e0;
        dyr = 13.753427808963997e0;
        pmsdnut2(mjd, dx, dy);
        REQUIRE(std::abs(dx - dxr) < PRECISION);
        REQUIRE(std::abs(dy - dyr) < PRECISION);
      }
      SUBCASE("Example  90: Checking for descripancies  > PRECISION") {
        mjd = 69291.195312500000000e0;
        dxr = -2.557822292507256e0;
        dyr = -22.756618728025106e0;
        pmsdnut2(mjd, dx, dy);
        REQUIRE(std::abs(dx - dxr) < PRECISION);
        REQUIRE(std::abs(dy - dyr) < PRECISION);
      }
      SUBCASE("Example  91: Checking for descripancies  > PRECISION") {
        mjd = 46773.734375000000000e0;
        dxr = 4.065085129745244e0;
        dyr = -7.239554047629208e0;
        pmsdnut2(mjd, dx, dy);
        REQUIRE(std::abs(dx - dxr) < PRECISION);
        REQUIRE(std::abs(dy - dyr) < PRECISION);
      }
      SUBCASE("Example  92: Checking for descripancies  > PRECISION") {
        mjd = 49285.406250000000000e0;
        dxr = -0.903391923474148e0;
        dyr = 3.282453829563390e0;
        pmsdnut2(mjd, dx, dy);
        REQUIRE(std::abs(dx - dxr) < PRECISION);
        REQUIRE(std::abs(dy - dyr) < PRECISION);
      }
      SUBCASE("Example  93: Checking for descripancies  > PRECISION") {
        mjd = 50821.359375000000000e0;
        dxr = 17.415077183128243e0;
        dyr = 16.507971340769142e0;
        pmsdnut2(mjd, dx, dy);
        REQUIRE(std::abs(dx - dxr) < PRECISION);
        REQUIRE(std::abs(dy - dyr) < PRECISION);
      }
      SUBCASE("Example  94: Checking for descripancies  > PRECISION") {
        mjd = 56396.703125000000000e0;
        dxr = -20.360383101493117e0;
        dyr = 1.631323847895962e0;
        pmsdnut2(mjd, dx, dy);
        REQUIRE(std::abs(dx - dxr) < PRECISION);
        REQUIRE(std::abs(dy - dyr) < PRECISION);
      }
      SUBCASE("Example  95: Checking for descripancies  > PRECISION") {
        mjd = 59070.339843750000000e0;
        dxr = -2.481036832415603e0;
        dyr = -4.418815416174116e0;
        pmsdnut2(mjd, dx, dy);
        REQUIRE(std::abs(dx - dxr) < PRECISION);
        REQUIRE(std::abs(dy - dyr) < PRECISION);
      }
      SUBCASE("Example  96: Checking for descripancies  > PRECISION") {
        mjd = 67827.226562500000000e0;
        dxr = -9.378038925079327e0;
        dyr = -12.937342239309576e0;
        pmsdnut2(mjd, dx, dy);
        REQUIRE(std::abs(dx - dxr) < PRECISION);
        REQUIRE(std::abs(dy - dyr) < PRECISION);
      }
      SUBCASE("Example  97: Checking for descripancies  > PRECISION") {
        mjd = 52919.500000000000000e0;
        dxr = -9.117092953984242e0;
        dyr = 0.454560376956332e0;
        pmsdnut2(mjd, dx, dy);
        REQUIRE(std::abs(dx - dxr) < PRECISION);
        REQUIRE(std::abs(dy - dyr) < PRECISION);
      }
      SUBCASE("Example  98: Checking for descripancies  > PRECISION") {
        mjd = 67674.984375000000000e0;
        dxr = -2.370789919308554e0;
        dyr = -5.153603232417709e0;
        pmsdnut2(mjd, dx, dy);
        REQUIRE(std::abs(dx - dxr) < PRECISION);
        REQUIRE(std::abs(dy - dyr) < PRECISION);
      }
      SUBCASE("Example  99: Checking for descripancies  > PRECISION") {
        mjd = 48014.281250000000000e0;
        dxr = 9.605335526428581e0;
        dyr = 1.549345122942371e0;
        pmsdnut2(mjd, dx, dy);
        REQUIRE(std::abs(dx - dxr) < PRECISION);
        REQUIRE(std::abs(dy - dyr) < PRECISION);
      }
      SUBCASE("Example 100: Checking for descripancies  > PRECISION") {
        mjd = 53033.105468750000000e0;
        dxr = -5.435957385487296e0;
        dyr = 12.899641470466372e0;
        pmsdnut2(mjd, dx, dy);
        REQUIRE(std::abs(dx - dxr) < PRECISION);
        REQUIRE(std::abs(dy - dyr) < PRECISION);
      }
  }
