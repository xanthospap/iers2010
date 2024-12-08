#include "pole_tide.hpp"

using namespace dso;

int main(int argc, char *argv[]) {
  if (argc != 2) {
    fprintf(stderr, "Usage: %s [desaiscopolecoef.txt]\n", argv[0]);
    return 1;
  }

  /* storing coeffs in range [0,360], so size must be MAX_DEGREE_DESAI_2002+1 */
  CoeffMatrix2D<MatrixStorageType::LwTriangularColWise> A_real(
      pole_tide_details::MAX_DEGREE_DESAI_2002+1,
      pole_tide_details::MAX_DEGREE_DESAI_2002+1);
  CoeffMatrix2D<MatrixStorageType::LwTriangularColWise> A_imag(
      pole_tide_details::MAX_DEGREE_DESAI_2002+1,
      pole_tide_details::MAX_DEGREE_DESAI_2002+1);
  CoeffMatrix2D<MatrixStorageType::LwTriangularColWise> B_real(
      pole_tide_details::MAX_DEGREE_DESAI_2002+1,
      pole_tide_details::MAX_DEGREE_DESAI_2002+1);
  CoeffMatrix2D<MatrixStorageType::LwTriangularColWise> B_imag(
      pole_tide_details::MAX_DEGREE_DESAI_2002+1,
      pole_tide_details::MAX_DEGREE_DESAI_2002+1);

  assert(!pole_tide_details::parse_desai02_coeffs(
      argv[1], pole_tide_details::MAX_DEGREE_DESAI_2002,
      pole_tide_details::MAX_DEGREE_DESAI_2002, A_real, A_imag, B_real,
      B_imag));

  /* check random entries */
assert((std::abs(A_real(0,0)- 0e0)<1e-13) && (std::abs(A_imag(0,0)-0e0)<1e-13));
assert((std::abs(B_real(0,0)- 0e0)<1e-13) && (std::abs(B_imag(0,0)-0e0)<1e-13));
assert((std::abs(A_real(1,0)- +1.8736759805448e-02)<1e-13) && (std::abs(A_imag(1,0)- +2.9688884960424e-02)<1e-13));
assert((std::abs(B_real(1,0)- +0.0000000000000e+00)<1e-13) && (std::abs(B_imag(1,0)- +0.0000000000000e+00)<1e-13));
assert((std::abs(A_real(1,1)- +2.8258913146935e-02)<1e-13) && (std::abs(A_imag(1,1)- +2.3898264393684e-02)<1e-13));
assert((std::abs(B_real(1,1)- +2.1774643075236e-02)<1e-13) && (std::abs(B_imag(1,1)- +5.6771602236635e-02)<1e-13));
assert((std::abs(A_real(2,0)- -3.9555099024374e-03)<1e-13) && (std::abs(A_imag(2,0)- +6.8390464271953e-04)<1e-13));
assert((std::abs(B_real(2,0)- +0.0000000000000e+00)<1e-13) && (std::abs(B_imag(2,0)- +0.0000000000000e+00)<1e-13));
assert((std::abs(A_real(2,1)- -2.4325330521304e-01)<1e-13) && (std::abs(A_imag(2,1)- +5.4680741193318e-03)<1e-13));
assert((std::abs(B_real(2,1)- +5.4680741193318e-03)<1e-13) && (std::abs(B_imag(2,1)- -1.9252111185300e-01)<1e-13));
assert((std::abs(A_real(2,2)- +1.9102047023374e-02)<1e-13) && (std::abs(A_imag(2,2)- -1.5123770169928e-02)<1e-13));
assert((std::abs(B_real(2,2)- +1.1158297399424e-02)<1e-13) && (std::abs(B_imag(2,2)- -2.4857839911518e-04)<1e-13));
assert((std::abs(A_real(12,3)- -2.9726928754481e-03)<1e-13) && (std::abs(A_imag(12,3)- -9.5509318510813e-04)<1e-13));
assert((std::abs(B_real(12,3)- -4.7318070752733e-03)<1e-13) && (std::abs(B_imag(12,3)- -3.4045563384392e-03)<1e-13));
assert((std::abs(A_real(26,25)- +3.3733426746577e-04)<1e-13) && (std::abs(A_imag(26,25)- -1.2876197892468e-03)<1e-13));
assert((std::abs(B_real(26,25)- +1.6331957652733e-03)<1e-13) && (std::abs(B_imag(26,25)- +1.2804341141451e-03)<1e-13));
assert((std::abs(A_real(59,40)- +1.1862431591654e-04)<1e-13) && (std::abs(A_imag(59,40)- +2.3137977106912e-04)<1e-13));
assert((std::abs(B_real(59,40)- +1.1887476883959e-05)<1e-13) && (std::abs(B_imag(59,40)- -1.2275808014139e-04)<1e-13));
assert((std::abs(A_real(62,59)- -1.9883896499429e-04)<1e-13) && (std::abs(A_imag(62,59)- -7.0142310388714e-05)<1e-13));
assert((std::abs(B_real(62,59)- +1.2726881474527e-04)<1e-13) && (std::abs(B_imag(62,59)- +8.9533047008726e-05)<1e-13));
assert((std::abs(A_real(68,8)- -1.5634459171970e-04)<1e-13) && (std::abs(A_imag(68,8)- -3.4122274274536e-04)<1e-13));
assert((std::abs(B_real(68,8)- +3.7082524545487e-04)<1e-13) && (std::abs(B_imag(68,8)- +2.7231462304964e-04)<1e-13));
assert((std::abs(A_real(156,119)- +1.1850296243836e-04)<1e-13) && (std::abs(A_imag(156,119)- -1.0196015957226e-04)<1e-13));
assert((std::abs(B_real(156,119)- +2.0459680023565e-05)<1e-13) && (std::abs(B_imag(156,119)- -1.1340095237695e-05)<1e-13));
assert((std::abs(A_real(157,86)- +1.1159614213762e-04)<1e-13) && (std::abs(A_imag(157,86)- +6.2768150097535e-06)<1e-13));
assert((std::abs(B_real(157,86)- +1.5037245139985e-04)<1e-13) && (std::abs(B_imag(157,86)- +1.4135356418098e-04)<1e-13));
assert((std::abs(A_real(192,81)- +5.8203193726805e-05)<1e-13) && (std::abs(A_imag(192,81)- -1.8268604245171e-04)<1e-13));
assert((std::abs(B_real(192,81)- -9.6989830944624e-05)<1e-13) && (std::abs(B_imag(192,81)- +1.8177201838200e-06)<1e-13));
assert((std::abs(A_real(193,138)- -6.2920487442486e-05)<1e-13) && (std::abs(A_imag(193,138)- +6.2657124574944e-05)<1e-13));
assert((std::abs(B_real(193,138)- -6.2109742928680e-05)<1e-13) && (std::abs(B_imag(193,138)- +1.0893152959543e-04)<1e-13));
assert((std::abs(A_real(193,176)- +4.5866450774686e-05)<1e-13) && (std::abs(A_imag(193,176)- +2.2623386719646e-05)<1e-13));
assert((std::abs(B_real(193,176)- +5.1523588449046e-05)<1e-13) && (std::abs(B_imag(193,176)- +6.2347285907247e-05)<1e-13));
assert((std::abs(A_real(195,108)- -5.1097560280963e-06)<1e-13) && (std::abs(A_imag(195,108)- -7.0457313034205e-05)<1e-13));
assert((std::abs(B_real(195,108)- -1.2163749695989e-05)<1e-13) && (std::abs(B_imag(195,108)- -6.8887463784543e-05)<1e-13));
assert((std::abs(A_real(197,32)- +8.9047633096001e-05)<1e-13) && (std::abs(A_imag(197,32)- +1.0201954531690e-04)<1e-13));
assert((std::abs(B_real(197,32)- +7.5894981652281e-05)<1e-13) && (std::abs(B_imag(197,32)- +4.4048730590866e-06)<1e-13));
assert((std::abs(A_real(209,2)- -1.1024834393314e-04)<1e-13) && (std::abs(A_imag(209,2)- +4.6863098489110e-05)<1e-13));
assert((std::abs(B_real(209,2)- -1.4544343459458e-05)<1e-13) && (std::abs(B_imag(209,2)- +1.3718881608936e-04)<1e-13));
assert((std::abs(A_real(211,44)- -4.2920401442671e-05)<1e-13) && (std::abs(A_imag(211,44)- +1.3323062710954e-05)<1e-13));
assert((std::abs(B_real(211,44)- +5.7950158814853e-05)<1e-13) && (std::abs(B_imag(211,44)- -3.5812186007931e-06)<1e-13));
assert((std::abs(A_real(221,27)- -4.8053395370844e-05)<1e-13) && (std::abs(A_imag(221,27)- -1.9044750107751e-05)<1e-13));
assert((std::abs(B_real(221,27)- -3.0797963569041e-05)<1e-13) && (std::abs(B_imag(221,27)- -1.0426587565770e-04)<1e-13));
assert((std::abs(A_real(223,6)- -2.3356671566528e-05)<1e-13) && (std::abs(A_imag(223,6)- +2.5737265329997e-05)<1e-13));
assert((std::abs(B_real(223,6)- +5.3369710510738e-06)<1e-13) && (std::abs(B_imag(223,6)- +4.9445087386581e-06)<1e-13));
assert((std::abs(A_real(225,198)- -5.0018288921820e-05)<1e-13) && (std::abs(A_imag(225,198)- -4.1030385569721e-05)<1e-13));
assert((std::abs(B_real(225,198)- -3.0637880257398e-05)<1e-13) && (std::abs(B_imag(225,198)- -5.9290211052641e-05)<1e-13));
assert((std::abs(A_real(234,170)- +1.8254661505124e-05)<1e-13) && (std::abs(A_imag(234,170)- +1.0966739643121e-04)<1e-13));
assert((std::abs(B_real(234,170)- +4.2236174384841e-05)<1e-13) && (std::abs(B_imag(234,170)- -4.9791073898843e-05)<1e-13));
assert((std::abs(A_real(242,240)- -6.1967733539319e-06)<1e-13) && (std::abs(A_imag(242,240)- -1.4611475170792e-05)<1e-13));
assert((std::abs(B_real(242,240)- -1.6566040766178e-05)<1e-13) && (std::abs(B_imag(242,240)- +8.8794598184049e-06)<1e-13));
assert((std::abs(A_real(247,101)- +6.1642402938237e-05)<1e-13) && (std::abs(A_imag(247,101)- -4.2853269731336e-05)<1e-13));
assert((std::abs(B_real(247,101)- -1.3185184503829e-04)<1e-13) && (std::abs(B_imag(247,101)- -5.7258185790330e-05)<1e-13));
assert((std::abs(A_real(248,102)- -8.9424603475032e-05)<1e-13) && (std::abs(A_imag(248,102)- -8.3896565245545e-05)<1e-13));
assert((std::abs(B_real(248,102)- -9.7513974690810e-05)<1e-13) && (std::abs(B_imag(248,102)- -6.4301851502633e-05)<1e-13));
assert((std::abs(A_real(248,133)- +1.8814953095876e-05)<1e-13) && (std::abs(A_imag(248,133)- +7.4673700903776e-05)<1e-13));
assert((std::abs(B_real(248,133)- +5.0245503944798e-06)<1e-13) && (std::abs(B_imag(248,133)- +4.4929188804904e-05)<1e-13));
assert((std::abs(A_real(255,214)- +7.4430993947414e-05)<1e-13) && (std::abs(A_imag(255,214)- +7.5321780117042e-05)<1e-13));
assert((std::abs(B_real(255,214)- -3.1500264549154e-06)<1e-13) && (std::abs(B_imag(255,214)- -4.0196403155838e-05)<1e-13));
assert((std::abs(A_real(258,203)- -9.3881116312068e-06)<1e-13) && (std::abs(A_imag(258,203)- +1.9932954288481e-05)<1e-13));
assert((std::abs(B_real(258,203)- -2.4146760169155e-05)<1e-13) && (std::abs(B_imag(258,203)- -8.2073681083930e-05)<1e-13));
assert((std::abs(A_real(259,165)- -2.0645583224183e-05)<1e-13) && (std::abs(A_imag(259,165)- -7.1043424513931e-05)<1e-13));
assert((std::abs(B_real(259,165)- +1.5494515373727e-05)<1e-13) && (std::abs(B_imag(259,165)- +4.8939017562677e-05)<1e-13));
assert((std::abs(A_real(266,171)- +2.2024365295433e-06)<1e-13) && (std::abs(A_imag(266,171)- -1.5078815161234e-05)<1e-13));
assert((std::abs(B_real(266,171)- -5.2453104452136e-05)<1e-13) && (std::abs(B_imag(266,171)- -2.3374745488810e-05)<1e-13));
assert((std::abs(A_real(273,85)- +6.1812814205470e-05)<1e-13) && (std::abs(A_imag(273,85)- -3.3333476666358e-06)<1e-13));
assert((std::abs(B_real(273,85)- +1.8489622494794e-05)<1e-13) && (std::abs(B_imag(273,85)- +9.5376087883704e-05)<1e-13));
assert((std::abs(A_real(276,111)- -4.4363733236357e-05)<1e-13) && (std::abs(A_imag(276,111)- +6.3390409421826e-05)<1e-13));
assert((std::abs(B_real(276,111)- +3.8443232694455e-05)<1e-13) && (std::abs(B_imag(276,111)- +3.0867677269273e-05)<1e-13));
assert((std::abs(A_real(280,179)- +2.1381450059962e-05)<1e-13) && (std::abs(A_imag(280,179)- -1.3020970009949e-05)<1e-13));
assert((std::abs(B_real(280,179)- +1.6695140211051e-05)<1e-13) && (std::abs(B_imag(280,179)- -1.6785509822463e-05)<1e-13));
assert((std::abs(A_real(281,202)- -1.8272735427407e-05)<1e-13) && (std::abs(A_imag(281,202)- -2.4193235969309e-06)<1e-13));
assert((std::abs(B_real(281,202)- +3.5182180716948e-05)<1e-13) && (std::abs(B_imag(281,202)- -9.3548039744294e-05)<1e-13));
assert((std::abs(A_real(282,179)- +5.9553573920183e-05)<1e-13) && (std::abs(A_imag(282,179)- -4.5393946806448e-05)<1e-13));
assert((std::abs(B_real(282,179)- -5.3395530399713e-06)<1e-13) && (std::abs(B_imag(282,179)- +8.8742958406099e-06)<1e-13));
assert((std::abs(A_real(284,184)- -9.5636868762348e-06)<1e-13) && (std::abs(A_imag(284,184)- -3.8084336237706e-06)<1e-13));
assert((std::abs(B_real(284,184)- +5.9811772809930e-05)<1e-13) && (std::abs(B_imag(284,184)- -2.9914923366106e-05)<1e-13));
assert((std::abs(A_real(286,169)- -5.0353947478491e-05)<1e-13) && (std::abs(A_imag(286,169)- +1.2199897540553e-05)<1e-13));
assert((std::abs(B_real(286,169)- +7.2164007695434e-05)<1e-13) && (std::abs(B_imag(286,169)- -4.2465971396716e-05)<1e-13));
assert((std::abs(A_real(293,253)- -6.9403861564037e-06)<1e-13) && (std::abs(A_imag(293,253)- -2.2553673019415e-05)<1e-13));
assert((std::abs(B_real(293,253)- -1.4126061258067e-05)<1e-13) && (std::abs(B_imag(293,253)- +1.4456555836944e-05)<1e-13));
assert((std::abs(A_real(294,282)- +2.9202294430404e-05)<1e-13) && (std::abs(A_imag(294,282)- -2.1223948737217e-05)<1e-13));
assert((std::abs(B_real(294,282)- +5.7578352413617e-07)<1e-13) && (std::abs(B_imag(294,282)- +2.6439566310109e-05)<1e-13));
assert((std::abs(A_real(299,53)- +1.1859333204977e-04)<1e-13) && (std::abs(A_imag(299,53)- -6.0674157101685e-05)<1e-13));
assert((std::abs(B_real(299,53)- +2.1908175795047e-05)<1e-13) && (std::abs(B_imag(299,53)- +4.1415101915445e-06)<1e-13));
assert((std::abs(A_real(303,15)- -6.5352720983904e-05)<1e-13) && (std::abs(A_imag(303,15)- +3.9233321297081e-05)<1e-13));
assert((std::abs(B_real(303,15)- +4.3811187612298e-05)<1e-13) && (std::abs(B_imag(303,15)- +3.4934309782791e-06)<1e-13));
assert((std::abs(A_real(316,27)- +2.1590591160658e-05)<1e-13) && (std::abs(A_imag(316,27)- -5.5164192794481e-05)<1e-13));
assert((std::abs(B_real(316,27)- +4.2158714717916e-05)<1e-13) && (std::abs(B_imag(316,27)- +4.8799590162213e-06)<1e-13));
assert((std::abs(A_real(316,213)- +1.6860387500490e-05)<1e-13) && (std::abs(A_imag(316,213)- -1.3312734755364e-05)<1e-13));
assert((std::abs(B_real(316,213)- -4.5995550219559e-05)<1e-13) && (std::abs(B_imag(316,213)- +9.8812715106239e-06)<1e-13));
assert((std::abs(A_real(325,213)- -2.9834881710284e-05)<1e-13) && (std::abs(A_imag(325,213)- +5.3176061343428e-05)<1e-13));
assert((std::abs(B_real(325,213)- +1.0534514701684e-05)<1e-13) && (std::abs(B_imag(325,213)- +3.3762086823339e-05)<1e-13));
assert((std::abs(A_real(328,326)- -6.9455501341986e-07)<1e-13) && (std::abs(A_imag(328,326)- -5.0813541502241e-06)<1e-13));
assert((std::abs(B_real(328,326)- -1.3013224043257e-06)<1e-13) && (std::abs(B_imag(328,326)- -3.9504430131316e-07)<1e-13));
assert((std::abs(A_real(330,183)- -2.5579764007251e-05)<1e-13) && (std::abs(A_imag(330,183)- +1.2499408252660e-05)<1e-13));
assert((std::abs(B_real(330,183)- +2.5096932083424e-05)<1e-13) && (std::abs(B_imag(330,183)- -3.5158942816871e-05)<1e-13));
assert((std::abs(A_real(335,247)- +1.3674213201010e-05)<1e-13) && (std::abs(A_imag(335,247)- -2.0152423344766e-05)<1e-13));
assert((std::abs(B_real(335,247)- +5.5763214218694e-05)<1e-13) && (std::abs(B_imag(335,247)- +2.6553173821598e-05)<1e-13));
assert((std::abs(A_real(336,283)- -9.0782729635102e-06)<1e-13) && (std::abs(A_imag(336,283)- -1.3160397410413e-05)<1e-13));
assert((std::abs(B_real(336,283)- -7.5962113363615e-06)<1e-13) && (std::abs(B_imag(336,283)- +5.4841444182297e-06)<1e-13));
assert((std::abs(A_real(338,179)- +1.7567852458397e-05)<1e-13) && (std::abs(A_imag(338,179)- +2.6362106643185e-06)<1e-13));
assert((std::abs(B_real(338,179)- -5.7629992191566e-06)<1e-13) && (std::abs(B_imag(338,179)- +1.1185175903847e-05)<1e-13));
assert((std::abs(A_real(346,202)- -7.6039573030539e-06)<1e-13) && (std::abs(A_imag(346,202)- -2.4425614970081e-06)<1e-13));
assert((std::abs(B_real(346,202)- +2.3496098966992e-07)<1e-13) && (std::abs(B_imag(346,202)- -3.6483503121274e-06)<1e-13));
assert((std::abs(A_real(349,16)- -1.4317128436480e-05)<1e-13) && (std::abs(A_imag(349,16)- -5.4325046402803e-06)<1e-13));
assert((std::abs(B_real(349,16)- +2.7663429047511e-05)<1e-13) && (std::abs(B_imag(349,16)- -4.0276444811660e-05)<1e-13));
assert((std::abs(A_real(349,152)- +2.9638650554077e-05)<1e-13) && (std::abs(A_imag(349,152)- -3.3444502188950e-05)<1e-13));
assert((std::abs(B_real(349,152)- -5.1877320928000e-05)<1e-13) && (std::abs(B_imag(349,152)- -1.8724583484464e-05)<1e-13));
assert((std::abs(A_real(356,129)- +8.1177730936507e-06)<1e-13) && (std::abs(A_imag(356,129)- +1.7183136112674e-05)<1e-13));
assert((std::abs(B_real(356,129)- -7.4077261468461e-06)<1e-13) && (std::abs(B_imag(356,129)- -1.1897804839549e-05)<1e-13));
assert((std::abs(A_real(360,0)- +0.0000000000000e+00)<1e-13) && (std::abs(A_imag(360,0)- +0.0000000000000e+00)<1e-13));
assert((std::abs(B_real(360,0)- +0.0000000000000e+00)<1e-13) && (std::abs(B_imag(360,0)- +0.0000000000000e+00)<1e-13));
assert((std::abs(A_real(360,359)- +2.0759190735672e-06)<1e-13) && (std::abs(A_imag(360,359)- +2.8814996095783e-06)<1e-13));
assert((std::abs(B_real(360,359)- -5.5461121517690e-06)<1e-13) && (std::abs(B_imag(360,359)- -2.4942012749576e-06)<1e-13));
assert((std::abs(A_real(360,360)- +0.0000000000000e+00)<1e-13) && (std::abs(A_imag(360,360)- +0.0000000000000e+00)<1e-13));
assert((std::abs(B_real(360,360)- +0.0000000000000e+00)<1e-13) && (std::abs(B_imag(360,360)- +0.0000000000000e+00)<1e-13));

return 0;
}
