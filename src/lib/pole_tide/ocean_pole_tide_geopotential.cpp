#include "geodesy/units.hpp"
#include "iersconst.hpp"
#include "pole_tide.hpp"

namespace {
constexpr const double kn[dso::pole_tide_details::MAX_DEGREE_DESAI_2002 + 1] = {
    /*0*/ -1.0000000000e+00,   /*1*/ -0.1000000000e+01,
    /*2*/ -0.3054020195e+00,   /*3*/ -0.1960294041e+00,
    /*4*/ -0.1336652689e+00,   /*5*/ -0.1047066267e+00,
    /*6*/ -0.9033564429e-01,   /*7*/ -0.8206984804e-01,
    /*8*/ -0.7655494644e-01,   /*9*/ -0.7243844815e-01,
    /*10*/ -0.6913401466e-01,  /*11*/ -0.6635869819e-01,
    /*12*/ -0.6395689877e-01,  /*13*/ -0.6183296641e-01,
    /*14*/ -0.5992172201e-01,  /*15*/ -0.5817772516e-01,
    /*16*/ -0.5656704205e-01,  /*17*/ -0.5506474110e-01,
    /*18*/ -0.5365205058e-01,  /*19*/ -0.5231479422e-01,
    /*20*/ -0.5104204279e-01,  /*21*/ -0.4982562670e-01,
    /*22*/ -0.4865919131e-01,  /*23*/ -0.4753764652e-01,
    /*24*/ -0.4645728556e-01,  /*25*/ -0.4541483359e-01,
    /*26*/ -0.4440818528e-01,  /*27*/ -0.4343487603e-01,
    /*28*/ -0.4249347823e-01,  /*29*/ -0.4158232061e-01,
    /*30*/ -0.4070034557e-01,  /*31*/ -0.3984645832e-01,
    /*32*/ -0.3901937759e-01,  /*33*/ -0.3821827509e-01,
    /*34*/ -0.3744217649e-01,  /*35*/ -0.3669034206e-01,
    /*36*/ -0.3596186474e-01,  /*37*/ -0.3525594412e-01,
    /*38*/ -0.3457187576e-01,  /*39*/ -0.3390882631e-01,
    /*40*/ -0.3326609909e-01,  /*41*/ -0.3264299923e-01,
    /*42*/ -0.3203883174e-01,  /*43*/ -0.3145293280e-01,
    /*44*/ -0.3088463896e-01,  /*45*/ -0.3033334191e-01,
    /*46*/ -0.2979842108e-01,  /*47*/ -0.2927929373e-01,
    /*48*/ -0.2877539004e-01,  /*49*/ -0.2828615868e-01,
    /*50*/ -0.2781108192e-01,  /*51*/ -0.2734963353e-01,
    /*52*/ -0.2690133637e-01,  /*53*/ -0.2646570945e-01,
    /*54*/ -0.2604229418e-01,  /*55*/ -0.2563066630e-01,
    /*56*/ -0.2523039452e-01,  /*57*/ -0.2484107774e-01,
    /*58*/ -0.2446233127e-01,  /*59*/ -0.2409377795e-01,
    /*60*/ -0.2373506405e-01,  /*61*/ -0.2338584530e-01,
    /*62*/ -0.2304579309e-01,  /*63*/ -0.2271459030e-01,
    /*64*/ -0.2239193633e-01,  /*65*/ -0.2207753919e-01,
    /*66*/ -0.2177111937e-01,  /*67*/ -0.2147240908e-01,
    /*68*/ -0.2118115170e-01,  /*69*/ -0.2089710139e-01,
    /*70*/ -0.2062002059e-01,  /*71*/ -0.2034968220e-01,
    /*72*/ -0.2008586807e-01,  /*73*/ -0.1982836895e-01,
    /*74*/ -0.1957698319e-01,  /*75*/ -0.1933151728e-01,
    /*76*/ -0.1909178544e-01,  /*77*/ -0.1885760922e-01,
    /*78*/ -0.1862881706e-01,  /*79*/ -0.1840524275e-01,
    /*80*/ -0.1818672781e-01,  /*81*/ -0.1797312161e-01,
    /*82*/ -0.1776427277e-01,  /*83*/ -0.1756004164e-01,
    /*84*/ -0.1736029175e-01,  /*85*/ -0.1716489163e-01,
    /*86*/ -0.1697371498e-01,  /*87*/ -0.1678663942e-01,
    /*88*/ -0.1660354744e-01,  /*89*/ -0.1642432560e-01,
    /*90*/ -0.1624886478e-01,  /*91*/ -0.1607705911e-01,
    /*92*/ -0.1590880688e-01,  /*93*/ -0.1574400976e-01,
    /*94*/ -0.1558257306e-01,  /*95*/ -0.1542440486e-01,
    /*96*/ -0.1526941673e-01,  /*97*/ -0.1511752309e-01,
    /*98*/ -0.1496864146e-01,  /*99*/ -0.1482269168e-01,
    /*100*/ -0.1467959656e-01, /*101*/ -0.1453928135e-01,
    /*102*/ -0.1440167387e-01, /*103*/ -0.1426670401e-01,
    /*104*/ -0.1413430411e-01, /*105*/ -0.1400440861e-01,
    /*106*/ -0.1387695416e-01, /*107*/ -0.1375187918e-01,
    /*108*/ -0.1362912416e-01, /*109*/ -0.1350863142e-01,
    /*110*/ -0.1339034516e-01, /*111*/ -0.1327421107e-01,
    /*112*/ -0.1316017669e-01, /*113*/ -0.1304819106e-01,
    /*114*/ -0.1293820488e-01, /*115*/ -0.1283017014e-01,
    /*116*/ -0.1272404038e-01, /*117*/ -0.1261977050e-01,
    /*118*/ -0.1251731678e-01, /*119*/ -0.1241663664e-01,
    /*120*/ -0.1231768886e-01, /*121*/ -0.1222043336e-01,
    /*122*/ -0.1212483127e-01, /*123*/ -0.1203084476e-01,
    /*124*/ -0.1193843710e-01, /*125*/ -0.1184757260e-01,
    /*126*/ -0.1175821983e-01, /*127*/ -0.1167033892e-01,
    /*128*/ -0.1158389997e-01, /*129*/ -0.1149887118e-01,
    /*130*/ -0.1141522155e-01, /*131*/ -0.1133292103e-01,
    /*132*/ -0.1125194019e-01, /*133*/ -0.1117225053e-01,
    /*134*/ -0.1109382428e-01, /*135*/ -0.1101663454e-01,
    /*136*/ -0.1094065488e-01, /*137*/ -0.1086585976e-01,
    /*138*/ -0.1079222424e-01, /*139*/ -0.1071972412e-01,
    /*140*/ -0.1064833570e-01, /*141*/ -0.1057803597e-01,
    /*142*/ -0.1050880250e-01, /*143*/ -0.1044061352e-01,
    /*144*/ -0.1037344765e-01, /*145*/ -0.1030728418e-01,
    /*146*/ -0.1024210288e-01, /*147*/ -0.1017788410e-01,
    /*148*/ -0.1011460854e-01, /*149*/ -0.1005225749e-01,
    /*150*/ -0.9990812687e-02, /*151*/ -0.9930256346e-02,
    /*152*/ -0.9870571027e-02, /*153*/ -0.9811739803e-02,
    /*154*/ -0.9753746115e-02, /*155*/ -0.9696573861e-02,
    /*156*/ -0.9640207245e-02, /*157*/ -0.9584630905e-02,
    /*158*/ -0.9529829829e-02, /*159*/ -0.9475790189e-02,
    /*160*/ -0.9422496073e-02, /*161*/ -0.9369934305e-02,
    /*162*/ -0.9318091242e-02, /*163*/ -0.9266953580e-02,
    /*164*/ -0.9216508278e-02, /*165*/ -0.9166742633e-02,
    /*166*/ -0.9117644217e-02, /*167*/ -0.9069200890e-02,
    /*168*/ -0.9021400853e-02, /*169*/ -0.8974232435e-02,
    /*170*/ -0.8927684326e-02, /*171*/ -0.8881745447e-02,
    /*172*/ -0.8836405029e-02, /*173*/ -0.8791653407e-02,
    /*174*/ -0.8747481667e-02, /*175*/ -0.8703874181e-02,
    /*176*/ -0.8660824198e-02, /*177*/ -0.8618321971e-02,
    /*178*/ -0.8576358051e-02, /*179*/ -0.8534923154e-02,
    /*180*/ -0.8494008256e-02, /*181*/ -0.8453604422e-02,
    /*182*/ -0.8413702992e-02, /*183*/ -0.8374295452e-02,
    /*184*/ -0.8335373508e-02, /*185*/ -0.8296928978e-02,
    /*186*/ -0.8258953895e-02, /*187*/ -0.8221440446e-02,
    /*188*/ -0.8184381010e-02, /*189*/ -0.8147768058e-02,
    /*190*/ -0.8111594272e-02, /*191*/ -0.8075852456e-02,
    /*192*/ -0.8040535598e-02, /*193*/ -0.8005636771e-02,
    /*194*/ -0.7971149225e-02, /*195*/ -0.7937066335e-02,
    /*196*/ -0.7903381636e-02, /*197*/ -0.7870088750e-02,
    /*198*/ -0.7837181442e-02, /*199*/ -0.7804653597e-02,
    /*200*/ -0.7772499254e-02, /*201*/ -0.7740712516e-02,
    /*202*/ -0.7709287634e-02, /*203*/ -0.7678218959e-02,
    /*204*/ -0.7647500968e-02, /*205*/ -0.7617128214e-02,
    /*206*/ -0.7587095379e-02, /*207*/ -0.7557397242e-02,
    /*208*/ -0.7528028665e-02, /*209*/ -0.7498984668e-02,
    /*210*/ -0.7470260263e-02, /*211*/ -0.7441850636e-02,
    /*212*/ -0.7413751027e-02, /*213*/ -0.7385956807e-02,
    /*214*/ -0.7358463369e-02, /*215*/ -0.7331266234e-02,
    /*216*/ -0.7304360991e-02, /*217*/ -0.7277743343e-02,
    /*218*/ -0.7251409005e-02, /*219*/ -0.7225353835e-02,
    /*220*/ -0.7199573711e-02, /*221*/ -0.7174064656e-02,
    /*222*/ -0.7148822685e-02, /*223*/ -0.7123843938e-02,
    /*224*/ -0.7099129566e-02, /*225*/ -0.7074666551e-02,
    /*226*/ -0.7050455306e-02, /*227*/ -0.7026492490e-02,
    /*228*/ -0.7002774540e-02, /*229*/ -0.6979298005e-02,
    /*230*/ -0.6956059431e-02, /*231*/ -0.6933055469e-02,
    /*232*/ -0.6910282814e-02, /*233*/ -0.6887738228e-02,
    /*234*/ -0.6865418504e-02, /*235*/ -0.6843320529e-02,
    /*236*/ -0.6821441187e-02, /*237*/ -0.6799777491e-02,
    /*238*/ -0.6778326425e-02, /*239*/ -0.6757085071e-02,
    /*240*/ -0.6736050549e-02, /*241*/ -0.6715220049e-02,
    /*242*/ -0.6694590762e-02, /*243*/ -0.6674159952e-02,
    /*244*/ -0.6653925199e-02, /*245*/ -0.6633883871e-02,
    /*246*/ -0.6614032563e-02, /*247*/ -0.6594369232e-02,
    /*248*/ -0.6574891356e-02, /*249*/ -0.6555596461e-02,
    /*250*/ -0.6536482218e-02, /*251*/ -0.6517546015e-02,
    /*252*/ -0.6498785614e-02, /*253*/ -0.6480198683e-02,
    /*254*/ -0.6461782835e-02, /*255*/ -0.6443535971e-02,
    /*256*/ -0.6425455846e-02, /*257*/ -0.6407540263e-02,
    /*258*/ -0.6389786969e-02, /*259*/ -0.6372194033e-02,
    /*260*/ -0.6354759310e-02, /*261*/ -0.6337480778e-02,
    /*262*/ -0.6320356275e-02, /*263*/ -0.6303384010e-02,
    /*264*/ -0.6286561980e-02, /*265*/ -0.6269888285e-02,
    /*266*/ -0.6253360871e-02, /*267*/ -0.6236978078e-02,
    /*268*/ -0.6220737994e-02, /*269*/ -0.6204638851e-02,
    /*270*/ -0.6188678733e-02, /*271*/ -0.6172856062e-02,
    /*272*/ -0.6157169047e-02, /*273*/ -0.6141616004e-02,
    /*274*/ -0.6126195168e-02, /*275*/ -0.6110905001e-02,
    /*276*/ -0.6095743854e-02, /*277*/ -0.6080714363e-02,
    /*278*/ -0.6065808661e-02, /*279*/ -0.6051025499e-02,
    /*280*/ -0.6036365150e-02, /*281*/ -0.6021826113e-02,
    /*282*/ -0.6007406839e-02, /*283*/ -0.5993105878e-02,
    /*284*/ -0.5978921847e-02, /*285*/ -0.5964853301e-02,
    /*286*/ -0.5950898825e-02, /*287*/ -0.5937057006e-02,
    /*288*/ -0.5923326579e-02, /*289*/ -0.5909706151e-02,
    /*290*/ -0.5896194416e-02, /*291*/ -0.5882789995e-02,
    /*292*/ -0.5869491718e-02, /*293*/ -0.5856298304e-02,
    /*294*/ -0.5843208448e-02, /*295*/ -0.5830220904e-02,
    /*296*/ -0.5817334546e-02, /*297*/ -0.5804548161e-02,
    /*298*/ -0.5791860577e-02, /*299*/ -0.5779270522e-02,
    /*300*/ -0.5766776986e-02, /*301*/ -0.5754378816e-02,
    /*302*/ -0.5742074899e-02, /*303*/ -0.5729864025e-02,
    /*304*/ -0.5717745257e-02, /*305*/ -0.5705717481e-02,
    /*306*/ -0.5693779662e-02, /*307*/ -0.5681930652e-02,
    /*308*/ -0.5670169557e-02, /*309*/ -0.5658495321e-02,
    /*310*/ -0.5646906981e-02, /*311*/ -0.5635403426e-02,
    /*312*/ -0.5623983856e-02, /*313*/ -0.5612647214e-02,
    /*314*/ -0.5601392619e-02, /*315*/ -0.5590219039e-02,
    /*316*/ -0.5579125640e-02, /*317*/ -0.5568111507e-02,
    /*318*/ -0.5557175755e-02, /*319*/ -0.5546317422e-02,
    /*320*/ -0.5535535712e-02, /*321*/ -0.5524829768e-02,
    /*322*/ -0.5514198739e-02, /*323*/ -0.5503641734e-02,
    /*324*/ -0.5493157962e-02, /*325*/ -0.5482746632e-02,
    /*326*/ -0.5472406939e-02, /*327*/ -0.5462138054e-02,
    /*328*/ -0.5451939170e-02, /*329*/ -0.5441809594e-02,
    /*330*/ -0.5431748532e-02, /*331*/ -0.5421755212e-02,
    /*332*/ -0.5411828857e-02, /*333*/ -0.5401975569e-02,
    /*334*/ -0.5392181719e-02, /*335*/ -0.5382453176e-02,
    /*336*/ -0.5372788210e-02, /*337*/ -0.5363186748e-02,
    /*338*/ -0.5353648059e-02, /*339*/ -0.5344171504e-02,
    /*340*/ -0.5334756310e-02, /*341*/ -0.5325401926e-02,
    /*342*/ -0.5316107680e-02, /*343*/ -0.5306872920e-02,
    /*344*/ -0.5297696936e-02, /*345*/ -0.5288579215e-02,
    /*346*/ -0.5279519432e-02, /*347*/ -0.5270516389e-02,
    /*348*/ -0.5261569543e-02, /*349*/ -0.5252678535e-02,
    /*350*/ -0.5243842709e-02, /*351*/ -0.5235061516e-02,
    /*352*/ -0.5226334273e-02, /*353*/ -0.5217660524e-02,
    /*354*/ -0.5209039665e-02, /*355*/ -0.5200471169e-02,
    /*356*/ -0.5191954431e-02, /*357*/ -0.5183488963e-02,
    /*358*/ -0.5175074221e-02, /*359*/ -0.5166709698e-02,
    /*360*/ -0.5158394787e-02};
} /* unnamed namespace */

int dso::OceanPoleTide::stokes_coeffs(const dso::MjdEpoch &t, double xp,
                                      double yp, int max_degree, int max_order,
                                      double OmegaEarth, double G,
                                      double ge) noexcept {
  if (max_degree < 0) max_degree = this->max_degree();
  if (max_order  < 0) max_order = this->max_order();

  if ((max_degree < max_order) || (max_degree > this->max_degree()) ||
      (max_order > this->max_order())) {
    fprintf(stderr,
            "[ERROR] Invalid degree/order for Ocean Pole Tide expansion "
            "(traceback: %s)\n",
            __func__);
    return 1;
  }

  const double Re = mcs.Re();
  const double GM = mcs.GM();

  /* density of sea water in [kgm^−3] */
  constexpr const double rhow = 1025e0;
  constexpr const double g2_real = 0.6870e0;
  constexpr const double g2_imag = 0.0036e0;
  
  const auto m12 = dso::pole_tide_details::mcoeffs(t, xp, yp);
  const double m1 = dso::sec2rad(m12.m1); /* m1 in [rad] */
  const double m2 = dso::sec2rad(m12.m2); /* m2 in [rad] */
  
  const double freal = m1 * g2_real + m2 * g2_imag;
  const double fimag = m2 * g2_real - m1 * g2_imag;

  double knfac_arr[pole_tide_details::MAX_DEGREE_DESAI_2002];
  double *__restrict__ knfac = &(knfac_arr[0]);

  /* order: m = 0; compute kn factors */
  int m = 0;
  for (int l = m; l <= max_degree; l++) {
    knfac[l] = (1e0 + kn[l]) / (2 * l + 1);
    mcs.C(l, m) =
        knfac[l] * (A_real(l, m) * freal + A_imag(l, m) * fimag);
    mcs.S(l, m) =
        knfac[l] * (B_real(l, m) * freal + B_imag(l, m) * fimag);
  }
  
  /* all other orders (m) except m=0 */
  for (m = 1; m <= max_order; m++) {
    for (int l = m; l <= max_degree; l++) {
      mcs.C(l, m) =
          knfac[l] * (A_real(l, m) * freal + A_imag(l, m) * fimag);
      mcs.S(l, m) =
          knfac[l] * (B_real(l, m) * freal + B_imag(l, m) * fimag);
    //if (l==2) {
    //  printf("C(%d,%d) += %.13e * %.15f + %.13e * %.15f\n", l,m, A_real(l, m), knfac[l]*freal*fac2, A_imag(l, m), knfac[l]*fimag*fac2);
    //}
    }
  }
  
  /* scale (i.e. multiply) with constant term */
  const double fac = ((OmegaEarth * OmegaEarth) * std::pow(Re, 4) / GM) *
                     (2e0 * D2PI * G * rhow) / ge;
  mcs.scale(fac);

  /* all done */
  return 0;
}
