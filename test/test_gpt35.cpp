#include "tropo.hpp"
#include <cassert>
#include <cmath>
#include <cstdio>

/// --------------------------------------------------------------------------
/// Note
/// This test is pretty much automatically created using the matlab
/// script make_tests_35, which should be located in the var/gpt3/matlab
/// dir.
/// For instructions see the script.
/// It checks that the results for gpt3 (using the gpt3_5.grd grid file) are
/// the same to within PRECISION, when using the implementations from iers2010
/// C++ library and the original matlab source code.
/// --------------------------------------------------------------------------

constexpr const double PRECISION = 1e-12;

int main(int argc, char *argv[]) {
  if (argc != 2) {
    fprintf(stderr, "Usage: %s [gpt3_5 gride file]\n", argv[0]);
    return 1;
  }

  // Paste the output of make_tests_35.m withi the '//' block
  ////////////////////////////////////////////////////////////////////////////
  constexpr const int num_tests = 5;
dso::modified_julian_day mjd{45203};
dso::nanoseconds nanosec{60602224292792L};
dso::datetime<dso::nanoseconds> t{mjd, nanosec};
const double lats[] = {
	-0.11852614892908808208,
	+0.94744897313100162606,
	+1.44588545693157000116,
	+0.48890974171220147326,
-1.47743444412358981843};
const double lons[] = {
	-0.11570188376856194168,
	-2.42640168543301015802,
	+0.95865417532489871633,
	-0.85305715290699213327,
-0.07596389295925254359};
const double hgts[] = {
	+165.05962591308667697376,
	+1688.63784334831166233926,
	+1579.10162106680218130350,
	+1893.77759588068261109584,
+2155.55015299565820896532};
const dso::gpt3_result res35it0[] = {
{998.74938292689625995990,21.00954748668183569293,-9.43139015250267576107,287.33136543967907527986,18.74926034507060634837,0.00127524291389034838,0.00053813568201086567,3.66780429501609095411,16.26339260920092755214,0.00000399885611230212,-0.00000644983562123541,0.00011810865996845097,0.00003166097572714890},
{830.52640709515105754690,1.76780723696219332908,-6.51500760121828292881,275.71023425234875503520,6.09828194820626823258,0.00124287440319618167,0.00058407474012442363,2.70695200223316190247,-2.22730711328744357402,-0.00019255414166097353,0.00004985504016526740,-0.00004437414844702558,0.00001425774641950164},
{830.38528644674283896165,-5.15664975945731995921,-3.06538915937093081610,266.04069798833376125913,3.10492512420901034531,0.00121852255670867101,0.00056562099123538844,2.13999392925577547331,16.78332386699947775810,-0.00007795228150391820,-0.00009230657276061030,-0.00002310376852034346,0.00000334485788807370},
{822.73248080788584957190,8.59322509879115337128,-9.37577661088478997442,287.83914393418058352836,10.33851067942581281045,0.00127368696174704659,0.00058648050529708793,3.05745023078738675792,-15.86511564589401146463,0.00004541023402874210,0.00008590337507907085,-0.00000076669013304750,-0.00006055457239893031},
{741.42698013871699913580,-39.11432675978882400614,-1.64943909262694043960,230.81821183613806169888,0.12465026639509274586,0.00109991504769876326,0.00047943194332183807,1.45224758746786797303,-15.98708459315884766738,0.00005544903104689129,0.00000971691551718601,0.00000029631944738992,-0.00000184613295355759}};
const dso::gpt3_result res35it1[] = {
{996.32268727435553046234,22.90970532335541420821,-9.44575318491824944545,287.98780125546676345039,21.53561364445598513839,0.00127627671008407554,0.00056463697667319423,3.29706762810369458094,16.26339260920092755214,0.00000838889007982624,-0.00000473921305863234,0.00020780917893557950,0.00004997575891553748},
{823.06848273860578046879,-4.58622698134561179728,-7.16839823189727365360,270.11283147998796039246,3.91656363290034459368,0.00122307201760208360,0.00055728155557138964,3.02386995638889777283,-2.22730711328744357402,-0.00015785375175686654,0.00009235387168294147,-0.00003126351096877020,0.00001721820976937581},
{825.89497456834919830726,-13.69008320419446178562,-1.70098208191635702491,256.59645440655560832965,1.64070733029487758614,0.00118517278670451486,0.00056470035773775935,2.06755697967618345956,16.78332386699947775810,-0.00003420496309467471,-0.00009085615602314772,-0.00001243607027907334,-0.00000458051427868302},
{819.87724184941339444777,5.43100625847463014395,-9.26925211422357975266,285.00947450667382554457,8.00638409328598221748,0.00126336128127362154,0.00056752747551611625,3.29804198811871707164,-15.86511564589401146463,-0.00006171462117153775,0.00007950904937071968,-0.00001228307304106621,-0.00002170900095918119},
{745.96031733896165860642,-31.56693830957890156697,-2.65342296174124125940,235.73459031599270474544,0.38080625165173914715,0.00113999624113587201,0.00047396785118855665,2.12074644748551977003,-15.98708459315884766738,0.00004753233644294096,0.00004761953949219343,0.00000079943968842947,-0.00000107402732629839}};
const double zd[] = {
	+0.65350776710112723489,
	+0.92780590228259396390,
	+0.10414005723476622844,
	+0.17086827920614433562,
+0.74429259469149977591};
const double resvmf3[][2] = {
	{1.25858097602313789842,1.25910033664729370528},
	{1.66421455583841337855,1.66615699986320220383},
	{1.00543425340819747582,1.00544100489251220232},
	{1.01473977851154395502,1.01476061120904637036},
	{1.35819128987430248756,1.35895080701129433720},
};
  ////////////////////////////////////////////////////////////////////////////

  /* combine lon,lat,hgt into arrays[3] */
  std::vector<std::array<double, 3>> input;
  for (int i = 0; i < num_tests; i++)
    input.emplace_back(std::array<double, 3>{lons[i], lats[i], hgts[i]});

  /* read in the grid */
  dso::Gpt3Grid grid(argv[1]);
  std::vector<dso::gpt3_result> gout;

  /* compute gpt3 using it=0 */
  if (dso::gpt3_fast(t, input, 0, grid, gout)) {
    fprintf(stderr, "ERROR. Failed call to gpt3_fast\n");
    return 9;
  }
  /* check results vs matlab */
  for (int i = 0; i < num_tests; i++) {
    assert(std::abs(gout[i].p - res35it0[i].p) < PRECISION);
    assert(std::abs(gout[i].T - res35it0[i].T) < PRECISION);
    assert(std::abs(gout[i].dT - res35it0[i].dT) < PRECISION);
    assert(std::abs(gout[i].Tm - res35it0[i].Tm) < PRECISION);
    assert(std::abs(gout[i].e - res35it0[i].e) < PRECISION);
    assert(std::abs(gout[i].ah - res35it0[i].ah) < PRECISION);
    assert(std::abs(gout[i].aw - res35it0[i].aw) < PRECISION);
    assert(std::abs(gout[i].la - res35it0[i].la) < PRECISION);
    assert(std::abs(gout[i].undu - res35it0[i].undu) < PRECISION);
    assert(std::abs(gout[i].Gn_h - res35it0[i].Gn_h) < PRECISION);
    assert(std::abs(gout[i].Ge_h - res35it0[i].Ge_h) < PRECISION);
    assert(std::abs(gout[i].Gn_w - res35it0[i].Gn_w) < PRECISION);
    assert(std::abs(gout[i].Ge_w - res35it0[i].Ge_w) < PRECISION);
  }

  /* compute gpt3 using it=1 */
  if (dso::gpt3_fast(t, input, 1, grid, gout)) {
    fprintf(stderr, "ERROR. Failed call to gpt3_fast\n");
    return 9;
  }
  /* check results vs matlab */
  for (int i = 0; i < num_tests; i++) {
    assert(std::abs(gout[i].p - res35it1[i].p) < PRECISION);
    assert(std::abs(gout[i].T - res35it1[i].T) < PRECISION);
    assert(std::abs(gout[i].dT - res35it1[i].dT) < PRECISION);
    assert(std::abs(gout[i].Tm - res35it1[i].Tm) < PRECISION);
    assert(std::abs(gout[i].e - res35it1[i].e) < PRECISION);
    assert(std::abs(gout[i].ah - res35it1[i].ah) < PRECISION);
    assert(std::abs(gout[i].aw - res35it1[i].aw) < PRECISION);
    assert(std::abs(gout[i].la - res35it1[i].la) < PRECISION);
    assert(std::abs(gout[i].undu - res35it1[i].undu) < PRECISION);
    assert(std::abs(gout[i].Gn_h - res35it1[i].Gn_h) < PRECISION);
    assert(std::abs(gout[i].Ge_h - res35it1[i].Ge_h) < PRECISION);
    assert(std::abs(gout[i].Gn_w - res35it1[i].Gn_w) < PRECISION);
    assert(std::abs(gout[i].Ge_w - res35it1[i].Ge_w) < PRECISION);
  }

  /* check vmf3 results */
  double mfh, mfw;
  for (int i = 0; i < num_tests; i++) {
    if (dso::vmf3(gout[i].ah, gout[i].aw, t, input[i][1], input[i][2], zd[i],
                  mfh, mfw)) {
      fprintf(stderr, "ERROR. Failed call to vmf3\n");
      return 10;
    }
    printf("\tExpecting: mfh=%+.12f mfw=%+.12f\n", resvmf3[i][0], resvmf3[i][1]);
    printf("\tGot      : mfh=%+.12f mfw=%+.12f\n", mfh, mfw);
    assert(std::abs(mfh - resvmf3[i][0]) < PRECISION);
    assert(std::abs(mfw - resvmf3[i][1]) < PRECISION);
  }

  return 0;
}
