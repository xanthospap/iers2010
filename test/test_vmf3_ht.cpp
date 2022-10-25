#include "tropo.hpp"
#include <cstdio>
#include <cmath>
#include <cassert>

/// --------------------------------------------------------------------------
/// Note
/// This test is pretty much automatically created using the matlab 
/// script make_tests_31, which should be located in the var/gpt3/matlab
/// dir.
/// For instructions see the script.
/// It checks that the results for gpt3 (using the gpt3_1.grd grid file) are
/// the same to within PRECISION, when using the implementations from iers2010
/// C++ library and the original matlab source code.
//constexpr const int num_tests = 10;

constexpr const double PRECISION = 1e-12;

int main() {
  [[maybe_unused]] constexpr const int num_tests = 10;
const dso::datetime<dso::nanoseconds> ts[] = {
	dso::datetime<dso::nanoseconds>(dso::modified_julian_day(49910), dso::nanoseconds(63840313720400L)),
	dso::datetime<dso::nanoseconds>(dso::modified_julian_day(47252), dso::nanoseconds(85454786300822L)),
	dso::datetime<dso::nanoseconds>(dso::modified_julian_day(52573), dso::nanoseconds(45968822588632L)),
	dso::datetime<dso::nanoseconds>(dso::modified_julian_day(60309), dso::nanoseconds(74895554286777L)),
	dso::datetime<dso::nanoseconds>(dso::modified_julian_day(47176), dso::nanoseconds(14428517423314L)),
	dso::datetime<dso::nanoseconds>(dso::modified_julian_day(56552), dso::nanoseconds(37377501566196L)),
	dso::datetime<dso::nanoseconds>(dso::modified_julian_day(58950), dso::nanoseconds(21841862341063L)),
	dso::datetime<dso::nanoseconds>(dso::modified_julian_day(45303), dso::nanoseconds(22580717394152L)),
	dso::datetime<dso::nanoseconds>(dso::modified_julian_day(48817), dso::nanoseconds(3412187595642L)),
	dso::datetime<dso::nanoseconds>(dso::modified_julian_day(59043), dso::nanoseconds(40460378009873L))};
const double lats[] = {
	-1.00655839762322463216,
	-0.99961005735638530645,
	-0.72809675130317486058,
	-0.14342150316881907557,
	-1.43112307333969557277,
	-0.02487747751083047021,
	-0.71993925109076728486,
	+0.82448949290365325382,
	-0.58926989005100749441,
+1.27979393855173473327};
const double lons[] = {
	-1.26984600639700895286,
	-1.53555962473150975356,
	+0.64657972158758836301,
	+0.99219874199870705667,
	-1.37089139055061060013,
	-1.12132861360957836183,
	+2.00668485003957730584,
	-1.80861638303667837668,
	-1.43854519911185629866,
-1.58509209250454596862};
const double hgts[] = {
	+829.89198810515767945617,
	+2024.08432814885327388765,
	+1200.01280238938329603116,
	+1009.32952967655137399561,
	+2814.83230872419426304987,
	+1504.49462636883822597156,
	+1679.09578652563459399971,
	+83.82284862202868680470,
	+2010.43036092576539886068,
+1146.66644991213684079412};
const double zds[] = {
	+1.47974689287359240097,
	+1.00938291119670298102,
	+0.19754704707866257984,
	+0.20592471138694495858,
	+1.08709711214142656921,
	+0.52481859559835963847,
	+1.13420562474041242140,
	+0.33307940529340013258,
	+1.56713618809026944412,
+0.36094074669212516948};
const double ahs[] = {
	+0.00357300986953864620,
	+0.00331715125402642485,
	+0.01632044127182305515,
	+0.01885650717138715252,
	+0.01096636564494461036,
	+0.01387855180187660617,
	+0.00282556259640217352,
	+0.00631341240425861532,
	+0.01422581960880303295,
+0.00280761457543542235};
const double aws[] = {
	+0.00031829758849169666,
	+0.00027656306962914589,
	+0.00169622903295432812,
	+0.00140194053314290413,
	+0.00085074300404099242,
	+0.00148759202241975056,
	+0.00077474758675946903,
	+0.00151078488330494485,
	+0.00141396019100691289,
+0.00167577124632660316};

const double mfhr[] = {
	+8.12189908207207977853,
	+1.86318272595860490171,
	+1.01918370014659420697,
	+1.02076474149959128823,
	+2.07055200013497930911,
	+1.15028248968598201252,
	+2.33558828564583320642,
	+1.05736669532487703727,
	+6.66270689167481577897,
+1.06845328343232970347};
const double mfwr[] = {
	+10.63468585983741654388,
	+1.87704170912981260777,
	+1.01976581586607806251,
	+1.02152145342352418389,
	+2.14371381767069602020,
	+1.15494183279158391109,
	+2.35658184533182124554,
	+1.05796573151444350813,
	+23.01069694070993065793,
+1.06861877217342082780};
double _mfh,_mfw;
for (int i=0; i<10; i++) {
	dso::vmf3_ht(ahs[i], aws[i], ts[i], lats[i], lons[i], hgts[i], zds[i], _mfh, _mfw);
	printf("DeltaMfh=%+15.17e (aka %.12f %.12f) DeltaMfw=%+15.17e (aka %.12f %.12f)\n", std::abs(_mfh-mfhr[i]), _mfh, mfhr[i], std::abs(_mfw-mfwr[i]), _mfw, mfwr[i]);

}
return 0;
}
