#include "iers2010.hpp"
#include <cassert>

int main() {
  const double PRECISION = 1e-11;
  double fundarg[5];
  // testing fundarg
     const double ijc[] = {
.000000000000000E+00,
.100000000000000E+00,
.100000000000000E-01,
.100000000000000E-02,
.100000000000000E-03,
.100000000000000E-04,
.100000000000000E-05,
.100000000000000E-06,
.100000000000000E-07,
.100000000000000E-08,
.100000000000000E-09,
.100000000000000E-10,
.100000000000000E-11,
.100000000000000E-12,
.100000000000000E-13,
.100000000000000E-14,
.123456789000100E+00,
.112345678900010E+01,
.212345678900010E+01,
.412345678900010E+01,
.512345678900010E+01,
.199999999900010E+01,
.199999999900000E+01,
.199999999999990E+01   };
iers2010::fundarg(ijc[  0], fundarg);
  assert(std::abs(fundarg[0]-( 0.235555574349388E+01 ))<PRECISION);
  assert(std::abs(fundarg[1]-( 0.624006012692298E+01 ))<PRECISION);
  assert(std::abs(fundarg[2]-( 0.162790508153752E+01 ))<PRECISION);
  assert(std::abs(fundarg[3]-( 0.519846658865050E+01 ))<PRECISION);
  assert(std::abs(fundarg[4]-( 0.218243919661567E+01 ))<PRECISION);
iers2010::fundarg(ijc[  1], fundarg);
  assert(std::abs(fundarg[0]-( 0.584423931349453E+01 ))<PRECISION);
  assert(std::abs(fundarg[1]-( 0.623840254544757E+01 ))<PRECISION);
  assert(std::abs(fundarg[2]-( 0.302768899290967E+01 ))<PRECISION);
  assert(std::abs(fundarg[3]-( 0.322120274892972E+01 ))<PRECISION);
  assert(std::abs(fundarg[4]-(-0.119326503644763E+01 ))<PRECISION);
iers2010::fundarg(ijc[  2], fundarg);
  assert(std::abs(fundarg[0]-( 0.396106102280584E+01 ))<PRECISION);
  assert(std::abs(fundarg[1]-( 0.623989437118916E+01 ))<PRECISION);
  assert(std::abs(fundarg[2]-( 0.428115765118467E+01 ))<PRECISION);
  assert(std::abs(fundarg[3]-( 0.123082904816449E+01 ))<PRECISION);
  assert(std::abs(fundarg[4]-( 0.184486874070202E+01 ))<PRECISION);
iers2010::fundarg(ijc[  3], fundarg);
  assert(std::abs(fundarg[0]-( 0.440106186218793E+01 ))<PRECISION);
  assert(std::abs(fundarg[1]-( 0.585176774912109E+00 ))<PRECISION);
  assert(std::abs(fundarg[2]-( 0.377818593121249E+01 ))<PRECISION);
  assert(std::abs(fundarg[3]-( 0.403473119854158E+00 ))<PRECISION);
  assert(std::abs(fundarg[4]-( 0.214868215069827E+01 ))<PRECISION);
iers2010::fundarg(ijc[  4], fundarg);
  assert(std::abs(fundarg[0]-( 0.318842488606733E+01 ))<PRECISION);
  assert(std::abs(fundarg[1]-( 0.197050152605070E-01 ))<PRECISION);
  assert(std::abs(fundarg[2]-( 0.247125169722854E+01 ))<PRECISION);
  assert(std::abs(fundarg[3]-( 0.597560430320957E+01 ))<PRECISION);
  assert(std::abs(fundarg[4]-( 0.217906349202067E+01 ))<PRECISION);
iers2010::fundarg(ijc[  5], fundarg);
  assert(std::abs(fundarg[0]-( 0.243884265775109E+01 ))<PRECISION);
  assert(std::abs(fundarg[1]-( 0.624634314647469E+01 ))<PRECISION);
  assert(std::abs(fundarg[2]-( 0.171223974310668E+01 ))<PRECISION);
  assert(std::abs(fundarg[3]-( 0.527618036010644E+01 ))<PRECISION);
  assert(std::abs(fundarg[4]-( 0.218210162615614E+01 ))<PRECISION);
iers2010::fundarg(ijc[  6], fundarg);
  assert(std::abs(fundarg[0]-( 0.236388443491960E+01 ))<PRECISION);
  assert(std::abs(fundarg[1]-( 0.624068842887815E+01 ))<PRECISION);
  assert(std::abs(fundarg[2]-( 0.163633854769444E+01 ))<PRECISION);
  assert(std::abs(fundarg[3]-( 0.520623796579610E+01 ))<PRECISION);
  assert(std::abs(fundarg[4]-( 0.218240543956972E+01 ))<PRECISION);
iers2010::fundarg(ijc[  7], fundarg);
  assert(std::abs(fundarg[0]-( 0.235638861263645E+01 ))<PRECISION);
  assert(std::abs(fundarg[1]-( 0.624012295711850E+01 ))<PRECISION);
  assert(std::abs(fundarg[2]-( 0.162874842815321E+01 ))<PRECISION);
  assert(std::abs(fundarg[3]-( 0.519924372636506E+01 ))<PRECISION);
  assert(std::abs(fundarg[4]-( 0.218243582091108E+01 ))<PRECISION);
iers2010::fundarg(ijc[  8], fundarg);
  assert(std::abs(fundarg[0]-( 0.235563903040814E+01 ))<PRECISION);
  assert(std::abs(fundarg[1]-( 0.624006640994253E+01 ))<PRECISION);
  assert(std::abs(fundarg[2]-( 0.162798941619909E+01 ))<PRECISION);
  assert(std::abs(fundarg[3]-( 0.519854430242196E+01 ))<PRECISION);
  assert(std::abs(fundarg[4]-( 0.218243885904521E+01 ))<PRECISION);
iers2010::fundarg(ijc[  9], fundarg);
  assert(std::abs(fundarg[0]-( 0.235556407218531E+01 ))<PRECISION);
  assert(std::abs(fundarg[1]-( 0.624006075522493E+01 ))<PRECISION);
  assert(std::abs(fundarg[2]-( 0.162791351500368E+01 ))<PRECISION);
  assert(std::abs(fundarg[3]-( 0.519847436002765E+01 ))<PRECISION);
  assert(std::abs(fundarg[4]-( 0.218243916285862E+01 ))<PRECISION);
iers2010::fundarg(ijc[ 10], fundarg);
  assert(std::abs(fundarg[0]-( 0.235555657636302E+01 ))<PRECISION);
  assert(std::abs(fundarg[1]-( 0.624006018975318E+01 ))<PRECISION);
  assert(std::abs(fundarg[2]-( 0.162790592488413E+01 ))<PRECISION);
  assert(std::abs(fundarg[3]-( 0.519846736578822E+01 ))<PRECISION);
  assert(std::abs(fundarg[4]-( 0.218243919323997E+01 ))<PRECISION);
iers2010::fundarg(ijc[ 11], fundarg);
  assert(std::abs(fundarg[0]-( 0.235555582678079E+01 ))<PRECISION);
  assert(std::abs(fundarg[1]-( 0.624006013320600E+01 ))<PRECISION);
  assert(std::abs(fundarg[2]-( 0.162790516587218E+01 ))<PRECISION);
  assert(std::abs(fundarg[3]-( 0.519846666636428E+01 ))<PRECISION);
  assert(std::abs(fundarg[4]-( 0.218243919627810E+01 ))<PRECISION);
iers2010::fundarg(ijc[ 12], fundarg);
  assert(std::abs(fundarg[0]-( 0.235555575182257E+01 ))<PRECISION);
  assert(std::abs(fundarg[1]-( 0.624006012755128E+01 ))<PRECISION);
  assert(std::abs(fundarg[2]-( 0.162790508997099E+01 ))<PRECISION);
  assert(std::abs(fundarg[3]-( 0.519846659642188E+01 ))<PRECISION);
  assert(std::abs(fundarg[4]-( 0.218243919658191E+01 ))<PRECISION);
iers2010::fundarg(ijc[ 13], fundarg);
  assert(std::abs(fundarg[0]-( 0.235555574432675E+01 ))<PRECISION);
  assert(std::abs(fundarg[1]-( 0.624006012698581E+01 ))<PRECISION);
  assert(std::abs(fundarg[2]-( 0.162790508238087E+01 ))<PRECISION);
  assert(std::abs(fundarg[3]-( 0.519846658942764E+01 ))<PRECISION);
  assert(std::abs(fundarg[4]-( 0.218243919661230E+01 ))<PRECISION);
iers2010::fundarg(ijc[ 14], fundarg);
  assert(std::abs(fundarg[0]-( 0.235555574357717E+01 ))<PRECISION);
  assert(std::abs(fundarg[1]-( 0.624006012692926E+01 ))<PRECISION);
  assert(std::abs(fundarg[2]-( 0.162790508162185E+01 ))<PRECISION);
  assert(std::abs(fundarg[3]-( 0.519846658872822E+01 ))<PRECISION);
  assert(std::abs(fundarg[4]-( 0.218243919661533E+01 ))<PRECISION);
iers2010::fundarg(ijc[ 15], fundarg);
  assert(std::abs(fundarg[0]-( 0.235555574350221E+01 ))<PRECISION);
  assert(std::abs(fundarg[1]-( 0.624006012692361E+01 ))<PRECISION);
  assert(std::abs(fundarg[2]-( 0.162790508154595E+01 ))<PRECISION);
  assert(std::abs(fundarg[3]-( 0.519846658865828E+01 ))<PRECISION);
  assert(std::abs(fundarg[4]-( 0.218243919661564E+01 ))<PRECISION);
iers2010::fundarg(ijc[ 16], fundarg);
  assert(std::abs(fundarg[0]-( 0.146667714113844E+00 ))<PRECISION);
  assert(std::abs(fundarg[1]-( 0.212679300065726E+01 ))<PRECISION);
  assert(std::abs(fundarg[2]-( 0.607098032857748E+01 ))<PRECISION);
  assert(std::abs(fundarg[3]-( 0.330038262325206E+01 ))<PRECISION);
  assert(std::abs(fundarg[4]-(-0.198509675073301E+01 ))<PRECISION);
iers2010::fundarg(ijc[ 17], fundarg);
  assert(std::abs(fundarg[0]-( 0.361775448923671E+01 ))<PRECISION);
  assert(std::abs(fundarg[1]-( 0.211021411073158E+01 ))<PRECISION);
  assert(std::abs(fundarg[2]-( 0.121919261202679E+01 ))<PRECISION);
  assert(std::abs(fundarg[3]-( 0.237726476923700E+01 ))<PRECISION);
  assert(std::abs(fundarg[4]-(-0.432617094504672E+01 ))<PRECISION);
iers2010::fundarg(ijc[ 18], fundarg);
  assert(std::abs(fundarg[0]-( 0.805966733716738E+00 ))<PRECISION);
  assert(std::abs(fundarg[1]-( 0.209362986031654E+01 ))<PRECISION);
  assert(std::abs(fundarg[2]-( 0.265046652998818E+01 ))<PRECISION);
  assert(std::abs(fundarg[3]-( 0.145408535696608E+01 ))<PRECISION);
  assert(std::abs(fundarg[4]-(-0.383987132920292E+00 ))<PRECISION);
iers2010::fundarg(ijc[ 19], fundarg);
  assert(std::abs(fundarg[0]-( 0.146651465411506E+01 ))<PRECISION);
  assert(std::abs(fundarg[1]-( 0.206044528382486E+01 ))<PRECISION);
  assert(std::abs(fundarg[2]-( 0.551264323088849E+01 ))<PRECISION);
  assert(std::abs(fundarg[3]-( 0.589072790433524E+01 ))<PRECISION);
  assert(std::abs(fundarg[4]-(-0.506577118086145E+01 ))<PRECISION);
iers2010::fundarg(ijc[ 20], fundarg);
  assert(std::abs(fundarg[0]-( 0.493885315617898E+01 ))<PRECISION);
  assert(std::abs(fundarg[1]-( 0.204384495730885E+01 ))<PRECISION);
  assert(std::abs(fundarg[2]-( 0.660360649354012E+00 ))<PRECISION);
  assert(std::abs(fundarg[3]-( 0.496736491732693E+01 ))<PRECISION);
  assert(std::abs(fundarg[4]-(-0.112336802165395E+01 ))<PRECISION);
iers2010::fundarg(ijc[ 21], fundarg);
  assert(std::abs(fundarg[0]-( 0.301476972330248E+01 ))<PRECISION);
  assert(std::abs(fundarg[1]-( 0.620689768198652E+01 ))<PRECISION);
  assert(std::abs(fundarg[2]-( 0.449059869350620E+01 ))<PRECISION);
  assert(std::abs(fundarg[3]-( 0.335217675427253E+01 ))<PRECISION);
  assert(std::abs(fundarg[4]-(-0.249965440598928E+01 ))<PRECISION);
iers2010::fundarg(ijc[ 22], fundarg);
  assert(std::abs(fundarg[0]-( 0.301476972246793E+01 ))<PRECISION);
  assert(std::abs(fundarg[1]-( 0.620689768192366E+01 ))<PRECISION);
  assert(std::abs(fundarg[2]-( 0.449059869266240E+01 ))<PRECISION);
  assert(std::abs(fundarg[3]-( 0.335217675349577E+01 ))<PRECISION);
  assert(std::abs(fundarg[4]-(-0.249965440598589E+01 ))<PRECISION);
iers2010::fundarg(ijc[ 23], fundarg);
  assert(std::abs(fundarg[0]-( 0.301477805032755E+01 ))<PRECISION);
  assert(std::abs(fundarg[1]-( 0.620689831016282E+01 ))<PRECISION);
  assert(std::abs(fundarg[2]-( 0.449060712528458E+01 ))<PRECISION);
  assert(std::abs(fundarg[3]-( 0.335218452409587E+01 ))<PRECISION);
  assert(std::abs(fundarg[4]-(-0.249965443973942E+01 ))<PRECISION);
  return 0;
}
