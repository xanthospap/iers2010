#include "eop.hpp"
#include <cassert>
#include <cstdio>

int main(int argc, char *argv[]) {
  if (argc != 2) {
    fprintf(stderr, "Usage %s [EOP C04 FILE]\n", argv[0]);
    return 1;
  }

  dso::MjdEpoch t1(57750);
  dso::MjdEpoch t2(57759);
  dso::EopSeries eop;

  if (dso::parse_iers_C04(argv[1], t1, t2, eop)) {
    fprintf(stderr, "ERROR Failed parsing eop file\n");
    return 1;
  }
  assert(eop.num_entries() == 9);

  /* make a copy of the original series */
  dso::EopSeries eop_copy(eop);

  /* regularize the series */
  eop_copy.regularize();
  for (int i = 0; i < eop_copy.num_entries(); i++) {
    assert(eop_copy.give_me_the_vector()[i].t() ==
           eop.give_me_the_vector()[i].t());
    assert(eop_copy.give_me_the_vector()[i].xp() ==
           eop.give_me_the_vector()[i].xp());
    assert(eop_copy.give_me_the_vector()[i].yp() ==
           eop.give_me_the_vector()[i].yp());
    assert(eop_copy.give_me_the_vector()[i].dX() ==
           eop.give_me_the_vector()[i].dX());
    assert(eop_copy.give_me_the_vector()[i].dY() ==
           eop.give_me_the_vector()[i].dY());
    assert(eop_copy.give_me_the_vector()[i].dut() !=
                    eop.give_me_the_vector()[i].dut());
    assert(eop_copy.give_me_the_vector()[i].lod() !=
                    eop.give_me_the_vector()[i].lod());
  }

  /* de-regularize series */
  eop_copy.deregularize();

  /* make sure we are back at the original series */
  for (int i = 0; i < eop_copy.num_entries(); i++) {
    assert(eop_copy.give_me_the_vector()[i].t() ==
           eop.give_me_the_vector()[i].t());
    assert(eop_copy.give_me_the_vector()[i].xp() ==
           eop.give_me_the_vector()[i].xp());
    assert(eop_copy.give_me_the_vector()[i].yp() ==
           eop.give_me_the_vector()[i].yp());
    assert(eop_copy.give_me_the_vector()[i].dX() ==
           eop.give_me_the_vector()[i].dX());
    assert(eop_copy.give_me_the_vector()[i].dY() ==
           eop.give_me_the_vector()[i].dY());
    //printf("Original %.12f %.12f\n", eop.give_me_the_vector()[i].dut(),
    //       eop.give_me_the_vector()[i].lod());
    //printf("         %.12f %.12f\n", eop_copy.give_me_the_vector()[i].dut(),
    //       eop_copy.give_me_the_vector()[i].lod());
    //printf("         %.1e %.1e\n",
    //       std::abs(eop_copy.give_me_the_vector()[i].dut() -
    //                eop.give_me_the_vector()[i].dut()),
    //       std::abs(eop_copy.give_me_the_vector()[i].lod() -
    //                eop.give_me_the_vector()[i].lod()));
    assert(std::abs(eop_copy.give_me_the_vector()[i].dut() -
                    eop.give_me_the_vector()[i].dut()) < 1e-15);
    assert(std::abs(eop_copy.give_me_the_vector()[i].lod() -
                    eop.give_me_the_vector()[i].lod()) < 1e-15);
  }

  return 0;
}
