#include "doodson.hpp"
#include <cassert>
#include <cstdio>
#include <string>
#include <vector>

struct Ref {
  std::string s;
  int ar[6];
};

const std::vector<Ref> dstr = {
    Ref{std::string{"155.555"}, {1, 0, 0, 0, 0, 0}},
    Ref{std::string{"255.555"}, {2, 0, 0, 0, 0, 0}},
    Ref{std::string{"256.789"}, {2, 0, 1, 2, 3, 4}},
    Ref{std::string{"098.765"}, {0, 4, 3, 2, 1, 0}},
    Ref{std::string{"0a8.765"}, {0, 5, 3, 2, 1, 0}}, /* #5 */
    Ref{std::string{"0A8.765"}, {0, 5, 3, 2, 1, 0}},
    Ref{std::string{"0m8.765"}, {0, 22 - 5, 3, 2, 1, 0}},
    Ref{std::string{"0m8.765           "}, {0, 22 - 5, 3, 2, 1, 0}},
    Ref{std::string{"098.7f5"}, {0, 4, 3, 2, 15 - 5, 0}},
    Ref{std::string{"0ab.cde"},
        {0, 10 - 5, 11 - 5, 12 - 5, 13 - 5, 14 - 5}}, /* #10 */
    Ref{std::string{"0no.pqr"},
        {0, -13 - 5, -12 - 5, -11 - 5, -10 - 5, -9 - 5}},
    Ref{std::string{"0no.pqr  "},
        {0, -13 - 5, -12 - 5, -11 - 5, -10 - 5, -9 - 5}},
    Ref{std::string{"0NO.PQR  "},
        {0, -13 - 5, -12 - 5, -11 - 5, -10 - 5, -9 - 5}}, /* #13 */
    /* erronuous */
    Ref{std::string{"0-10.555"}, {0, 0, 0, 0, 0, 0}},
    Ref{std::string{"010.55-5"}, {0, 0, 0, 0, 0, 0}},
    Ref{std::string{"010555"}, {0, 0, 0, 0, 0, 0}},
    Ref{std::string{"010,555"}, {0, 0, 0, 0, 0, 0}}};

int main() {
  int num = 0;
  [[maybe_unused]] char buf1[8] = "\0";
  [[maybe_unused]] char buf2[8] = "\0";
  for (const auto &r : dstr) {
    try {
      /* Doodson from string */
      dso::DoodsonConstituent d1 =
          dso::DoodsonConstituent::from_chars(r.s.c_str());
      
      /* Doodson from int array */
      dso::DoodsonConstituent d2(r.ar);

      /* the two should match exactly */
      assert(d1 == d2);

      /* assert resolved ints match the ones given */
      for (int i = 0; i < 6; i++)
        assert((d1.int_array()[i] == d2.int_array()[i]) &&
               (d1.int_array()[i] == r.ar[i]));
    } catch (std::exception &) {
      if (num < 13) {
        fprintf(stderr, "[ERROR] Failed converting Doodson wave: [%s]\n",
                r.s.c_str());
        return 5;
      }
    }
    ++num;
  }

  return 0;
}
