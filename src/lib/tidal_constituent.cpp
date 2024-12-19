#include "tide_atlas.hpp"

dso::TidalConstituent::TidalConstituent(const dso::detail::TidalConstituentArrayEntry *wave,
                      int max_degree, int max_order) noexcept
    : mwave(*wave), mCosCs(max_degree, max_order),
      mSinCs(max_degree, max_order){};

dso::TidalConstituent::TidalConstituent(const dso::TidalWave &wave, int max_degree,
                      int max_order) noexcept
    : mwave(wave), mCosCs(max_degree, max_order),
      mSinCs(max_degree, max_order){};
