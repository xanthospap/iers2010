/** @file
 * Decleration/definition of a template Aod1bDataStream<T>, where T is 
 * any AOD1B coefficient type; this class enables easy, sequential parsing
 * of the data stored in AOD1B files, in chronological order.
 */

#ifndef __DSO_AOD1B_IN_DATA_STREAM_HPP__
#define __DSO_AOD1B_IN_DATA_STREAM_HPP__

#include "aod1b.hpp"

namespace dso {

template <AOD1BCoefficientType T> class Aod1bDataStream {
  static constexpr const int AR_SIZE = 3;
private:
  Aod1bDataBlockIterator<T> mit;
  std::array<Aod1bIn::Aod1bBlockHeader, AR_SIZE> mhdrvec;
  StokesCoeffs *mcsvec{nullptr};

  void swap(int first, int second) noexcept {
    using std::swap;
    swap(mhdrvec[first], mhdrvec[second]);
    swap(mcsvec[first], mcsvec[second]);
  }

  int hunt_range(const Datetime<nanoseconds> &t) noexcept {
    /* quick return */
    if (t >= mhdrvec[0].mepoch && t < mhdrvec[1].mepoch) return 0;
    if (t >= mhdrvec[1].mepoch && t < mhdrvec[2].mepoch) return 1;
   
    /* shit, we don;t have the bounding coeffs!
     * First, read in one more header/block, and see if we are there yet.
     * Else, keep on reading untill we find a header that is > t (storing
     * coeffs allong the way). 
     * Both are only possible if the iterator is not marked as EOF!
     */
    if (t >= mhdrvec[2].mepoch && !mit.is_eof()) {
      /* advance to next header/block */
      if (!mit.advance()) {
        if (t >= mhdrvec[2].mepoch && t < mit.header().mepoch) {
          /* cool, got it! */
          swap(0,1);
          swap(1,2);
          mhdrvec[2] = mit.header();
          int j = mit.collect(mcsvec[2]);
          mit.advance();
          return (j>0)?(-1):(1);
        } else {
          /* this block doesn't fit either! keep on reading */
          int error = 0;
          while (t <= mit.header().mepoch && (!error)) {
            /* store coeffs at 0 index */
            error = mit.collect(mcsvec[0]);
            mhdrvec[0] = mit.header();
            error = mit.advance();
          }
          if (error) {
            fprintf(stderr,
                    "[ERROR] Failed collecting coeffs for MJD %ld+%.6f from "
                    "AOD1B file %s (traceback: %s)\n",
                    t.imjd(), t.fractional_days(), mit.aod1b().fn().c_str(),
                    __func__);
            return -2;
          }
          /* read next block so that we have the bounding interval */
          error = mit.collect(mcsvec[1]);
          mhdrvec[1] = mit.header();
          /* since we are here, read also the next block */
          int j = mit.advance();
          if (j < 0) {
            mit.set_eof();
            mhdrvec[2] = mit.header();
          } else if (j > 0) {
            fprintf(stderr,
                    "[ERROR] Failed collecting coeffs for MJD %ld+%.6f from "
                    "AOD1B file %s (traceback: %s)\n",
                    t.imjd(), t.fractional_days(), mit.aod1b().fn().c_str(),
                    __func__);
            return -3;
          } else {
            error = mit.collect(mcsvec[2]);
            mhdrvec[2] = mit.header();
            mit.advance();
          }
          /* return the index */
          return 0;
        }
      } else {
        fprintf(stderr,
                "[ERROR] Failed initializing Aod1bDataStream (traceback: %s)\n",
                __func__);
        return -1;
      }
    } else if (t < mhdrvec[0].mepoch) {
      /* before anything else, check if the given t is after the first dataset
       * in the file
       */
      if (t < mit.aod().first_epoch()) {
        fprintf(stderr, "[ERROR] Failed collecting coeffs for MJD %ld+%.6f from AOD1B file %s (traceback: %s)\n",
        t.imjd(), t.fractional_days(), mit.aod1b().fn().c_str(),__func__);
        return -1;
      }
      /* shit, we need to go back in the file ! */
      if (initialize()) {
        fprintf(stderr,
                "[ERROR] Failed initializing Aod1bDataStream (traceback: %s)\n",
                __func__);
        return -1;
      }
      /* hunt interval from the top of the file */
      return hunt_range(t);
    } else if (t >= mhdrvec[2].mepoch && mit.is_eof()) {
      /* reached end-of-file, cannot collect any more data sets! */
        fprintf(stderr, "[ERROR] Failed collecting coeffs for MJD %ld+%.6f from AOD1B file %s (traceback: %s)\n",
        t.imjd(), t.fractional_days(), mit.aod1b().fn().c_str(),__func__);
        return -1;
    } else {
      /* don't know why i am here! */
      assert(1==2);
    }
  }

public:

  int initialize() noexcept {
    mit.set_begin();
    if (mit.collect(mcsvec[0])) {
      fprintf(stderr, "[ERROR] Failed initializing Aod1bDataStream (traceback: %s)\n", __func__);
      return 1;
    }
    mhdrvec[0] = mit.header();

    int j=0;
    for (int i=1; i<AR_SIZE && (!j); i++) {
      j = mit.advance();
      j+= mit.collect(mcsvec[i]);
      mhdrvec[i] = mit.header();
    }

    if (j) {
      fprintf(stderr, "[ERROR] Failed initializing Aod1bDataStream (traceback: %s)\n", __func__);
      return 1;
    }
    return 0;
  }
  
  int coefficients_at(const Datetime<nanoseconds> &t, StokesCoeffs &cs) noexcept {
    /* find suitable range (i.e. [0,1] of [1,2] */
    int index;
    if (t >= mhdrvec[0].mepoch && t < mhdrvec[1].mepoch) {
      index = 0;
    } else {
      if ((index = hunt_range(t)) < 0) {
        fprintf(stderr, "[ERROR] Failed getting Stokes coefficients from AOD1B file (traceback: %s)\n", __func__);
        return 1;
      }
    }
    /* linear interpolation */
    const auto t1 = mhdrvec[index].mepoch;
    const auto t2 = mhdrvec[index+1].mepoch;
    const double t2ti = t2.diff<DateTimeDifferenceType::FractionalDays>(t);
    const double t2t1 = t2.diff<DateTimeDifferenceType::FractionalDays>(t1);
    const double tit1 = t.diff<DateTimeDifferenceType::FractionalDays>(t1);
    const double fac1 = t2ti/t2t1;
    const double fac2 = tit1/t2t1;
    for (int m=0; m<=cs.max_order(); m++) {
      for (int l=m; l<=cs.max_degree(); l++) {
        cs.C(l,m) = fac1 * mcsvec[index].C(l,m) + fac2 * mcsvec[index+1].C(l,m);
        cs.S(l,m) = fac1 * mcsvec[index].S(l,m) + fac2 * mcsvec[index+1].S(l,m);
      }
    }
    return 0;
  }

  /** Constructor given the (file)name of an AOD1B file */
  Aod1bDataStream(const Aod1bIn &aod1b)
      : mit(aod1b), mcsvec(static_cast<StokesCoeffs *>(operator new[](
                        AR_SIZE * sizeof(StokesCoeffs)))) {
    for (int i = 0; i < AR_SIZE; i++)
      new (&mcsvec[i])
          StokesCoeffs(mit.aod1b().max_degree(), mit.aod1b().max_degree(),
                       mit.aod1b().GM(), mit.aod1b().Re());
  }
  ~Aod1bDataStream() noexcept {
    for (int i = 0; i < AR_SIZE; i++) mcsvec[i].~StokesCoeffs();
    operator delete[](mcsvec);
  }
  Aod1bDataStream(const Aod1bDataStream &other) noexcept = delete;
  Aod1bDataStream &operator=(const Aod1bDataStream &other) noexcept = delete;
  Aod1bDataStream(Aod1bDataStream &&other) noexcept
      : mit(std::move(other.mit)), mhdrvec(std::move(other.mhdrvec)),
        mcsvec(other.mcsvec) {
    other.mcsvec = nullptr;
  }
  Aod1bDataStream& operator=(Aod1bDataStream &&other) noexcept {
    if (this != &other) {
      mit = std::move(other.mit);
      mhdrvec = std::move(other.mhdrvec);
      mcsvec = other.mcsvec;
      other.mcsvec = nullptr;
    }
    return *this;
  }
}; /* Aod1bDataStream */

} /* namespace dso */

#endif
