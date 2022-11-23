#include "tropo.hpp"
#include <algorithm>
#include <cstring>
#if defined(__GNUC__) && (__GNUC__ >= 11)
#  include <charconv>
#else
#  include <cstdlib>
#  include <cerrno>
#endif 
#include <stdexcept>

dso::vmf3_details::SiteVMF3GRRecord::SiteVMF3GRRecord(
    const char *site_name, int site_length) noexcept {
  if (site_name) {
    if (site_length < 0) {
      std::strcpy(site, site_name);
    } else {
      std::strncpy(site, site_name, site_length);
      site[site_length] = '\0';
    }
  }
}

int dso::vmf3_details::parse_v3gr_line(
    const char *line, dso::vmf3_details::SiteVMF3GRRecord &rec) noexcept {
  // copy site name
#if defined(__GNUC__) && (__GNUC__ >= 8)
#  pragma GCC diagnostic push
#  pragma GCC diagnostic ignored "-Wstringop-truncation"
  std::strncpy(rec.site, line, 10);
#  pragma GCC diagnostic pop
#else
  std::strncpy(rec.site, line, 10);
#endif
  {
    // probably there are some extra trailing whitespace characters ...
    int i = 9;
    while (i > 0 && rec.site[i] == ' ')
      rec.site[i--] = '\0';
  }

#if defined(__GNUC__) && (__GNUC__ < 11)
  // for this version we are using errno to signal errors in string to double
  // conversions. Clear it
  errno = 0;
#endif

  double nums[12];
#if defined(__GNUC__) && (__GNUC__ >= 11)
  const auto sz = std::strlen(line);
  const char *end = line + sz;
#else
  char *end;
#endif
  const char *start = line + 10;
  int error = 0;

  // resolve numeric values
  for (int i = 0; i < 12; i++) {
    // note that <charconv> is only available after gcc 7.7
#if defined(__GNUC__) && (__GNUC__ >= 11)
    const auto res = std::from_chars(start, end, nums[i]);
    error += (end == start) + (res.ec != std::errc{});
#else
    nums[i] = std::strtod(start, &end);
    error += (end == start) && (errno == 0);
#endif

#ifdef DEBUG
    if (error) {
      fprintf(stderr,
              "[ERROR] Failed parsing numeric value nr %d from line \"%s\" "
              "(traceback: %s)\n",
              i + 1, line, __func__);
#if defined(__GNUC__) && (__GNUC__ < 11)
      // for this version we are using errno to signal errors in string to
      // double conversions. Clear it
      errno = 0;
#endif
    }
#endif

#if defined(__GNUC__) && (__GNUC__ >= 11)
    start = res.ptr;
#else
    start = end;
#endif
    while (*start && *start == ' ')
      ++start;
  }

  if (!error) {
    rec.mrec.mjd = nums[0];
    rec.mrec.ah = nums[1];   ///< hydrostatic "a" coefficient
    rec.mrec.aw = nums[2];   ///< wet "a" coefficient
    rec.mrec.zhd = nums[3];  ///< zenith hydrostatic delay [m]
    rec.mrec.zwd = nums[4];  ///< zenith wet delay [m]
    rec.mrec.pres = nums[5]; ///< pressure at the site [hPa]
    rec.mrec.temp = nums[6]; ///< temperature at the site [°C]
    rec.mrec.wvp = nums[7];  ///< water vapor pressure at the site [hPa]
    rec.mrec.hng = nums[8];  ///< hydrostatic north gradient Gn_h [mm]
    rec.mrec.heg = nums[9];  ///< hydrostatic east gradient Ge_h [mm]
    rec.mrec.wng = nums[10]; ///< wet north gradient Gn_w [mm]
    rec.mrec.weg = nums[11]; ///< wet east gradient Ge_w [mm]
  }

  return error;
}

int dso::vmf3_details::SiteVMF3Records::interpolate(
    const dso::datetime<dso::nanoseconds> &t,
    dso::vmf3_details::SiteVMF3GRMeteoRecord &irecs) noexcept {
  const double cmjd = t.as_mjd();

  // it might happen that we have not collecetd this site, i.e. site is
  // missing
  // assert(cmjd >= meteo_arr[0].mjd && cmjd <= meteo_arr[1].mjd);
  if (!(cmjd >= meteo_arr[0].mjd && cmjd <= meteo_arr[1].mjd)) {
    fprintf(stderr,
            "[ERROR] Failed to interpolate VMF3 for MJD=%.9f (traceback: %s)\n",
            cmjd, __func__);
    return 1;
  }

  const double x1x0 = meteo_arr[1].mjd - meteo_arr[0].mjd;
  const double xx0 = -meteo_arr[0].mjd + cmjd;
  const double x1x = meteo_arr[1].mjd - cmjd;

  irecs.mjd = cmjd;
  irecs.ah = (x1x * meteo_arr[0].ah + xx0 * meteo_arr[1].ah) / x1x0;
  irecs.aw = (x1x * meteo_arr[0].aw + xx0 * meteo_arr[1].aw) / x1x0;
  irecs.zhd = (x1x * meteo_arr[0].zhd + xx0 * meteo_arr[1].zhd) / x1x0;
  irecs.zwd = (x1x * meteo_arr[0].zwd + xx0 * meteo_arr[1].zwd) / x1x0;
  irecs.pres = (x1x * meteo_arr[0].pres + xx0 * meteo_arr[1].pres) / x1x0;
  irecs.temp = (x1x * meteo_arr[0].temp + xx0 * meteo_arr[1].temp) / x1x0;
  irecs.wvp = (x1x * meteo_arr[0].wvp + xx0 * meteo_arr[1].wvp) / x1x0;
  irecs.hng = (x1x * meteo_arr[0].hng + xx0 * meteo_arr[1].hng) / x1x0;
  irecs.heg = (x1x * meteo_arr[0].heg + xx0 * meteo_arr[1].heg) / x1x0;
  irecs.wng = (x1x * meteo_arr[0].wng + xx0 * meteo_arr[1].wng) / x1x0;
  irecs.weg = (x1x * meteo_arr[0].weg + xx0 * meteo_arr[1].weg) / x1x0;

  return 0;
}

int dso::SiteVMF3Feed::find_site(const char *site) const noexcept {
  auto it = std::find_if(recs.begin(), recs.end(),
                         [&](const vmf3_details::SiteVMF3Records &r) {
                           return !std::strncmp(site, r.site, 4);
                         });
  return (it == recs.end()) ? -1 : std::distance(recs.begin(), it);
}

int dso::SiteVMF3Feed::interpolate(
    const char *site, const dso::datetime<dso::nanoseconds> &t,
    dso::vmf3_details::SiteVMF3GRMeteoRecord &imrec) noexcept {

  // already in interval
  if (t.as_mjd() >= ct1 && t.as_mjd() <= ct2) {
    int idx = find_site(site);
    assert(idx >= 0 && idx < (int)recs.size());
    if (recs[idx].interpolate(t, imrec)) {
      fprintf(stderr,
              "[ERROR] Failed to interpolate for site: %s (traceback: %s)\n",
              site, __func__);
      return 2;
    }
    return 0;
  }

  // search for interval and then interpolate
  if (this->feed(t) > 0) {
    fprintf(stderr,
            "[ERROR] Failed to find suitable interval for interpolation! "
            "mjd=%.9f (traceback: %s)\n",
            t.as_mjd(), __func__);
    return 1;
  }

  int idx = find_site(site);
  assert(idx >= 0 && idx < (int)recs.size());
  if (recs[idx].interpolate(t, imrec)) {
    fprintf(stderr,
            "[ERROR] Failed to interpolate for site: %s (traceback: %s)\n",
            site, __func__);
    return 2;
  }
  return 0;
}

dso::SiteVMF3Feed::SiteVMF3Feed(const char *fn,
                                std::vector<const char *> &sites)
    : fin(fn) {
  cline[0] = '\0';
  std::strcpy(filename, fn);
  if (!fin.is_open()) {
    fprintf(stderr,
            "[ERROR] Failed locating VMF3 grid file %s (traceback: %s)\n", fn,
            __func__);
    throw std::runtime_error(__func__);
  }

  int sz = sites.size();
  recs.reserve(sz);

  for (int i = 0; i < sz; i++)
    recs.emplace_back(sites[i]);
}

int dso::SiteVMF3Feed::get_sites_for_current_epoch(int store_at_index,
                                                   int &sites_collected,
                                                   double &next_mjd) noexcept {
  const int idx = store_at_index;
  int error = 0;
  dso::vmf3_details::SiteVMF3GRRecord rec;
  sites_collected = 0;

  // we need current epoch; either get it from the current (already read)
  // line, or read next line from stream
  if (!cline[0])
    while (fin.getline(cline, 256) && cline[0] == '#')
      ;
  if (parse_v3gr_line(cline, rec)) {
    fprintf(stderr,
            "[ERROR] Failed to parse input line \"%s\" (traceback: %s)\n",
            cline, __func__);
    return 1;
  }
  const double current_epoch = rec.mrec.mjd;

  // keep on reading untill we come across a line with a different epoch
  while (!error) {
    // resolve line
    if (parse_v3gr_line(cline, rec)) {
      ++error;
      fprintf(stderr,
              "[ERROR] Failed to parse input line \"%s\" (traceback: %s)\n",
              cline, __func__);
    }
    // are we ok with the epoch?
    if (rec.mrec.mjd > current_epoch)
      break;
    // search for site and add
    if (int j = find_site(rec.site); j >= 0) {
      recs[j].meteo_arr[idx] = rec.mrec;
      ++sites_collected;
    }
    // read next line
    if (!fin.getline(cline, 256))
      ++error;
  }

  // next mjd
  next_mjd = rec.mrec.mjd;

  return error;
}

/// Mjd intervals in the file are every .25 of day, that is every 6 hours
int dso::SiteVMF3Feed::interval(
    const dso::datetime<dso::nanoseconds> &t) noexcept {
  const double fhours = t.sec().to_fractional_hours();
  return std::floor(fhours / 4e0);
}

int dso::SiteVMF3Feed::feed(const dso::datetime<dso::nanoseconds> &t) noexcept {
  // already in interval
  if (t.as_mjd() >= ct1 && t.as_mjd() <= ct2)
    return 0;

  // do we need to go back in the file ?
  // seconds check is to avoid infinit looping on first feed
  if (t.as_mjd() < ct1 && ct1 != std::numeric_limits<double>::max()) {
    cline[0] = '\0';
    fin.seekg(0);
    return feed(t);
  }

  bool sites_missing = false;

  // we need to go forward ...
  int num_sites = recs.size();
  int error = 0;
  while (!error) {
    // swap contents of the meteo arrays, so that the most recent meteo
    // records are the initial entries in the arrays
    for (auto it = recs.begin(); it != recs.end(); ++it) {
      it->meteo_arr[0] = it->meteo_arr[1];
    }
    ct1 = recs.begin()->meteo_arr[0].mjd;
    // read in current epoch from file, store them at index 1 of the
    // respective arrays
    double next_mjd;
    int sites_collected;
    error = get_sites_for_current_epoch(1, sites_collected, next_mjd);
    sites_missing += (sites_collected == num_sites) ? 0 : 100;
    if (sites_missing) {
      fprintf(stderr,
              "[ERROR] Failed to collect VMF3 info for all sites (%d/%d) "
              "(traceback: %s)\n",
              sites_collected, num_sites, __func__);
    }
    ct2 = recs.begin()->meteo_arr[1].mjd;
    // are we ok now ?
    if (t.as_mjd() >= ct1 && t.as_mjd() <= ct2)
      break;
  }

  if (error) {
    fprintf(stderr,
            "[ERROR] An error occured while reading VMF3 grid file \"%s\" "
            "(traceback: %s)\n",
            filename, __func__);
    return error;
  }

  return sites_missing ? -1 : 0;
}
