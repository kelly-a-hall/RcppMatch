#include <Rcpp.h>
#include <vector>
#include <algorithm>

using namespace Rcpp;

struct ControlUnit {
  int idx;
  double ps;
  ControlUnit(int i, double p) : idx(i), ps(p) {}
  bool operator<(const ControlUnit &other) const {
    return ps < other.ps;
  }
};

// [[Rcpp::export]]
DataFrame match_nearest(NumericVector ps,
                        IntegerVector treat,
                        IntegerVector id,
                        double caliper = -1.0) {
  const int n = ps.size();
  std::vector<int> treated_idx;
  std::vector<ControlUnit> controls;

  treated_idx.reserve(n / 2);
  controls.reserve(n / 2);

  for (int i = 0; i < n; ++i) {
    if (treat[i] == 1)
      treated_idx.push_back(i);
    else
      controls.emplace_back(i, ps[i]);
  }

  if (treated_idx.empty() || controls.empty()) {
    return DataFrame::create();  // No matches
  }

  std::sort(controls.begin(), controls.end());

  // Optimization: use bool vector instead of unordered_set
  std::vector<bool> used_controls(n, false);

  std::vector<int> result_id;
  std::vector<int> result_treat;
  std::vector<double> result_ps;
  result_id.reserve(2 * treated_idx.size());
  result_treat.reserve(2 * treated_idx.size());
  result_ps.reserve(2 * treated_idx.size());

  for (int ti : treated_idx) {
    double ps_t = ps[ti];
    ControlUnit target(-1, ps_t);
    auto it = std::lower_bound(controls.begin(), controls.end(), target);

    int best_idx = -1;
    double min_dist = std::numeric_limits<double>::max();

    auto try_match = [&](std::vector<ControlUnit>::const_iterator match_it) {
      if (match_it != controls.end()) {
        int idx = match_it->idx;
        if (!used_controls[idx]) {
          double dist = std::abs(ps_t - match_it->ps);
          if (dist < min_dist) {
            min_dist = dist;
            best_idx = idx;
          }
        }
      }
    };

    try_match(it);
    if (it != controls.begin()) try_match(std::prev(it));

    if (best_idx != -1 && (caliper < 0.0 || min_dist <= caliper)) {
      used_controls[best_idx] = true;

      result_id.push_back(id[ti]);
      result_treat.push_back(1);
      result_ps.push_back(ps[ti]);

      result_id.push_back(id[best_idx]);
      result_treat.push_back(0);
      result_ps.push_back(ps[best_idx]);
    }
  }

  return DataFrame::create(
    Named("id") = result_id,
    Named("treat") = result_treat,
    Named("ps") = result_ps
  );
}
