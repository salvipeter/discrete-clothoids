#include "discrete-lac.hh"

#include <algorithm>
#include <cmath>
#include <exception>

using namespace Geometry;

DiscreteLAC::DiscreteLAC()
  : closed_(false), alpha_(-1), sampling_rate_(20), tolerance_(epsilon)
{
}

bool DiscreteLAC::closed() const {
  return closed_;
}

void DiscreteLAC::setClosed(bool b) {
  closed_ = b;
}

double DiscreteLAC::alpha() const {
  return alpha_;
}

void DiscreteLAC::setAlpha(double x) {
  alpha_ = x;
}

size_t DiscreteLAC::samplingRate() const {
  return sampling_rate_;
}

void DiscreteLAC::setSamplingRate(size_t n) {
  sampling_rate_ = n;
}

double DiscreteLAC::tolerance() const {
  return tolerance_;
}

void DiscreteLAC::setTolerance(double x) {
  tolerance_ = x;
}

const Point2DVector &DiscreteLAC::points() const {
  return input_;
}

void DiscreteLAC::setPoints(const Point2DVector &points) {
  input_ = points;
}

const Point2DVector &DiscreteLAC::polyline() const {
  return polyline_;
}


// Interesting functions

static double approximateCurvature(const Point2D &prev, const Point2D &p, const Point2D &next) {
  auto denom = (p - prev).norm() * (next - p).norm() * (next - prev).norm();
  auto a = p - prev, b = next - p;
  return 2 * (a[0] * b[1] - a[1] * b[0]) / denom;
}

static Point2D updatePoint(const Point2D &prev, const Point2D &p, const Point2D &next,
                           double target) {
  target *= (p - prev).norm() * (next - p).norm() * (next - prev).norm() / 2;
  target += prev[0] * next[1] - prev[1] * next[0];
  double beta_x = next[1] - prev[1], beta_y = prev[0] - next[0];
  auto start = (next + prev) / 2;
  auto dir = Vector2D(prev[1] - next[1], next[0] - prev[0]).normalize();
  target -= start[0] * beta_x + start[1] * beta_y;
  return start + dir * target / (dir[0] * beta_x + dir[1] * beta_y);
}

void DiscreteLAC::fit() {
  if (input_.size() < 3)
    throw std::runtime_error("Not enough input points (needs at least 3)");

  // Subsampling
  // -----------

  // Compute the average distance & sampling density
  size_t n = input_.size();
  avg_distance_ = 0;
  DoubleVector distances;
  for (size_t i = 1; i < n; ++i) {
    distances.push_back((input_[i] - input_[i-1]).norm());
    avg_distance_ += distances.back();
  }
  if (closed_) {
    distances.push_back((input_[0] - input_.back()).norm());
    avg_distance_ = (avg_distance_ + distances.back()) / n;
  } else
    avg_distance_ /= n - 1;
  double density = sampling_rate_ / avg_distance_;
  avg_distance_ /= sampling_rate_;

  // Fill the polyline with sampled points
  polyline_.clear();
  polyline_.push_back(input_[0]);
  input_indices_.push_back(0);
  for (size_t i = 1; i < n; ++i) {
    size_t resolution = std::round(distances[i-1] * density);
    for (size_t j = 1; j < resolution; ++j) {
      double u = (double)j / resolution;
      polyline_.push_back(input_[i-1] * (1 - u) + input_[i] * u);
    }
    input_indices_.push_back(polyline_.size());
    polyline_.push_back(input_[i]);
  }
  if (closed_) {
    size_t resolution = std::round(distances.back() * density);
    for (size_t j = 1; j < resolution; ++j) {
      double u = (double)j / resolution;
      polyline_.push_back(input_.back() * (1 - u) + input_.front() * u);
    }
  }

  // Fitting iteration
  // -----------------

  DoubleVector curvature(polyline_.size());
  double max_change;
  do {
    // Curvature values at seed points
    if (closed_)
      curvature[0] = approximateCurvature(polyline_.back(), polyline_[0], polyline_[1]);
    else
      curvature[0] = 0;
    for (size_t i = 1; i < n - 1; ++i) {
      size_t j = input_indices_[i];
      curvature[j] = approximateCurvature(polyline_[j-1], polyline_[j], polyline_[j+1]);
    }
    if (closed_) {
      size_t j = input_indices_[n-1];
      curvature[j] = approximateCurvature(polyline_[j-1], polyline_[j], polyline_[j+1]);
    } else
      curvature.back() = 0;

    // Generate target curvature
    auto linearize = [&](double x) {
                       return std::pow(std::abs(x), -alpha_) * (x < 0 ? -1 : 1);
                     };
    auto delinearize = [&](double x) {
                         return std::pow(std::abs(x), -1 / alpha_) * (x < 0 ? -1 : 1);
                       };
    for (size_t i = 0; i < n; ++i)
      curvature[input_indices_[i]] = linearize(curvature[input_indices_[i]]);
    for (size_t i = 1; i < n; ++i) {
      size_t from = input_indices_[i-1], to = input_indices_[i];
      size_t resolution = to - from;
      for (size_t j = 1; j < resolution; ++j) {
        double u = (double)j / resolution;
        curvature[from+j] = curvature[from] * (1 - u) + curvature[to] * u;
      }
    }
    if (closed_) {
      size_t from = input_indices_[n-1], to = input_indices_[0];
      size_t resolution = curvature.size() - from;
      for (size_t j = 1; j < resolution; ++j) {
        double u = (double)j / resolution;
        curvature[from+j] = curvature[from] * (1 - u) + curvature[to] * u;
      }
    }
    std::transform(curvature.begin(), curvature.end(), curvature.begin(), delinearize);

    // Update points
    auto tmp = polyline_;
    for (size_t i = 1; i < n; ++i)
      for (size_t j = input_indices_[i-1] + 1; j < input_indices_[i]; ++j)
        polyline_[j] = updatePoint(tmp[j-1], tmp[j], tmp[j+1], curvature[j]);
    if (closed_) {
      for (size_t j = input_indices_[n-1] + 1; j < polyline_.size() - 1; ++j)
        polyline_[j] = updatePoint(tmp[j-1], tmp[j], tmp[j+1], curvature[j]);
      polyline_.back() = updatePoint(tmp.rbegin()[1], tmp.back(), tmp[0], curvature.back());
    }

    // Compute maximal change
    max_change = 0;
    for (size_t i = 0; i < tmp.size(); ++i) {
      double change = (tmp[i] - polyline_[i]).norm();
      if (change > max_change)
        max_change = change;
    }
  } while (max_change > avg_distance_ * tolerance_);
}

Vector2DVector DiscreteLAC::tangents() const {
  Vector2DVector result;
  Vector2D n;
  for (auto i : input_indices_) {
    if (i == 0) {
      if (closed_)
        n = polyline_[1] - polyline_.back();
      else
        n = polyline_[1] - polyline_[0];
    } else if (i == polyline_.size() - 1)
      n = polyline_[i] - polyline_[i-1];
    else
      n = polyline_[i+1] - polyline_[i-1];
    result.push_back(n.normalize());
  }
  return result;
}
