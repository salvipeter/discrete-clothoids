#include <geometry.hh>                             // https://github.com/salvipeter/libgeom

class DiscreteLAC {
public:
  DiscreteLAC();

  // Properties
  bool closed() const;                             // generate a closed curve?
  void setClosed(bool b);
  double alpha() const;                            // LAC alpha (-1 => Euler spiral etc.)
  void setAlpha(double x);
  size_t samplingRate() const;                     // average # of samples between input points
  void setSamplingRate(size_t n);
  double tolerance() const;                        // relative tol. for ending the iteration
  void setTolerance(double x);
  const Geometry::Point2DVector &points() const;   // input points
  void setPoints(const Geometry::Point2DVector &points);

  // Fitting
  void fit();

  // Results
  const Geometry::Point2DVector &polyline() const;
  Geometry::Vector2DVector tangents() const;

private:
  bool closed_;
  double alpha_;
  size_t sampling_rate_;
  double tolerance_;
  Geometry::Point2DVector input_, polyline_;
  std::vector<size_t> input_indices_;
  double avg_distance_;
};
