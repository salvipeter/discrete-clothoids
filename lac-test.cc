#include "discrete-lac.hh"

#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <iterator>

using namespace Geometry;

int main(int argc, char **argv) {
  Point2DVector points = { {250.0, 150.0}, {265.0, 277.0}, {178.0, 350.0},
                           {249.0, 391.0}, {353.0, 314.0}, {458.0, 354.0},
                           {457.0, 269.0}, {340.0, 221.0}, {321.0, 109.0} };
  DiscreteLAC fitter;
  fitter.setClosed(true);
  fitter.setAlpha(-1.0);
  fitter.setSamplingRate(20);
  fitter.setTolerance(1.0e-5);
  fitter.setPoints(points);
  fitter.fit();

  const auto &polyline = fitter.polyline();
  std::ofstream f("polyline.txt");
  std::copy(polyline.begin(), polyline.end(), std::ostream_iterator<Point2D>(f, "\n"));
  f << polyline.front() << std::endl; // closed polyline
  f.close();

  double scaling = 20;
  auto tangents = fitter.tangents();
  f.open("tangents.txt");
  for (size_t i = 0; i < points.size(); ++i) {
    const auto &p = points[i];
    const auto &t = tangents[i];
    f << p << std::endl;
    f << p + t * scaling << std::endl;
    f << p - t * scaling << std::endl;
    f << p << std::endl;
  }
  f.close();

  // Display the result with GNUPlot
  std::system("gnuplot lac-test.p");
}
