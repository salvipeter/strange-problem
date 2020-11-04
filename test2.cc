#include <fstream>

#include <Eigen/LU>

#include <geometry.hh>

using namespace Geometry;

BSSurface fitSimilarBSpline(BSSurface surface, const Point2DVector &params,
                            const PointVector &points) {
  // Retains the first V control row
  size_t n = points.size();
  auto [mu, mv] = surface.numControlPoints();
  const auto &base_u = surface.basisU();
  const auto &base_v = surface.basisV();

  Eigen::MatrixXd A = Eigen::MatrixXd::Zero(n, (mu - 1) * mv);
  Eigen::MatrixXd b(n, 3);
  using VecMap = Eigen::Map<const Eigen::Vector3d>;

  for (size_t i = 0; i < n; ++i) {
    DoubleVector coeff_u, coeff_v;
    size_t span_u = base_u.findSpan(params[i][0]);
    size_t span_v = base_v.findSpan(params[i][1]);
    base_u.basisFunctions(span_u, params[i][0], coeff_u);
    base_v.basisFunctions(span_v, params[i][1], coeff_v);
    size_t offset_u = span_u - base_u.degree();
    size_t offset_v = span_v - base_v.degree();
    for (size_t j = offset_u == 0 ? 1 : 0; j <= base_u.degree(); ++j)
      for (size_t k = 0; k <= base_v.degree(); ++k)
        A(i, (offset_u + j - 1) * mv + offset_v + k) = coeff_u[j] * coeff_v[k];

    b.block<1,3>(i, 0) = VecMap(points[i].data());
    if (offset_u == 0)
      for (size_t k = 0; k <= base_v.degree(); ++k) {
        auto p = surface.controlPoint(0, offset_v + k);
        b.block<1,3>(i, 0) -= VecMap((p * coeff_u[0] * coeff_v[k]).data());
      }
  }

  Eigen::MatrixXd x = A.fullPivLu().solve(b);
  for (size_t j = 1; j < mu; ++j)
    for (size_t k = 0; k < mv; ++k) {
      size_t i = (j - 1) * mv + k;
      surface.controlPoint(j, k) = { x(i, 0), x(i, 1), x(i, 2) };
    }

  return surface;
}

BSSurface readBSS(std::string filename) {
  std::ifstream f(filename);
  f.exceptions(std::ios::failbit | std::ios::badbit);

  size_t nu, nv, du, dv;
  double u, v, x, y, z;
  DoubleVector knots_u, knots_v;
  PointVector cpts;

  f >> du >> dv >> nu >> nv;

  size_t n_knots_u = nu + du + 1, n_knots_v = nv + dv + 1;
  for (size_t i = 0; i < n_knots_u; ++i) {
    f >> u;
    knots_u.push_back(u);
  }
  for (size_t j = 0; j < n_knots_v; ++j) {
    f >> v;
    knots_v.push_back(v);
  }

  for (size_t i = 0; i < nu; ++i)
    for (size_t j = 0; j < nv; ++j) {
      f >> x >> y >> z;
      cpts.emplace_back(x, y, z);
    }
  
  return { du, dv, knots_u, knots_v, cpts };
}

void writeBSS(const BSSurface &s, std::string filename) {
  std::ofstream f(filename);
  f.exceptions(std::ios::failbit | std::ios::badbit);

  f << s.basisU().degree() << ' ' << s.basisV().degree() << std::endl;
  auto num_cp = s.numControlPoints();
  f << num_cp[0] << ' ' << num_cp[1] << std::endl;

  for (double u : s.basisU().knots())
    f << u << ' ';
  f << std::endl;
  for (double v : s.basisV().knots())
    f << v << ' ';
  f << std::endl;

  for (const auto &p : s.controlPoints())
    f << p << std::endl;
}

void writeControlNet(const std::vector<BSSurface> &surfaces, std::string filename) {
  std::ofstream f(filename);
  f.exceptions(std::ios::failbit | std::ios::badbit);
  std::vector<size_t> start_indices;
  size_t index = 1;
  for (const auto &s : surfaces) {
    start_indices.push_back(index);
    for (const auto &p : s.controlPoints())
      f << "v " << p << std::endl;
    index += s.controlPoints().size();
  }
  for (size_t i = 0; i < surfaces.size(); ++i) {
    const auto &s = surfaces[i];
    size_t base = start_indices[i];
    auto [nu, nv] = s.numControlPoints();
    for (size_t j = 0; j < nu; ++j)
      for (size_t k = 1; k < nv; ++k)
        f << "l " << base + j * nv + k - 1 << ' ' << base + j * nv + k << std::endl;
    for (size_t j = 1; j < nu; ++j)
      for (size_t k = 0; k < nv; ++k)
        f << "l " << base + (j - 1) * nv + k << ' ' << base + j * nv + k << std::endl;
  }
}

void writePoints(const PointVector &points, std::string filename) {
  std::ofstream f(filename);
  f.exceptions(std::ios::failbit | std::ios::badbit);
  for (const auto &p : points)
    f << "v " << p[0] << ' ' << p[1] << ' ' << p[2] << std::endl;
  for (size_t i = 1; i <= points.size(); ++i)
    f << "p " << i << std::endl;
}

auto readParamPoints(std::string filename, int type) {
  Point2DVector params;
  PointVector points;

  std::ifstream f(filename);
  f.exceptions(std::ios::failbit | std::ios::badbit);

  size_t n;
  f >> n;
  for (size_t i = 0; i < n; ++i) {
    double u, v, x, y, z;
    f >> u >> v >> x >> y >> z;
    switch (type) {
    case 0: // default
      params.emplace_back(u, v);
      points.emplace_back(x, y, z);
      break;
    case 1: // symmetric
      if (y < 0) {
        params.emplace_back(u, v);
        points.emplace_back(x, y, z);
        params.emplace_back(u, 1 - v);
        points.emplace_back(x, -y, z);
      }
      break;
    case 2: // parabolic
      params.emplace_back(u, v);
      double u1 = (2 * v * v - 2 * v + 1) * u;
      points.emplace_back(-100 + u1 * 20, -100 + v * 200, u1 * 40);
      break;
    }
  }

  return std::make_pair(params, points);
}

int main(int argc, char **argv) {
  auto surface = readBSS("surface.bss");
  auto [params, points] = readParamPoints("parampoints.txt", 0);
  auto result = fitSimilarBSpline(surface, params, points);
  writePoints(points, "points.obj");
  writeControlNet({ result }, "result.obj");
}
