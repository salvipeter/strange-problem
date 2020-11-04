#include <fstream>

#include <Eigen/LU>

#include <Point.hh>
#include <Vector.hh>
#include <BSplineSurface.hh>

namespace S {
  using Point2D = Point<2, double>;
  using PointType = Point<3, double>;
  using VectorType = Vector<3, double>;
  using BSplineSurfaceType = BSplineSurface<PointType>;
}

S::BSplineSurfaceType fitSimilarBSpline(S::BSplineSurfaceType surface,
                                        const std::vector<S::Point2D> &params,
                                        const std::vector<S::PointType> &points) {
  // Retains the first V control row
  size_t n = points.size();
  size_t mu = surface.theControlNet().sizes()[0];
  size_t mv = surface.theControlNet().sizes()[1];
  const auto &base_u = surface.BaseFunctionsU();
  const auto &base_v = surface.BaseFunctionsV();

  Eigen::MatrixXd A = Eigen::MatrixXd::Zero(n, (mu - 1) * mv);
  Eigen::MatrixXd b(n, 3);
  using VecMap = Eigen::Map<const Eigen::Vector3d>;

  for (size_t i = 0; i < n; ++i) {
    const double *coeff_u = base_u.Eval(params[i][0]);
    const double *coeff_v = base_v.Eval(params[i][1]);
    size_t offset_u = base_u.Offset();
    size_t offset_v = base_v.Offset();
    for (size_t j = offset_u == 0 ? 1 : 0; j < base_u.NrValues(); ++j)
      for (size_t k = 0; k < base_v.NrValues(); ++k)
        A(i, (offset_u + j - 1) * mv + offset_v + k) = coeff_u[j] * coeff_v[k];

    b.block<1,3>(i, 0) = VecMap(points[i].Base());
    if (offset_u == 0)
      for (size_t k = 0; k < base_v.NrValues(); ++k) {
        auto p = S::VectorType(surface.theControlNet()[0][offset_v+k]);
        b.block<1,3>(i, 0) -= VecMap((p * coeff_u[0] * coeff_v[k]).Base());
      }
  }

  Eigen::MatrixXd x = A.fullPivLu().solve(b);
  for (size_t j = 1; j < mu; ++j)
    for (size_t k = 0; k < mv; ++k) {
      size_t i = (j - 1) * mv + k;
      surface.theControlNet()[j][k] = S::PointType(x(i, 0), x(i, 1), x(i, 2));
    }

  return surface;
}

S::BSplineSurfaceType readBSS(std::string filename) {
  std::ifstream f(filename);
  f.exceptions(std::ios::failbit | std::ios::badbit);

  UInt nu, nv, du, dv;
  double u, v, x, y, z;
  std::vector<double> knots_u, knots_v;

  f >> du >> dv >> nu >> nv;
  size_t size[2] = { nu, nv };
  S::BSplineSurfaceType::ControlNetType cpts(size);

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
      cpts[i][j] = { x, y, z };
    }
  
  return {
    du, dv,
    S::BSplineSurfaceType::KnotVectorType(knots_u.begin(), knots_u.end()),
    S::BSplineSurfaceType::KnotVectorType(knots_v.begin(), knots_v.end()),
    cpts
  };
}

void writeBSS(const S::BSplineSurfaceType &s, std::string filename) {
  std::ofstream f(filename);
  f.exceptions(std::ios::failbit | std::ios::badbit);

  f << s.DegreeU() << ' ' << s.DegreeV() << std::endl;
  auto num_cp = s.theControlNet().sizes();
  f << num_cp[0] << ' ' << num_cp[1] << std::endl;

  for (double u : s.KnotVectorU())
    f << u << ' ';
  f << std::endl;
  for (double v : s.KnotVectorV())
    f << v << ' ';
  f << std::endl;

  for (size_t i = 0; i < num_cp[0]; ++i)
    for (size_t j = 0; j < num_cp[1]; ++j) {
      const auto &p = s.theControlNet()[i][j];
      f << p[0] << ' ' << p[1] << ' ' << p[2] << std::endl;
    }
}

void writeControlNet(const std::vector<S::BSplineSurfaceType> &surfaces, std::string filename) {
  std::ofstream f(filename);
  f.exceptions(std::ios::failbit | std::ios::badbit);
  std::vector<size_t> start_indices;
  size_t index = 1;
  for (const auto &s : surfaces) {
    start_indices.push_back(index);
    auto sizes = s.theControlNet().sizes();
    for (size_t i = 0; i < sizes[0]; ++i)
      for (size_t j = 0; j < sizes[1]; ++j) {
        const auto &p = s.theControlNet()[i][j];
        f << "v " << p[0] << ' ' << p[1] << ' ' << p[2] << std::endl;
      }
    index += sizes[0] * sizes[1];
  }
  for (size_t i = 0; i < surfaces.size(); ++i) {
    const auto &s = surfaces[i];
    size_t base = start_indices[i];
    auto sizes = s.theControlNet().sizes();
    size_t nu = sizes[0], nv = sizes[1];
    for (size_t j = 0; j < nu; ++j)
      for (size_t k = 1; k < nv; ++k)
        f << "l " << base + j * nv + k - 1 << ' ' << base + j * nv + k << std::endl;
    for (size_t j = 1; j < nu; ++j)
      for (size_t k = 0; k < nv; ++k)
        f << "l " << base + (j - 1) * nv + k << ' ' << base + j * nv + k << std::endl;
  }
}

void writePoints(const std::vector<S::PointType> &points, std::string filename) {
  std::ofstream f(filename);
  f.exceptions(std::ios::failbit | std::ios::badbit);
  for (const auto &p : points)
    f << "v " << p[0] << ' ' << p[1] << ' ' << p[2] << std::endl;
  for (size_t i = 1; i <= points.size(); ++i)
    f << "p " << i << std::endl;
}

auto readParamPoints(std::string filename, int type) {
  std::vector<S::Point2D> params;
  std::vector<S::PointType> points;

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
