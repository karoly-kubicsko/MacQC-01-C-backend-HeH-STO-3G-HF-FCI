#include "macqc/internal_coordinates.h"

#include <algorithm>
#include <cmath>
#include <numeric>
#include <stdexcept>

namespace macqc {

namespace {

constexpr double kPi = 3.14159265358979323846;
constexpr double kLinearOn = 175.0 * kPi / 180.0;   // a bend becomes a linear bend above this
constexpr double kLinearOff = 170.0 * kPi / 180.0;  // and a linear bend a bend again below this

using Vec3 = std::array<double, 3>;

Vec3 pos(const Molecule& m, int a) { return m.atoms[static_cast<std::size_t>(a)].position; }
Vec3 sub(const Vec3& a, const Vec3& b) { return {a[0] - b[0], a[1] - b[1], a[2] - b[2]}; }
double dot(const Vec3& a, const Vec3& b) { return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]; }
double norm(const Vec3& a) { return std::sqrt(dot(a, a)); }
Vec3 scale(const Vec3& a, double f) { return {a[0] * f, a[1] * f, a[2] * f}; }
Vec3 cross(const Vec3& a, const Vec3& b) {
    return {a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]};
}

double bend_angle(const Molecule& m, const InternalCoordinate& c) {
    const Vec3 u = sub(pos(m, c.atoms[0]), pos(m, c.atoms[1])), v = sub(pos(m, c.atoms[2]), pos(m, c.atoms[1]));
    const double cs = std::max(-1.0, std::min(1.0, dot(u, v) / (norm(u) * norm(v))));
    return std::acos(cs);
}

double value_of(const Molecule& m, const InternalCoordinate& c) {
    switch (c.type) {
        case InternalType::Stretch: return norm(sub(pos(m, c.atoms[1]), pos(m, c.atoms[0])));
        case InternalType::Bend: return bend_angle(m, c);
        case InternalType::LinearBend: {
            Vec3 u = sub(pos(m, c.atoms[0]), pos(m, c.atoms[1])), v = sub(pos(m, c.atoms[2]), pos(m, c.atoms[1]));
            u = scale(u, 1.0 / norm(u));
            v = scale(v, 1.0 / norm(v));
            return dot(u, c.axis) + dot(v, c.axis);
        }
    }
    return 0.0;
}

std::string label(const Molecule& m, int a) {
    return m.atoms[static_cast<std::size_t>(a)].symbol + std::to_string(a + 1);
}

}  // namespace

double covalent_radius(int Z) {
    // Pyykkö & Atsumi, Chem. Eur. J. 15, 186 (2009), single-bond radii in angstrom.
    static const double r[] = {0.32, 0.46, 1.33, 1.02, 0.85, 0.75, 0.71, 0.63, 0.64, 0.67};
    if (Z < 1 || Z > 10) return 1.5 / kBohrInAngstrom;
    return r[Z - 1] / kBohrInAngstrom;
}

std::string InternalCoordinate::name() const {
    switch (type) {
        case InternalType::Stretch: return "R";
        case InternalType::Bend: return "A";
        case InternalType::LinearBend: return component == 0 ? "La" : "Lb";
    }
    return "?";
}

std::string InternalCoordinate::definition(const Molecule& m) const {
    switch (type) {
        case InternalType::Stretch: return "R(" + label(m, atoms[0]) + "-" + label(m, atoms[1]) + ")";
        case InternalType::Bend:
            return "A(" + label(m, atoms[0]) + "-" + label(m, atoms[1]) + "-" + label(m, atoms[2]) + ")";
        case InternalType::LinearBend:
            return std::string("L(") + label(m, atoms[0]) + "-" + label(m, atoms[1]) + "-" + label(m, atoms[2]) + ", " +
                   (component == 0 ? "1" : "2") + ")";
    }
    return "";
}

InternalCoordinates::InternalCoordinates(const Molecule& m) {
    const int n = static_cast<int>(m.atoms.size());
    if (n < 2) throw std::invalid_argument("internal coordinates need at least two atoms");
    // bonds by covalent radii, then connect the fragments through their closest atoms
    std::vector<std::pair<int, int>> bonds;
    for (int i = 0; i < n; ++i)
        for (int j = i + 1; j < n; ++j) {
            const double r = norm(sub(pos(m, i), pos(m, j)));
            if (r < 1.3 * (covalent_radius(m.atoms[static_cast<std::size_t>(i)].Z) + covalent_radius(m.atoms[static_cast<std::size_t>(j)].Z)))
                bonds.emplace_back(i, j);
        }
    while (true) {
        std::vector<int> comp(static_cast<std::size_t>(n));
        std::iota(comp.begin(), comp.end(), 0);
        bool changed = true;
        while (changed) {
            changed = false;
            for (const auto& b : bonds) {
                const int c = std::min(comp[static_cast<std::size_t>(b.first)], comp[static_cast<std::size_t>(b.second)]);
                for (int x : {b.first, b.second})
                    if (comp[static_cast<std::size_t>(x)] != c) {
                        comp[static_cast<std::size_t>(x)] = c;
                        changed = true;
                    }
            }
        }
        double best = 1e300;
        std::pair<int, int> link{-1, -1};
        for (int i = 0; i < n; ++i)
            for (int j = i + 1; j < n; ++j)
                if (comp[static_cast<std::size_t>(i)] != comp[static_cast<std::size_t>(j)]) {
                    const double r = norm(sub(pos(m, i), pos(m, j)));
                    if (r < best) {
                        best = r;
                        link = {i, j};
                    }
                }
        if (link.first < 0) break;
        bonds.push_back(link);
    }
    std::sort(bonds.begin(), bonds.end());
    for (const auto& b : bonds) {
        InternalCoordinate c;
        c.type = InternalType::Stretch;
        c.atoms = std::vector<int>{b.first, b.second};  // (not "= {…}": GCC 13 -Wnonnull false positive)
        coords_.push_back(c);
    }
    // bends at every atom with two bonds
    for (int B = 0; B < n; ++B) {
        std::vector<int> nb;
        for (const auto& b : bonds) {
            if (b.first == B) nb.push_back(b.second);
            if (b.second == B) nb.push_back(b.first);
        }
        std::sort(nb.begin(), nb.end());
        for (std::size_t x = 0; x < nb.size(); ++x)
            for (std::size_t y = x + 1; y < nb.size(); ++y) {
                InternalCoordinate c;
                c.type = InternalType::Bend;
                c.atoms = std::vector<int>{nb[x], B, nb[y]};
                if (bend_angle(m, c) < kLinearOn) {
                    coords_.push_back(c);
                    continue;
                }
                // linear: two components perpendicular to the A...C axis
                Vec3 e0 = sub(pos(m, nb[y]), pos(m, nb[x]));
                e0 = scale(e0, 1.0 / norm(e0));
                int k = 0;
                for (int d = 1; d < 3; ++d)
                    if (std::fabs(e0[static_cast<std::size_t>(d)]) < std::fabs(e0[static_cast<std::size_t>(k)])) k = d;
                Vec3 w{0.0, 0.0, 0.0};
                w[static_cast<std::size_t>(k)] = 1.0;
                Vec3 e1 = sub(w, scale(e0, dot(w, e0)));
                e1 = scale(e1, 1.0 / norm(e1));
                const Vec3 e2 = cross(e0, e1);
                for (int comp_k = 0; comp_k < 2; ++comp_k) {
                    InternalCoordinate l = c;
                    l.type = InternalType::LinearBend;
                    l.component = comp_k;
                    l.axis = comp_k == 0 ? e1 : e2;
                    coords_.push_back(l);
                }
            }
    }
}

bool InternalCoordinates::needs_rebuild(const Molecule& m) const {
    for (const auto& c : coords_) {
        if (c.type == InternalType::Bend && bend_angle(m, c) > kLinearOn) return true;
        if (c.type == InternalType::LinearBend && bend_angle(m, c) < kLinearOff) return true;
    }
    return false;
}

std::vector<double> InternalCoordinates::values(const Molecule& m) const {
    std::vector<double> q;
    q.reserve(coords_.size());
    for (const auto& c : coords_) q.push_back(value_of(m, c));
    return q;
}

Matrix InternalCoordinates::wilson_b(const Molecule& m) const {
    const std::size_t n3 = 3 * m.atoms.size();
    Matrix B(coords_.size(), n3);
    for (std::size_t i = 0; i < coords_.size(); ++i) {
        const InternalCoordinate& c = coords_[i];
        if (c.type == InternalType::Stretch) {
            Vec3 d = sub(pos(m, c.atoms[1]), pos(m, c.atoms[0]));
            const double r = norm(d);
            for (std::size_t k = 0; k < 3; ++k) {
                B(i, 3 * static_cast<std::size_t>(c.atoms[0]) + k) = -d[k] / r;
                B(i, 3 * static_cast<std::size_t>(c.atoms[1]) + k) = d[k] / r;
            }
        } else if (c.type == InternalType::Bend) {
            Vec3 u = sub(pos(m, c.atoms[0]), pos(m, c.atoms[1])), v = sub(pos(m, c.atoms[2]), pos(m, c.atoms[1]));
            const double ru = norm(u), rv = norm(v);
            u = scale(u, 1.0 / ru);
            v = scale(v, 1.0 / rv);
            const double cs = std::max(-1.0, std::min(1.0, dot(u, v))), sn = std::sqrt(std::max(1e-300, 1.0 - cs * cs));
            for (std::size_t k = 0; k < 3; ++k) {
                const double da = (cs * u[k] - v[k]) / (ru * sn), dc = (cs * v[k] - u[k]) / (rv * sn);
                B(i, 3 * static_cast<std::size_t>(c.atoms[0]) + k) = da;
                B(i, 3 * static_cast<std::size_t>(c.atoms[2]) + k) = dc;
                B(i, 3 * static_cast<std::size_t>(c.atoms[1]) + k) = -(da + dc);
            }
        } else {
            const double h = 1.0e-6;
            for (std::size_t x = 0; x < n3; ++x) {
                Molecule mp = m, mm = m;
                mp.atoms[x / 3].position[x % 3] += h;
                mm.atoms[x / 3].position[x % 3] -= h;
                B(i, x) = (value_of(mp, c) - value_of(mm, c)) / (2.0 * h);
            }
        }
    }
    return B;
}

Matrix InternalCoordinates::b_derivative_contracted(const Molecule& m, const std::vector<double>& c) const {
    const std::size_t n3 = 3 * m.atoms.size();
    Matrix K(n3, n3);
    const double h = 1.0e-5;
    for (std::size_t j = 0; j < n3; ++j) {
        Molecule mp = m, mm = m;
        mp.atoms[j / 3].position[j % 3] += h;
        mm.atoms[j / 3].position[j % 3] -= h;
        const Matrix Bp = wilson_b(mp), Bm = wilson_b(mm);
        for (std::size_t k = 0; k < n3; ++k) {
            double s = 0.0;
            for (std::size_t i = 0; i < coords_.size(); ++i) s += c[i] * (Bp(i, k) - Bm(i, k)) / (2.0 * h);
            K(j, k) = s;
        }
    }
    for (std::size_t j = 0; j < n3; ++j)
        for (std::size_t k = 0; k < j; ++k) K(j, k) = K(k, j) = 0.5 * (K(j, k) + K(k, j));
    return K;
}

std::vector<double> InternalCoordinates::difference(const std::vector<double>& q1, const std::vector<double>& q0) const {
    std::vector<double> d(q1.size());
    for (std::size_t i = 0; i < q1.size(); ++i) d[i] = q1[i] - q0[i];
    return d;
}

Matrix generalized_inverse(const Matrix& G, double tol, int* rank) {
    const EigenResult e = symmetric_eigen(G);
    const std::size_t n = G.rows();
    double mx = 0.0;
    for (double v : e.values) mx = std::max(mx, std::fabs(v));
    Matrix inv(n, n);
    int r = 0;
    for (std::size_t k = 0; k < n; ++k) {
        if (e.values[k] <= tol * mx) continue;
        ++r;
        for (std::size_t i = 0; i < n; ++i)
            for (std::size_t j = 0; j < n; ++j) inv(i, j) += e.vectors(i, k) * e.vectors(j, k) / e.values[k];
    }
    if (rank) *rank = r;
    return inv;
}

std::vector<double> internal_gradient(const Matrix& B, const Matrix& Ginv, const std::vector<double>& gx) {
    const std::size_t nq = B.rows(), n3 = B.cols();
    std::vector<double> Bg(nq, 0.0), gq(nq, 0.0);
    for (std::size_t i = 0; i < nq; ++i)
        for (std::size_t x = 0; x < n3; ++x) Bg[i] += B(i, x) * gx[x];
    for (std::size_t i = 0; i < nq; ++i)
        for (std::size_t j = 0; j < nq; ++j) gq[i] += Ginv(i, j) * Bg[j];
    return gq;
}

BackTransform back_transform(const InternalCoordinates& ic, const Molecule& mol, const std::vector<double>& dq,
                             int max_iter) {
    BackTransform r;
    const std::vector<double> q0 = ic.values(mol);
    std::vector<double> target(q0.size());
    for (std::size_t i = 0; i < q0.size(); ++i) target[i] = q0[i] + dq[i];
    Molecule x = mol;
    const std::size_t n3 = 3 * mol.atoms.size(), nq = q0.size();
    Molecule first;  // first-order step, fallback
    double prev = 1e300;
    for (int it = 1; it <= max_iter; ++it) {
        const std::vector<double> q = ic.values(x);
        const std::vector<double> d = ic.difference(target, q);
        const Matrix B = ic.wilson_b(x);
        Matrix G(nq, nq);
        for (std::size_t i = 0; i < nq; ++i)
            for (std::size_t j = 0; j < nq; ++j) {
                double s = 0.0;
                for (std::size_t k = 0; k < n3; ++k) s += B(i, k) * B(j, k);
                G(i, j) = s;
            }
        const Matrix Gi = generalized_inverse(G);
        std::vector<double> t(nq, 0.0);
        for (std::size_t i = 0; i < nq; ++i)
            for (std::size_t j = 0; j < nq; ++j) t[i] += Gi(i, j) * d[j];
        double dmax = 0.0;
        for (std::size_t k = 0; k < n3; ++k) {
            double s = 0.0;
            for (std::size_t i = 0; i < nq; ++i) s += B(i, k) * t[i];
            x.atoms[k / 3].position[k % 3] += s;
            dmax = std::max(dmax, std::fabs(s));
        }
        r.iterations = it;
        if (it == 1) first = x;
        if (dmax > 10.0 * prev && it > 2) {  // diverging: keep the first-order step
            x = first;
            break;
        }
        prev = dmax;
        if (dmax < 1.0e-10) {
            r.converged = true;
            break;
        }
    }
    const std::vector<double> qf = ic.values(x);
    r.residual = 0.0;
    for (std::size_t i = 0; i < nq; ++i) r.residual = std::max(r.residual, std::fabs(qf[i] - target[i]));
    r.molecule = x;
    return r;
}

}  // namespace macqc
