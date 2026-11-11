#include "macqc/fci.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <map>
#include <stdexcept>

#include "macqc/report.h"

namespace macqc {

namespace {

// Spin orbital index: s = 2*p + sigma (sigma 0 = alpha, 1 = beta); bit s of a mask.
struct SpinOrbitalIntegrals {
    const Matrix& h;
    const Tensor4& g;  // (pq|rs) chemists', MO basis

    double one(int s, int t) const {
        if ((s & 1) != (t & 1)) return 0.0;
        return h(s >> 1, t >> 1);
    }
    // <st|uv> (physicists')
    double two(int s, int t, int u, int v) const {
        if ((s & 1) != (u & 1) || (t & 1) != (v & 1)) return 0.0;
        return g(s >> 1, u >> 1, t >> 1, v >> 1);
    }
    // <st||uv>
    double anti(int s, int t, int u, int v) const { return two(s, t, u, v) - two(s, t, v, u); }
};

int popcount(unsigned x) {
    int n = 0;
    for (; x; x &= x - 1) ++n;
    return n;
}

// Number of occupied spin orbitals of `mask` strictly between bits i and a.
int occupied_between(unsigned mask, int i, int a) {
    const int lo = std::min(i, a), hi = std::max(i, a);
    int n = 0;
    for (int k = lo + 1; k < hi; ++k)
        if (mask & (1u << k)) ++n;
    return n;
}

// Apply a^dagger_a a_i to a determinant; returns false if the result vanishes.
bool excite(unsigned mask, int a, int i, unsigned& result, int& sign) {
    if (!(mask & (1u << i))) return false;
    unsigned m = mask & ~(1u << i);
    if (m & (1u << a)) return false;
    sign = (occupied_between(m, i, a) % 2 == 0) ? 1 : -1;
    result = m | (1u << a);
    return true;
}

std::vector<int> occupied_list(unsigned mask, int nso) {
    std::vector<int> o;
    for (int k = 0; k < nso; ++k)
        if (mask & (1u << k)) o.push_back(k);
    return o;
}

// <bra|H|ket> by the Slater-Condon rules (electronic Hamiltonian, no nuclear repulsion).
double slater_condon(unsigned bra, unsigned ket, int nso, const SpinOrbitalIntegrals& ints) {
    const unsigned diff = bra ^ ket;
    const int nd = popcount(diff);
    if (nd == 0) {
        const auto occ = occupied_list(ket, nso);
        double e = 0.0;
        for (int i : occ) e += ints.one(i, i);
        for (int i : occ)
            for (int j : occ) e += 0.5 * ints.anti(i, j, i, j);
        return e;
    }
    if (nd == 2) {
        const unsigned from = ket & diff, to = bra & diff;
        int i = 0, a = 0;
        while (!(from & (1u << i))) ++i;
        while (!(to & (1u << a))) ++a;
        unsigned res;
        int sign;
        if (!excite(ket, a, i, res, sign)) return 0.0;
        double v = ints.one(a, i);
        for (int k : occupied_list(bra & ket, nso)) v += ints.anti(a, k, i, k);
        return sign * v;
    }
    if (nd == 4) {
        const unsigned from = ket & diff, to = bra & diff;
        const auto fo = occupied_list(from, nso), tt = occupied_list(to, nso);
        const int i = fo[0], j = fo[1], a = tt[0], b = tt[1];
        unsigned d1, d2;
        int s1, s2;
        if (!excite(ket, a, i, d1, s1)) return 0.0;
        if (!excite(d1, b, j, d2, s2)) return 0.0;
        return s1 * s2 * ints.anti(a, b, i, j);
    }
    return 0.0;
}

// <S^2> operator matrix in the determinant basis via S^2 = S- S+ + Sz(Sz+1).
Matrix spin_squared_matrix(const std::vector<unsigned>& dets, int nmo) {
    const std::size_t nd = dets.size();
    Matrix s2(nd, nd);
    for (std::size_t k = 0; k < nd; ++k) {
        const unsigned ket = dets[k];
        // S+ |ket> : sum_p a^dag_{p alpha} a_{p beta}
        std::map<unsigned, double> plus;
        for (int p = 0; p < nmo; ++p) {
            unsigned r;
            int sg;
            if (excite(ket, 2 * p, 2 * p + 1, r, sg)) plus[r] += sg;
        }
        std::map<unsigned, double> minus;
        for (const auto& [m, c] : plus)
            for (int p = 0; p < nmo; ++p) {
                unsigned r;
                int sg;
                if (excite(m, 2 * p + 1, 2 * p, r, sg)) minus[r] += c * sg;
            }
        int na = 0, nb = 0;
        for (int p = 0; p < nmo; ++p) {
            if (ket & (1u << (2 * p))) ++na;
            if (ket & (1u << (2 * p + 1))) ++nb;
        }
        const double sz = 0.5 * (na - nb);
        minus[ket] += sz * (sz + 1.0);
        for (std::size_t b = 0; b < nd; ++b) {
            auto it = minus.find(dets[b]);
            if (it != minus.end()) s2(b, k) += it->second;
        }
    }
    return s2;
}

std::string spin_label(double s2) {
    const double s = 0.5 * (-1.0 + std::sqrt(1.0 + 4.0 * std::max(s2, 0.0)));
    const int mult = static_cast<int>(std::lround(2.0 * s + 1.0));
    static const char* names[] = {"", "singlet", "doublet", "triplet", "quartet", "quintet",
                                  "sextet"};
    if (mult >= 1 && mult <= 6) return names[mult];
    return "multiplicity " + std::to_string(mult);
}

std::string orbital_name(int mo) { return std::to_string(mo + 1) + "sigma"; }

// ASCII form of a determinant label for the text report (alpha -> A, beta -> B).
std::string ascii_label(const std::string& label) {
    std::string ascii;
    for (std::size_t k = 0; k < label.size();) {
        if (label.compare(k, 2, "\u03b1") == 0) { ascii += 'A'; k += 2; }
        else if (label.compare(k, 2, "\u03b2") == 0) { ascii += 'B'; k += 2; }
        else { ascii += label[k]; ++k; }
    }
    return ascii;
}

}  // namespace

FCIResult run_fci(const Matrix& h_ao, const Tensor4& eri_ao, const Matrix& C,
                  double nuclear_repulsion, double hf_total_energy) {
    const std::size_t nmo = C.cols(), nao = C.rows();
    if (nmo < 1 || nmo > 8) throw std::invalid_argument("run_fci: unsupported number of MOs");
    if (h_ao.rows() != nao || eri_ao.extent() != nao)
        throw std::invalid_argument("run_fci: dimension mismatch");

    FCIResult res;

    // h_pq = sum C_mu,p H_mu,nu C_nu,q
    res.h_mo = matmul(matmul(C.transposed(), h_ao), C);

    // (pq|rs): four successive quarter transformations.
    Tensor4 t1(nao), t2(nao), t3(nao);
    auto transform = [&](const Tensor4& in, Tensor4& out, int index) {
        out = Tensor4(nao);
        for (std::size_t a = 0; a < nmo; ++a)
            for (std::size_t i = 0; i < nao; ++i)
                for (std::size_t j = 0; j < nao; ++j)
                    for (std::size_t k = 0; k < nao; ++k) {
                        double sum = 0.0;
                        for (std::size_t mu = 0; mu < nao; ++mu) {
                            const double c = C(mu, a);
                            switch (index) {
                                case 0: sum += c * in(mu, i, j, k); break;
                                case 1: sum += c * in(i, mu, j, k); break;
                                case 2: sum += c * in(i, j, mu, k); break;
                                default: sum += c * in(i, j, k, mu); break;
                            }
                        }
                        switch (index) {
                            case 0: out(a, i, j, k) = sum; break;
                            case 1: out(i, a, j, k) = sum; break;
                            case 2: out(i, j, a, k) = sum; break;
                            default: out(i, j, k, a) = sum; break;
                        }
                    }
    };
    transform(eri_ao, t1, 0);
    transform(t1, t2, 1);
    transform(t2, t3, 2);
    transform(t3, res.eri_mo, 3);
    // (only the first nmo indices of each dimension are meaningful; nmo == nao here)

    // Determinants: one alpha electron in MO p, one beta electron in MO q.
    std::vector<unsigned> dets;
    for (std::size_t p = 0; p < nmo; ++p)
        for (std::size_t q = 0; q < nmo; ++q) {
            dets.push_back((1u << (2 * p)) | (1u << (2 * q + 1)));
            res.determinants.push_back("|" + std::to_string(p + 1) + "α " +
                                       std::to_string(q + 1) + "β>");
        }
    const std::size_t nd = dets.size();
    const int nso = static_cast<int>(2 * nmo);

    SpinOrbitalIntegrals ints{res.h_mo, res.eri_mo};
    res.hamiltonian = Matrix(nd, nd);
    for (std::size_t i = 0; i < nd; ++i)
        for (std::size_t j = 0; j < nd; ++j)
            res.hamiltonian(i, j) = slater_condon(dets[i], dets[j], nso, ints);
    // enforce exact symmetry
    for (std::size_t i = 0; i < nd; ++i)
        for (std::size_t j = i + 1; j < nd; ++j) {
            const double m = 0.5 * (res.hamiltonian(i, j) + res.hamiltonian(j, i));
            res.hamiltonian(i, j) = res.hamiltonian(j, i) = m;
        }

    const Matrix s2 = spin_squared_matrix(dets, static_cast<int>(nmo));
    const EigenResult eig = jacobi_eigen(res.hamiltonian);

    for (std::size_t k = 0; k < nd; ++k) {
        std::vector<double> c(nd);
        for (std::size_t i = 0; i < nd; ++i) c[i] = eig.vectors(i, k);
        // Sign convention: largest component positive.
        std::size_t big = 0;
        for (std::size_t i = 1; i < nd; ++i)
            if (std::fabs(c[i]) > std::fabs(c[big]) + 1.0e-12) big = i;
        if (c[big] < 0.0)
            for (double& x : c) x = -x;

        double expect = 0.0;
        for (std::size_t i = 0; i < nd; ++i)
            for (std::size_t j = 0; j < nd; ++j) expect += c[i] * s2(i, j) * c[j];
        if (std::fabs(expect) < 1.0e-12) expect = 0.0;

        FCIState st;
        st.energy = eig.values[k] + nuclear_repulsion;
        st.S2 = expect;
        st.spin = spin_label(expect);
        res.states.push_back(st);
        res.vectors.push_back(c);
    }

    res.electronic_energy = eig.values[0];
    res.total_energy = res.states[0].energy;
    res.correlation_energy = res.total_energy - hf_total_energy;
    res.ground_state_vector = res.vectors[0];

    // Weights of the singlet configurations (p,q), p <= q, in the ground state.
    for (std::size_t p = 0; p < nmo; ++p)
        for (std::size_t q = p; q < nmo; ++q) {
            SingletWeight w;
            const double c_pq = res.ground_state_vector[p * nmo + q];
            const double c_qp = res.ground_state_vector[q * nmo + p];
            if (p == q) {
                w.configuration = "(" + orbital_name(static_cast<int>(p)) + ")^2";
                w.weight = c_pq * c_pq;
            } else {
                w.configuration = "(" + orbital_name(static_cast<int>(p)) + ")^1(" +
                                  orbital_name(static_cast<int>(q)) + ")^1";
                w.weight = c_pq * c_pq + c_qp * c_qp;
            }
            res.singlet_weights.push_back(w);
        }
    return res;
}

std::string fci_report(const FCIResult& r, int print_level) {
    std::string out;
    if (print_level == 0) return out;
    const std::size_t nd = r.determinants.size();
    auto upper = [](std::string s) {
        for (char& ch : s) ch = static_cast<char>(std::toupper(static_cast<unsigned char>(ch)));
        return s;
    };

    out += "\n\n\n    FULL CI CALCULATION (TWO ELECTRONS, M_S = 0, " + std::to_string(nd) +
           " DETERMINANTS)\n";
    out += "\n    DETERMINANTS (A = ALPHA, B = BETA; ORBITALS ARE THE HF MOLECULAR ORBITALS)\n\n";
    for (std::size_t i = 0; i < nd; ++i)
        out += "    " + integer(static_cast<long>(i) + 1, 3) + "    " +
               ascii_label(r.determinants[i]) + "\n";

    if (print_level >= 2) {
        matout(out, r.h_mo, "HMO ");
        matout(out, r.hamiltonian, "HCI ");
        out += "\n\n\n";
        const std::size_t n = r.eri_mo.extent();
        for (std::size_t i = 0; i < n; ++i)
            for (std::size_t j = 0; j < n; ++j)
                for (std::size_t k = 0; k < n; ++k)
                    for (std::size_t l = 0; l < n; ++l)
                        out += "   (" + integer(static_cast<long>(i) + 1, 2) +
                               integer(static_cast<long>(j) + 1, 2) +
                               integer(static_cast<long>(k) + 1, 2) +
                               integer(static_cast<long>(l) + 1, 2) + " )" +
                               fixed(r.eri_mo(i, j, k, l), 10, 6) + "   (MO BASIS)\n";
    }

    out += "\n\n\n    FULL CI EIGENVALUES (TOTAL ENERGIES, HARTREE)\n\n";
    out += "    STATE         TOTAL ENERGY          <S**2>   SPIN\n";
    for (std::size_t k = 0; k < r.states.size(); ++k)
        out += "    " + integer(static_cast<long>(k) + 1, 5) + "   " +
               fortran_d(r.states[k].energy, 20, 12) + "   " + fixed(r.states[k].S2, 8, 4) + "   " +
               upper(r.states[k].spin) + "\n";

    out += "\n\n    GROUND STATE CI VECTOR\n\n";
    for (std::size_t i = 0; i < nd; ++i)
        out += "    " + integer(static_cast<long>(i) + 1, 3) + "    " +
               ascii_label(r.determinants[i]) + "   " + fixed(r.ground_state_vector[i], 12, 8) + "\n";

    out += "\n\n    WEIGHTS OF THE SINGLET CONFIGURATIONS IN THE GROUND STATE\n\n";
    for (const auto& w : r.singlet_weights) {
        const std::string& name = w.configuration;
        out += "    " + name + std::string(name.size() < 22 ? 22 - name.size() : 1, ' ') +
               fixed(w.weight, 12, 8) + "\n";
    }

    out += "\n\n    FCI ELECTRONIC ENERGY =   " + fortran_d(r.electronic_energy, 20, 12) + "\n\n";
    out += "    FCI TOTAL ENERGY =        " + fortran_d(r.total_energy, 20, 12) + "\n\n";
    out += "    CORRELATION ENERGY =      " + fortran_d(r.correlation_energy, 20, 12) + "\n";
    return out;
}

}  // namespace macqc
