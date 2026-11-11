#include "macqc/scf.h"

#include <cmath>
#include <stdexcept>

#include "macqc/report.h"

namespace macqc {

namespace {
constexpr std::size_t kNbf = 2;  // stage 1: one s function per atom
}

TwoElectronDiatomicHF::TwoElectronDiatomicHF(const DiatomicBasis& basis, const HFOptions& options)
    : basis_(basis),
      opt_(options),
      S_(kNbf, kNbf), X_(kNbf, kNbf), XT_(kNbf, kNbf), H_(kNbf, kNbf), F_(kNbf, kNbf),
      G_(kNbf, kNbf), C_(kNbf, kNbf), Fprime_(kNbf, kNbf), Cprime_(kNbf, kNbf), P_(kNbf, kNbf),
      oldP_(kNbf, kNbf), E_(kNbf, kNbf), TT_(kNbf), PS_(kNbf, kNbf) {
    if (basis.n_gaussians < 1 || basis.n_gaussians > 3)
        throw std::invalid_argument("n_gaussians must be 1, 2 or 3");
}

void TwoElectronDiatomicHF::calculate_integrals() {
    ints_ = macqc::calculate_integrals(basis_, opt_.exact_erf);
    if (opt_.print_level == 0) return;

    const Integrals& v = ints_;
    auto f11 = [](double x) { return fixed(x, 11, 6); };
    out_ += "   STO-" + std::to_string(basis_.n_gaussians) + "G FOR ATOMIC NUMBERS " +
            fixed(basis_.za, 5, 2) + " AND " + fixed(basis_.zb, 5, 2) + "\n\n\n";
    out_ += "   R          ZETA1      ZETA2      S12        T11\n\n";
    out_ += f11(basis_.r) + f11(basis_.zeta1) + f11(basis_.zeta2) + f11(v.S12) + f11(v.T11) + "\n\n\n";
    out_ += "   T12        T22        V11A       V12A       V22A\n\n";
    out_ += f11(v.T12) + f11(v.T22) + f11(v.V11A) + f11(v.V12A) + f11(v.V22A) + "\n\n\n";
    out_ += "   V11B       V12B       V22B       V1111      V2111\n\n";
    out_ += f11(v.V11B) + f11(v.V12B) + f11(v.V22B) + f11(v.V1111) + f11(v.V2111) + "\n\n\n";
    out_ += "   V2121      V2211      V2221      V2222\n\n";
    out_ += f11(v.V2121) + f11(v.V2211) + f11(v.V2221) + f11(v.V2222) + "\n";
}

void TwoElectronDiatomicHF::collect_matrices() {
    const Integrals& v = ints_;
    H_(0, 0) = v.T11 + v.V11A + v.V11B;
    H_(0, 1) = v.T12 + v.V12A + v.V12B;
    H_(1, 0) = H_(0, 1);
    H_(1, 1) = v.T22 + v.V22A + v.V22B;

    S_(0, 0) = 1.0;
    S_(0, 1) = v.S12;
    S_(1, 0) = v.S12;
    S_(1, 1) = 1.0;

    // Symmetric (Lowdin-like) orthogonalisation for two functions.
    X_(0, 0) = 1.0 / std::sqrt(2.0 * (1.0 + v.S12));
    X_(1, 0) = X_(0, 0);
    X_(0, 1) = 1.0 / std::sqrt(2.0 * (1.0 - v.S12));
    X_(1, 1) = -X_(0, 1);
    XT_ = X_.transposed();

    TT_ = Tensor4(kNbf);
    TT_(0, 0, 0, 0) = v.V1111;
    TT_(1, 0, 0, 0) = v.V2111;
    TT_(0, 1, 0, 0) = v.V2111;
    TT_(0, 0, 1, 0) = v.V2111;
    TT_(0, 0, 0, 1) = v.V2111;
    TT_(1, 0, 1, 0) = v.V2121;
    TT_(0, 1, 1, 0) = v.V2121;
    TT_(1, 0, 0, 1) = v.V2121;
    TT_(0, 1, 0, 1) = v.V2121;
    TT_(1, 1, 0, 0) = v.V2211;
    TT_(0, 0, 1, 1) = v.V2211;
    TT_(1, 1, 1, 0) = v.V2221;
    TT_(1, 1, 0, 1) = v.V2221;
    TT_(1, 0, 1, 1) = v.V2221;
    TT_(0, 1, 1, 1) = v.V2221;
    TT_(1, 1, 1, 1) = v.V2222;

    if (opt_.print_level != 0) {
        matout(out_, S_, "S   ");
        matout(out_, X_, "X   ");
        matout(out_, H_, "H   ");
        out_ += "\n\n\n";
        for (std::size_t i = 0; i < kNbf; ++i)
            for (std::size_t j = 0; j < kNbf; ++j)
                for (std::size_t k = 0; k < kNbf; ++k)
                    for (std::size_t l = 0; l < kNbf; ++l)
                        out_ += "   (" + integer(static_cast<long>(i) + 1, 2) +
                                integer(static_cast<long>(j) + 1, 2) +
                                integer(static_cast<long>(k) + 1, 2) +
                                integer(static_cast<long>(l) + 1, 2) + " )" +
                                fixed(TT_(i, j, k, l), 10, 6) + "\n";
    }
}

void TwoElectronDiatomicHF::form_g() {
    G_.fill(0.0);
    for (std::size_t i = 0; i < kNbf; ++i)
        for (std::size_t j = 0; j < kNbf; ++j)
            for (std::size_t k = 0; k < kNbf; ++k)
                for (std::size_t l = 0; l < kNbf; ++l)
                    G_(i, j) += P_(k, l) * (TT_(i, j, k, l) - 0.5 * TT_(i, l, k, j));
}

// Analytic diagonalisation of a real symmetric 2x2 matrix, with the original
// eigenvalue-ordering rule (swap when E2 <= E1).
void TwoElectronDiatomicHF::diag_2x2(const Matrix& f, Matrix& c, Matrix& e) {
    c = Matrix(2, 2);
    e = Matrix(2, 2);
    double theta;
    if (std::fabs(f(0, 0) - f(1, 1)) <= 1.0e-20)
        theta = kPi / 4.0;
    else
        theta = 0.5 * std::atan(2.0 * f(0, 1) / (f(0, 0) - f(1, 1)));

    const double ct = std::cos(theta), st = std::sin(theta);
    c(0, 0) = ct;
    c(1, 0) = st;
    c(0, 1) = st;
    c(1, 1) = -ct;

    e(0, 0) = f(0, 0) * (ct * ct) + f(1, 1) * (st * st) + f(0, 1) * std::sin(2.0 * theta);
    e(1, 1) = f(1, 1) * (ct * ct) + f(0, 0) * (st * st) - f(0, 1) * std::sin(2.0 * theta);
    e(1, 0) = 0.0;
    e(0, 1) = 0.0;

    if (e(1, 1) <= e(0, 0)) {
        const double tmp = e(0, 0);
        e(0, 0) = e(1, 1);
        e(1, 1) = tmp;
        for (std::size_t i = 0; i < 2; ++i) {
            const double t = c(i, 0);
            c(i, 0) = c(i, 1);
            c(i, 1) = t;
        }
    }
}

void TwoElectronDiatomicHF::scf() {
    const int iop = opt_.print_level;
    P_.fill(0.0);
    if (iop >= 2) matout(out_, P_, "P   ");

    converged_ = false;
    iterations_ = 0;
    history_.clear();
    double electronic_energy = 0.0;

    for (int it = 1; it <= opt_.max_iter; ++it) {
        iterations_ = it;
        if (iop >= 2) out_ += "\n    START OF ITERATION NUMBER = " + integer(it, 2) + "\n";

        form_g();
        if (iop >= 2) matout(out_, G_, "G   ");

        for (std::size_t i = 0; i < kNbf; ++i)
            for (std::size_t j = 0; j < kNbf; ++j) F_(i, j) = H_(i, j) + G_(i, j);
        double sum = 0.0;
        for (std::size_t i = 0; i < kNbf; ++i)
            for (std::size_t j = 0; j < kNbf; ++j) sum += P_(i, j) * (H_(i, j) + F_(i, j));
        electronic_energy = 0.5 * sum;

        if (iop >= 2) {
            matout(out_, F_, "F   ");
            out_ += "\n\n\n    ELECTRONIC ENERGY = " + fortran_d(electronic_energy, 20, 12) + "\n";
        }

        const Matrix temp_g = matmul(F_, X_);
        Fprime_ = matmul(XT_, temp_g);
        diag_2x2(Fprime_, Cprime_, E_);
        C_ = matmul(X_, Cprime_);

        oldP_ = P_;
        P_.fill(0.0);
        for (std::size_t i = 0; i < kNbf; ++i)
            for (std::size_t j = 0; j < kNbf; ++j) P_(i, j) = 2.0 * C_(i, 0) * C_(j, 0);

        if (iop >= 2) {
            matout(out_, Fprime_, "F'  ");
            matout(out_, Cprime_, "C'  ");
            matout(out_, E_, "E   ");
            matout(out_, C_, "C   ");
            matout(out_, P_, "P   ");
        }

        double dsum = 0.0;
        for (std::size_t i = 0; i < kNbf; ++i)
            for (std::size_t j = 0; j < kNbf; ++j) {
                const double d = P_(i, j) - oldP_(i, j);
                dsum += d * d;
            }
        const double delta = std::sqrt(dsum / 4.0);
        history_.push_back({it, electronic_energy, delta});
        if (iop != 0)
            out_ += "\n    DELTA(CONVERGENCE OF DENSITY MATRIX) = " + fixed(delta, 10, 6) + "\n\n";

        if (delta < opt_.convergence) {
            converged_ = true;
            break;
        }
    }

    electronic_energy_ = electronic_energy;
    total_energy_ = electronic_energy + basis_.za * basis_.zb / basis_.r;

    if (!converged_) {
        out_ += "    NO CONVERGENCE IN SCF\n";
    } else {
        if (iop != 0) {
            out_ += "\n\n    CALCULATION CONVERGED\n\n";
            out_ += "    ELECTRONIC ENERGY = " + fortran_d(electronic_energy_, 20, 12) + "\n\n";
            out_ += "    TOTAL ENERGY =      " + fortran_d(total_energy_, 20, 12) + "\n";
        }
        if (iop == 1) {
            matout(out_, G_, "G   ");
            matout(out_, F_, "F   ");
            matout(out_, E_, "E   ");
            matout(out_, C_, "C   ");
            matout(out_, P_, "P   ");
        }
    }

    PS_ = matmul(P_, S_);
    if (converged_ && iop != 0) matout(out_, PS_, "PS  ");
}

HFResult TwoElectronDiatomicHF::run() {
    out_.clear();
    calculate_integrals();
    collect_matrices();
    scf();

    HFResult r;
    r.converged = converged_;
    r.iterations = iterations_;
    r.electronic_energy = electronic_energy_;
    r.nuclear_repulsion = basis_.za * basis_.zb / basis_.r;
    r.total_energy = total_energy_;
    for (std::size_t i = 0; i < kNbf; ++i) r.orbital_energies.push_back(E_(i, i));
    r.C = C_;
    r.P = P_;
    r.S = S_;
    r.X = X_;
    r.H = H_;
    r.F = F_;
    r.G = G_;
    r.PS = PS_;
    r.eri = TT_;
    r.integrals = ints_;
    r.history = history_;
    r.report = out_;
    return r;
}

}  // namespace macqc
