#include "macqc/linalg.h"

#include <algorithm>
#include <cmath>
#include <numeric>
#include <stdexcept>

namespace macqc {

void Matrix::fill(double value) { std::fill(data_.begin(), data_.end(), value); }

Matrix Matrix::transposed() const {
    Matrix t(cols_, rows_);
    for (std::size_t i = 0; i < rows_; ++i)
        for (std::size_t j = 0; j < cols_; ++j) t(j, i) = (*this)(i, j);
    return t;
}

Matrix matmul(const Matrix& a, const Matrix& b) {
    if (a.cols() != b.rows()) throw std::invalid_argument("matmul: dimension mismatch");
    Matrix c(a.rows(), b.cols());
    for (std::size_t i = 0; i < a.rows(); ++i)
        for (std::size_t j = 0; j < b.cols(); ++j) {
            double sum = 0.0;
            for (std::size_t k = 0; k < a.cols(); ++k) sum += a(i, k) * b(k, j);
            c(i, j) = sum;
        }
    return c;
}

EigenResult jacobi_eigen(const Matrix& symmetric) {
    const std::size_t n = symmetric.rows();
    if (symmetric.cols() != n) throw std::invalid_argument("jacobi_eigen: matrix not square");
    Matrix a = symmetric;
    Matrix v(n, n);
    for (std::size_t i = 0; i < n; ++i) v(i, i) = 1.0;

    for (int sweep = 0; sweep < 100; ++sweep) {
        double off = 0.0, diag = 0.0;
        for (std::size_t p = 0; p < n; ++p) {
            diag += a(p, p) * a(p, p);
            for (std::size_t q = p + 1; q < n; ++q) off += a(p, q) * a(p, q);
        }
        if (off <= 1.0e-32 * (diag + off) || off == 0.0) break;

        for (std::size_t p = 0; p + 1 < n; ++p) {
            for (std::size_t q = p + 1; q < n; ++q) {
                if (a(p, q) == 0.0) continue;
                const double theta = (a(q, q) - a(p, p)) / (2.0 * a(p, q));
                const double t = (theta >= 0.0 ? 1.0 : -1.0) /
                                 (std::fabs(theta) + std::sqrt(theta * theta + 1.0));
                const double c = 1.0 / std::sqrt(t * t + 1.0);
                const double s = t * c;
                for (std::size_t k = 0; k < n; ++k) {  // A <- A J (columns p, q)
                    const double akp = a(k, p), akq = a(k, q);
                    a(k, p) = c * akp - s * akq;
                    a(k, q) = s * akp + c * akq;
                }
                for (std::size_t k = 0; k < n; ++k) {  // A <- J^T A (rows p, q)
                    const double apk = a(p, k), aqk = a(q, k);
                    a(p, k) = c * apk - s * aqk;
                    a(q, k) = s * apk + c * aqk;
                }
                for (std::size_t k = 0; k < n; ++k) {  // V <- V J
                    const double vkp = v(k, p), vkq = v(k, q);
                    v(k, p) = c * vkp - s * vkq;
                    v(k, q) = s * vkp + c * vkq;
                }
            }
        }
    }

    std::vector<std::size_t> order(n);
    std::iota(order.begin(), order.end(), std::size_t{0});
    std::stable_sort(order.begin(), order.end(),
                     [&](std::size_t x, std::size_t y) { return a(x, x) < a(y, y); });

    EigenResult result;
    result.values.resize(n);
    result.vectors = Matrix(n, n);
    for (std::size_t k = 0; k < n; ++k) {
        result.values[k] = a(order[k], order[k]);
        for (std::size_t i = 0; i < n; ++i) result.vectors(i, k) = v(i, order[k]);
    }
    return result;
}

EigenResult symmetric_eigen(const Matrix& symmetric) {
    const std::size_t n = symmetric.rows();
    if (symmetric.cols() != n) throw std::invalid_argument("symmetric_eigen: matrix not square");
    EigenResult result;
    result.values.assign(n, 0.0);
    result.vectors = Matrix(n, n);
    if (n == 0) return result;
    if (n == 1) {
        result.values[0] = symmetric(0, 0);
        result.vectors(0, 0) = 1.0;
        return result;
    }
    // V(i, j) row-major; symmetrised copy of the input.
    std::vector<double> V(n * n);
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t j = 0; j < n; ++j) V[i * n + j] = 0.5 * (symmetric(i, j) + symmetric(j, i));
    auto v = [&](std::size_t i, std::size_t j) -> double& { return V[i * n + j]; };
    std::vector<double> d(n), e(n);

    // ---- tred2: Householder reduction to tridiagonal form.
    for (std::size_t j = 0; j < n; ++j) d[j] = v(n - 1, j);
    for (std::size_t i = n - 1; i > 0; --i) {
        double scale = 0.0, h = 0.0;
        for (std::size_t k = 0; k < i; ++k) scale += std::fabs(d[k]);
        if (scale == 0.0) {
            e[i] = d[i - 1];
            for (std::size_t j = 0; j < i; ++j) {
                d[j] = v(i - 1, j);
                v(i, j) = 0.0;
                v(j, i) = 0.0;
            }
        } else {
            for (std::size_t k = 0; k < i; ++k) {
                d[k] /= scale;
                h += d[k] * d[k];
            }
            double f = d[i - 1];
            double g = std::sqrt(h);
            if (f > 0) g = -g;
            e[i] = scale * g;
            h -= f * g;
            d[i - 1] = f - g;
            for (std::size_t j = 0; j < i; ++j) e[j] = 0.0;
            for (std::size_t j = 0; j < i; ++j) {
                f = d[j];
                v(j, i) = f;
                g = e[j] + v(j, j) * f;
                for (std::size_t k = j + 1; k <= i - 1; ++k) {
                    g += v(k, j) * d[k];
                    e[k] += v(k, j) * f;
                }
                e[j] = g;
            }
            f = 0.0;
            for (std::size_t j = 0; j < i; ++j) {
                e[j] /= h;
                f += e[j] * d[j];
            }
            const double hh = f / (h + h);
            for (std::size_t j = 0; j < i; ++j) e[j] -= hh * d[j];
            for (std::size_t j = 0; j < i; ++j) {
                f = d[j];
                g = e[j];
                for (std::size_t k = j; k <= i - 1; ++k) v(k, j) -= (f * e[k] + g * d[k]);
                d[j] = v(i - 1, j);
                v(i, j) = 0.0;
            }
        }
        d[i] = h;
    }
    for (std::size_t i = 0; i + 1 < n; ++i) {
        v(n - 1, i) = v(i, i);
        v(i, i) = 1.0;
        const double h = d[i + 1];
        if (h != 0.0) {
            for (std::size_t k = 0; k <= i; ++k) d[k] = v(k, i + 1) / h;
            for (std::size_t j = 0; j <= i; ++j) {
                double g = 0.0;
                for (std::size_t k = 0; k <= i; ++k) g += v(k, i + 1) * v(k, j);
                for (std::size_t k = 0; k <= i; ++k) v(k, j) -= g * d[k];
            }
        }
        for (std::size_t k = 0; k <= i; ++k) v(k, i + 1) = 0.0;
    }
    for (std::size_t j = 0; j < n; ++j) {
        d[j] = v(n - 1, j);
        v(n - 1, j) = 0.0;
    }
    v(n - 1, n - 1) = 1.0;
    e[0] = 0.0;

    // ---- tql2 on W = V^T (rotations then act on contiguous rows).
    std::vector<double> W(n * n);
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t j = 0; j < n; ++j) W[j * n + i] = V[i * n + j];
    for (std::size_t i = 1; i < n; ++i) e[i - 1] = e[i];
    e[n - 1] = 0.0;
    double f = 0.0, tst1 = 0.0;
    const double eps = std::ldexp(1.0, -52);
    for (std::size_t l = 0; l < n; ++l) {
        tst1 = std::max(tst1, std::fabs(d[l]) + std::fabs(e[l]));
        std::size_t m = l;
        while (m < n) {
            if (std::fabs(e[m]) <= eps * tst1) break;
            ++m;
        }
        if (m == n) m = n - 1;  // e[n-1] == 0, so this cannot be reached in exact arithmetic
        if (m > l) {
            int iter = 0;
            do {
                if (++iter > 200) throw std::runtime_error("symmetric_eigen: no convergence");
                double g = d[l];
                double p = (d[l + 1] - g) / (2.0 * e[l]);
                double r = std::hypot(p, 1.0);
                if (p < 0) r = -r;
                d[l] = e[l] / (p + r);
                d[l + 1] = e[l] * (p + r);
                const double dl1 = d[l + 1];
                double h = g - d[l];
                for (std::size_t i = l + 2; i < n; ++i) d[i] -= h;
                f += h;
                p = d[m];
                double c = 1.0, c2 = c, c3 = c;
                const double el1 = e[l + 1];
                double s = 0.0, s2 = 0.0;
                for (std::size_t ii = m; ii-- > l;) {
                    c3 = c2;
                    c2 = c;
                    s2 = s;
                    g = c * e[ii];
                    h = c * p;
                    r = std::hypot(p, e[ii]);
                    e[ii + 1] = s * r;
                    s = e[ii] / r;
                    c = p / r;
                    p = c * d[ii] - s * g;
                    d[ii + 1] = h + s * (c * g + s * d[ii]);
                    double* w0 = &W[ii * n];
                    double* w1 = &W[(ii + 1) * n];
                    for (std::size_t k = 0; k < n; ++k) {
                        const double hk = w1[k];
                        w1[k] = s * w0[k] + c * hk;
                        w0[k] = c * w0[k] - s * hk;
                    }
                }
                p = -s * s2 * c3 * el1 * e[l] / dl1;
                e[l] = s * p;
                d[l] = c * p;
                if (!std::isfinite(d[l])) throw std::runtime_error("symmetric_eigen: non-finite values");
            } while (std::fabs(e[l]) > eps * tst1);
        }
        d[l] += f;
        e[l] = 0.0;
    }
    std::vector<std::size_t> order(n);
    std::iota(order.begin(), order.end(), std::size_t{0});
    std::stable_sort(order.begin(), order.end(), [&](std::size_t x, std::size_t y) { return d[x] < d[y]; });
    for (std::size_t k = 0; k < n; ++k) {
        result.values[k] = d[order[k]];
        const double* w = &W[order[k] * n];
        for (std::size_t i = 0; i < n; ++i) result.vectors(i, k) = w[i];
    }
    return result;
}

}  // namespace macqc
