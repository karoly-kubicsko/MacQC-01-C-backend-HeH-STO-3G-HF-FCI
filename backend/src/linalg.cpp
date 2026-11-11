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

}  // namespace macqc
