// Small runtime-sized dense linear algebra types used throughout MacQC.
#pragma once

#include <cstddef>
#include <vector>

namespace macqc {

// Dense row-major matrix of doubles, sized at runtime.
class Matrix {
public:
    Matrix() = default;
    Matrix(std::size_t rows, std::size_t cols, double value = 0.0)
        : rows_(rows), cols_(cols), data_(rows * cols, value) {}

    std::size_t rows() const { return rows_; }
    std::size_t cols() const { return cols_; }
    double& operator()(std::size_t i, std::size_t j) { return data_[i * cols_ + j]; }
    double operator()(std::size_t i, std::size_t j) const { return data_[i * cols_ + j]; }
    void fill(double value);
    Matrix transposed() const;
    double* data() { return data_.data(); }
    const double* data() const { return data_.data(); }

private:
    std::size_t rows_ = 0, cols_ = 0;
    std::vector<double> data_;
};

// Dense 4-index array with extent n in every dimension, stored row-major.
// Used for two-electron integrals (ij|kl) in chemists' notation.
class Tensor4 {
public:
    Tensor4() = default;
    explicit Tensor4(std::size_t n, double value = 0.0) : n_(n), data_(n * n * n * n, value) {}

    std::size_t extent() const { return n_; }
    double& operator()(std::size_t i, std::size_t j, std::size_t k, std::size_t l) {
        return data_[((i * n_ + j) * n_ + k) * n_ + l];
    }
    double operator()(std::size_t i, std::size_t j, std::size_t k, std::size_t l) const {
        return data_[((i * n_ + j) * n_ + k) * n_ + l];
    }

private:
    std::size_t n_ = 0;
    std::vector<double> data_;
};

// C = A * B with a fixed, sequential summation order.
Matrix matmul(const Matrix& a, const Matrix& b);

// Eigen-decomposition of a real symmetric matrix by cyclic Jacobi rotations.
// Eigenvalues are returned in ascending order; eigenvectors are the columns of
// `vectors`.
struct EigenResult {
    std::vector<double> values;
    Matrix vectors;
};
EigenResult jacobi_eigen(const Matrix& symmetric);

// Eigen-decomposition of a real symmetric matrix by Householder tridiagonalisation
// and the implicit QL algorithm (EISPACK tred2/tql2, as in the public-domain JAMA).
// O(n^3); used for the general molecules. Eigenvalues ascending, vectors in columns.
EigenResult symmetric_eigen(const Matrix& symmetric);

// General matrix product on row-major arrays:
//   C = alpha * op(A) * op(B) + beta * C,  op(X) = X or X^T,
// with op(A) m x k, op(B) k x n, C m x n and leading dimensions lda, ldb, ldc.
// Cache-blocked with packed panels; `threads` > 1 splits the rows of C.
// "Apple Accelerate (BLAS)" or "portable blocked GEMM": how gemm() multiplies (for the output file).
const char* gemm_backend_name();

void gemm(bool trans_a, bool trans_b, std::size_t m, std::size_t n, std::size_t k, double alpha, const double* a,
          std::size_t lda, const double* b, std::size_t ldb, double beta, double* c, std::size_t ldc,
          int threads = 1);

}  // namespace macqc
