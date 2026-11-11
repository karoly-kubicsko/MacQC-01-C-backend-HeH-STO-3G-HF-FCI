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

}  // namespace macqc
