// Dense tensors of arbitrary rank with an einsum-style contraction (ΨLab 5.0), used for the
// coupled-cluster Lambda equations and density matrices. A contraction of two tensors is
// mapped to a matrix product (gemm) after transposing the operands; indices that appear in
// both operands and the result are batch indices (one product per value).
//
//   Tensor c = einsum("ikca,jkcb->aijb", l2, t2);
//   einsum_acc("ia,jb->ijab", a, b, c, -0.5);      // c += -0.5 * einsum(...)
//
// Single-operand forms transpose and sum ("aijb->iabj", "dljd->jl" with a repeated index =
// diagonal).
#pragma once

#include <cstddef>
#include <initializer_list>
#include <string>
#include <vector>

namespace macqc {

class Tensor {
public:
    Tensor() = default;
    explicit Tensor(std::vector<std::size_t> dims, double value = 0.0);

    const std::vector<std::size_t>& dims() const { return dims_; }
    std::size_t dim(std::size_t k) const { return dims_[k]; }
    std::size_t rank() const { return dims_.size(); }
    std::size_t size() const { return data_.size(); }
    double* data() { return data_.data(); }
    const double* data() const { return data_.data(); }
    std::vector<double>& vec() { return data_; }
    const std::vector<double>& vec() const { return data_; }

    double& operator[](std::size_t i) { return data_[i]; }
    double operator[](std::size_t i) const { return data_[i]; }
    double& operator()(std::size_t i, std::size_t j) { return data_[i * dims_[1] + j]; }
    double operator()(std::size_t i, std::size_t j) const { return data_[i * dims_[1] + j]; }
    double& operator()(std::size_t i, std::size_t j, std::size_t k, std::size_t l) {
        return data_[((i * dims_[1] + j) * dims_[2] + k) * dims_[3] + l];
    }
    double operator()(std::size_t i, std::size_t j, std::size_t k, std::size_t l) const {
        return data_[((i * dims_[1] + j) * dims_[2] + k) * dims_[3] + l];
    }

    // numpy transpose: result.dims[k] = dims[perm[k]].
    Tensor transpose(std::initializer_list<int> perm) const;
    Tensor transpose(const std::vector<int>& perm) const;

    Tensor& operator+=(const Tensor& o);
    Tensor& operator-=(const Tensor& o);
    Tensor& operator*=(double f);
    Tensor operator+(const Tensor& o) const;
    Tensor operator-(const Tensor& o) const;
    Tensor operator*(double f) const;
    // this += f * o
    Tensor& axpy(double f, const Tensor& o);
    void fill(double v);
    double dot(const Tensor& o) const;
    double norm() const;

private:
    std::vector<std::size_t> dims_;
    std::vector<double> data_;
};

// Two-operand contraction "AB,CD->EF" (single letters as indices).
Tensor einsum(const std::string& spec, const Tensor& a, const Tensor& b, double alpha = 1.0);
// c += alpha * einsum(spec, a, b) (c must have the result shape, or be empty).
void einsum_acc(const std::string& spec, const Tensor& a, const Tensor& b, Tensor& c, double alpha = 1.0);
// Single-operand transpose / partial trace / sum.
Tensor einsum(const std::string& spec, const Tensor& a, double alpha = 1.0);
// Three operands: contracted left to right through the intermediate of the first two.
Tensor einsum3(const std::string& spec, const Tensor& a, const Tensor& b, const Tensor& c, double alpha = 1.0);

// Number of threads used by the tensor products (0: all cores).
void set_tensor_threads(int threads);

}  // namespace macqc
