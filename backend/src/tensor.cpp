#include "macqc/tensor.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>

#include "macqc/gto_integrals.h"
#include "macqc/linalg.h"

namespace macqc {

namespace {

int g_threads = 1;

std::size_t product(const std::vector<std::size_t>& d) {
    std::size_t p = 1;
    for (std::size_t x : d) p *= x;
    return p;
}

[[noreturn]] void bad(const std::string& spec, const std::string& why) {
    throw std::invalid_argument("einsum '" + spec + "': " + why);
}

// Generic permutation: out[idx] = in[...] with out.dims[k] = in.dims[perm[k]].
void permute(const std::vector<std::size_t>& dims, const double* in, const std::vector<int>& perm, double* out) {
    const std::size_t r = dims.size();
    if (r == 0) {
        out[0] = in[0];
        return;
    }
    std::vector<std::size_t> istride(r, 1);
    for (std::size_t k = r - 1; k-- > 0;) istride[k] = istride[k + 1] * dims[k + 1];
    std::vector<std::size_t> odims(r), ostride_in(r);
    for (std::size_t k = 0; k < r; ++k) {
        odims[k] = dims[static_cast<std::size_t>(perm[k])];
        ostride_in[k] = istride[static_cast<std::size_t>(perm[k])];
    }
    const std::size_t total = product(dims);
    if (total == 0) return;
    // Iterate over the output; the last output dimension is the inner loop.
    const std::size_t inner = odims[r - 1], step = ostride_in[r - 1];
    std::vector<std::size_t> idx(r, 0);
    std::size_t pos = 0, base = 0;
    while (pos < total) {
        const double* src = in + base;
        if (step == 1)
            std::copy(src, src + inner, out + pos);
        else
            for (std::size_t i = 0; i < inner; ++i) out[pos + i] = src[i * step];
        pos += inner;
        // increment the multi-index over dimensions 0 .. r-2
        std::size_t k = r - 1;
        while (k-- > 0) {
            ++idx[k];
            base += ostride_in[k];
            if (idx[k] < odims[k]) break;
            base -= idx[k] * ostride_in[k];
            idx[k] = 0;
        }
    }
}

struct Spec {
    std::vector<std::string> in;
    std::string out;
};

Spec parse(const std::string& spec, std::size_t operands) {
    const auto arrow = spec.find("->");
    if (arrow == std::string::npos) bad(spec, "missing '->'");
    Spec s;
    s.out = spec.substr(arrow + 2);
    std::string left = spec.substr(0, arrow);
    std::size_t start = 0;
    while (true) {
        const auto comma = left.find(',', start);
        s.in.push_back(left.substr(start, comma == std::string::npos ? std::string::npos : comma - start));
        if (comma == std::string::npos) break;
        start = comma + 1;
    }
    if (s.in.size() != operands) bad(spec, "wrong number of operands");
    for (std::size_t i = 0; i < s.out.size(); ++i)
        if (s.out.find(s.out[i]) != i) bad(spec, "repeated index in the result");
    return s;
}

// Unary einsum: diagonal over repeated letters, sum over letters not in `out`.
Tensor unary(const std::string& spec, const std::string& in, const std::string& out, const Tensor& a, double alpha) {
    if (in.size() != a.rank()) bad(spec, "operand rank does not match '" + in + "'");
    // letter sizes
    std::vector<std::size_t> odims;
    for (char c : out) {
        const auto p = in.find(c);
        if (p == std::string::npos) bad(spec, std::string("index '") + c + "' not in the operand");
        odims.push_back(a.dim(p));
    }
    for (std::size_t i = 0; i < in.size(); ++i)
        for (std::size_t j = 0; j < i; ++j)
            if (in[i] == in[j] && a.dim(i) != a.dim(j)) bad(spec, "diagonal of unequal dimensions");
    // fast path: pure transpose
    bool unique = true;
    for (std::size_t i = 0; i < in.size(); ++i)
        if (in.find(in[i]) != i) unique = false;
    if (unique && in.size() == out.size()) {
        std::vector<int> perm;
        for (char c : out) perm.push_back(static_cast<int>(in.find(c)));
        Tensor r(odims);
        permute(a.dims(), a.data(), perm, r.data());
        if (alpha != 1.0) r *= alpha;
        return r;
    }
    // general path: loop over the unique letters of the operand
    std::string letters;
    for (char c : in)
        if (letters.find(c) == std::string::npos) letters += c;
    std::vector<std::size_t> ldim(letters.size());
    for (std::size_t k = 0; k < letters.size(); ++k) ldim[k] = a.dim(in.find(letters[k]));
    std::vector<std::size_t> astride(in.size(), 1);
    for (std::size_t k = in.size(); k-- > 1;) astride[k - 1] = astride[k] * a.dim(k);
    std::vector<std::size_t> lstride_a(letters.size(), 0), lstride_o(letters.size(), 0);
    for (std::size_t i = 0; i < in.size(); ++i) lstride_a[letters.find(in[i])] += astride[i];
    std::vector<std::size_t> ostride(out.size(), 1);
    for (std::size_t k = out.size(); k-- > 1;) ostride[k - 1] = ostride[k] * odims[k];
    for (std::size_t k = 0; k < out.size(); ++k) lstride_o[letters.find(out[k])] = ostride[k];
    Tensor r(odims);
    const std::size_t total = product(ldim);
    std::vector<std::size_t> idx(letters.size(), 0);
    for (std::size_t n = 0; n < total; ++n) {
        std::size_t ia = 0, io = 0;
        for (std::size_t k = 0; k < letters.size(); ++k) {
            ia += idx[k] * lstride_a[k];
            io += idx[k] * lstride_o[k];
        }
        r[io] += alpha * a[ia];
        for (std::size_t k = letters.size(); k-- > 0;) {
            if (++idx[k] < ldim[k]) break;
            idx[k] = 0;
        }
    }
    return r;
}

// Reduces an operand to the letters needed (by the other operand or the result).
Tensor reduce_operand(const std::string& spec, const std::string& in, const std::string& keep_from, const Tensor& a,
                      std::string& reduced) {
    reduced.clear();
    for (char c : in)
        if (reduced.find(c) == std::string::npos && keep_from.find(c) != std::string::npos) reduced += c;
    bool same = reduced == in;
    if (same) return a;
    return unary(spec, in, reduced, a, 1.0);
}

}  // namespace

Tensor::Tensor(std::vector<std::size_t> dims, double value) : dims_(std::move(dims)) {
    data_.assign(product(dims_), value);
}

Tensor Tensor::transpose(std::initializer_list<int> perm) const { return transpose(std::vector<int>(perm)); }

Tensor Tensor::transpose(const std::vector<int>& perm) const {
    if (perm.size() != dims_.size()) throw std::invalid_argument("Tensor::transpose: wrong permutation length");
    std::vector<std::size_t> od;
    for (int p : perm) od.push_back(dims_[static_cast<std::size_t>(p)]);
    Tensor r(od);
    permute(dims_, data_.data(), perm, r.data());
    return r;
}

Tensor& Tensor::operator+=(const Tensor& o) { return axpy(1.0, o); }
Tensor& Tensor::operator-=(const Tensor& o) { return axpy(-1.0, o); }
Tensor& Tensor::operator*=(double f) {
    for (double& x : data_) x *= f;
    return *this;
}
Tensor Tensor::operator+(const Tensor& o) const {
    Tensor r = *this;
    return r.axpy(1.0, o);
}
Tensor Tensor::operator-(const Tensor& o) const {
    Tensor r = *this;
    return r.axpy(-1.0, o);
}
Tensor Tensor::operator*(double f) const {
    Tensor r = *this;
    return r *= f;
}
Tensor& Tensor::axpy(double f, const Tensor& o) {
    if (dims_.empty() && data_.empty()) {
        *this = o;
        return *this *= f;
    }
    if (o.dims_ != dims_) throw std::invalid_argument("Tensor: shapes differ");
    for (std::size_t i = 0; i < data_.size(); ++i) data_[i] += f * o.data_[i];
    return *this;
}
void Tensor::fill(double v) { std::fill(data_.begin(), data_.end(), v); }
double Tensor::dot(const Tensor& o) const {
    if (o.data_.size() != data_.size()) throw std::invalid_argument("Tensor::dot: sizes differ");
    double s = 0.0;
    for (std::size_t i = 0; i < data_.size(); ++i) s += data_[i] * o.data_[i];
    return s;
}
double Tensor::norm() const { return std::sqrt(dot(*this)); }

void set_tensor_threads(int threads) { g_threads = resolve_threads(threads); }

Tensor einsum(const std::string& spec, const Tensor& a, double alpha) {
    const Spec s = parse(spec, 1);
    return unary(spec, s.in[0], s.out, a, alpha);
}

void einsum_acc(const std::string& spec, const Tensor& a0, const Tensor& b0, Tensor& c, double alpha) {
    const Spec s = parse(spec, 2);
    if (s.in[0].size() != a0.rank() || s.in[1].size() != b0.rank()) bad(spec, "operand rank does not match");
    std::string sa, sb;
    const Tensor a = reduce_operand(spec, s.in[0], s.in[1] + s.out, a0, sa);
    const Tensor b = reduce_operand(spec, s.in[1], s.in[0] + s.out, b0, sb);
    auto dim_of = [&](char ch) -> std::size_t {
        const auto pa = sa.find(ch);
        if (pa != std::string::npos) return a.dim(pa);
        const auto pb = sb.find(ch);
        if (pb != std::string::npos) return b.dim(pb);
        bad(spec, std::string("index '") + ch + "' of the result not in the operands");
    };
    for (std::size_t i = 0; i < sa.size(); ++i) {
        const auto pb = sb.find(sa[i]);
        if (pb != std::string::npos && b.dim(pb) != a.dim(i)) bad(spec, std::string("dimension mismatch of '") + sa[i] + "'");
    }
    std::string batch, contract, freeA, freeB;
    for (char ch : sa) {
        const bool inB = sb.find(ch) != std::string::npos, inO = s.out.find(ch) != std::string::npos;
        if (inB && inO) batch += ch;
        else if (inB) contract += ch;
        else freeA += ch;
    }
    for (char ch : sb)
        if (sa.find(ch) == std::string::npos) freeB += ch;
    auto size_of = [&](const std::string& l) {
        std::size_t p = 1;
        for (char ch : l) p *= dim_of(ch);
        return p;
    };
    const std::size_t nb = size_of(batch), m = size_of(freeA), n = size_of(freeB), k = size_of(contract);
    // permuted operands
    std::vector<int> pa, pb;
    for (char ch : batch + freeA + contract) pa.push_back(static_cast<int>(sa.find(ch)));
    for (char ch : batch + contract + freeB) pb.push_back(static_cast<int>(sb.find(ch)));
    std::vector<double> A(a.size()), B(b.size());
    permute(a.dims(), a.data(), pa, A.data());
    permute(b.dims(), b.data(), pb, B.data());
    std::vector<double> C(nb * m * n, 0.0);
    for (std::size_t q = 0; q < nb; ++q)
        if (m > 0 && n > 0 && k > 0)
            gemm(false, false, m, n, k, alpha, A.data() + q * m * k, k, B.data() + q * k * n, n, 0.0,
                 C.data() + q * m * n, n, g_threads);
    // C has letters batch + freeA + freeB; permute to the result order
    const std::string cl = batch + freeA + freeB;
    std::vector<std::size_t> cdims;
    for (char ch : cl) cdims.push_back(dim_of(ch));
    std::vector<int> pc;
    std::vector<std::size_t> odims;
    for (char ch : s.out) {
        const auto p = cl.find(ch);
        if (p == std::string::npos) bad(spec, std::string("index '") + ch + "' of the result not in the operands");
        pc.push_back(static_cast<int>(p));
        odims.push_back(dim_of(ch));
    }
    if (cl.size() != s.out.size()) bad(spec, "internal: result letters");
    if (c.rank() == 0 && c.size() == 0) c = Tensor(odims);
    if (c.dims() != odims) bad(spec, "accumulator has the wrong shape");
    std::vector<double> R(C.size());
    permute(cdims, C.data(), pc, R.data());
    for (std::size_t i = 0; i < R.size(); ++i) c[i] += R[i];
}

Tensor einsum(const std::string& spec, const Tensor& a, const Tensor& b, double alpha) {
    Tensor c;
    einsum_acc(spec, a, b, c, alpha);
    return c;
}

Tensor einsum3(const std::string& spec, const Tensor& a, const Tensor& b, const Tensor& c, double alpha) {
    const Spec s = parse(spec, 3);
    // intermediate letters: those of a and b needed by c or the result
    std::string mid;
    for (char ch : s.in[0] + s.in[1])
        if (mid.find(ch) == std::string::npos &&
            (s.in[2].find(ch) != std::string::npos || s.out.find(ch) != std::string::npos))
            mid += ch;
    const Tensor ab = einsum(s.in[0] + "," + s.in[1] + "->" + mid, a, b);
    return einsum(mid + "," + s.in[2] + "->" + s.out, ab, c, alpha);
}

}  // namespace macqc
