// Cancellation and progress of a running calculation (used by the C API's job handle).
//
// The calculation thread installs its Job with JobScope; worker threads started through
// run_parallel (and the RHF Fock builder) inherit it. Long loops call check_cancelled(),
// which throws Cancelled once the job has been cancelled, and report their progress.
#pragma once

#include <atomic>
#include <stdexcept>

namespace macqc {

enum class JobStage : int {
    Idle = 0,
    Integrals = 1,
    SCF = 2,
    Transformation = 3,
    FullCI = 4,
    Done = 5,
    CoupledCluster = 6,  // iteration, total energy, largest amplitude change
    Triples = 7,         // (T) correction; iteration = percent done
    Gradient = 8,        // analytic gradient (iteration = part)
    Hessian = 9,         // analytic Hessian: derivative integrals (iteration = part)
    CPHF = 10,           // coupled-perturbed HF / orbital response (iteration, residual)
    Lambda = 11          // coupled-cluster Lambda equations (iteration, residual)
};

struct Job {
    std::atomic<bool> cancelled{false};
    std::atomic<int> stage{0};
    std::atomic<int> iteration{0};
    std::atomic<double> energy{0.0};    // latest total energy (E_h)
    std::atomic<double> residual{0.0};  // SCF orbital gradient, CI residual norm or CC amplitude change
};

struct Cancelled : std::runtime_error {
    Cancelled() : std::runtime_error("calculation cancelled") {}
};

// The job of the calling thread (nullptr outside a job).
Job* current_job();
void set_current_job(Job* job);

// Throws Cancelled if the current job has been cancelled.
inline void check_cancelled() {
    const Job* j = current_job();
    if (j && j->cancelled.load(std::memory_order_relaxed)) throw Cancelled();
}

void report_progress(JobStage stage, int iteration = 0, double energy = 0.0, double residual = 0.0);

// Installs a job for the lifetime of the scope (restores the previous one).
class JobScope {
public:
    explicit JobScope(Job* job) : previous_(current_job()) { set_current_job(job); }
    ~JobScope() { set_current_job(previous_); }
    JobScope(const JobScope&) = delete;
    JobScope& operator=(const JobScope&) = delete;

private:
    Job* previous_;
};

}  // namespace macqc
