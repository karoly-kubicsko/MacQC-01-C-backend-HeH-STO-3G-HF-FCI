#include "macqc/job.h"

namespace macqc {

namespace {
thread_local Job* g_current_job = nullptr;
}

Job* current_job() { return g_current_job; }

void set_current_job(Job* job) { g_current_job = job; }

void report_progress(JobStage stage, int iteration, double energy, double residual) {
    Job* j = g_current_job;
    if (!j) return;
    j->stage.store(static_cast<int>(stage), std::memory_order_relaxed);
    j->iteration.store(iteration, std::memory_order_relaxed);
    j->energy.store(energy, std::memory_order_relaxed);
    j->residual.store(residual, std::memory_order_relaxed);
}

}  // namespace macqc
