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

void report_progress_text(const std::string& text) {
    if (Job* j = g_current_job) j->set_progress_text(text);
}

void Job::append_output(const std::string& text) {
    std::lock_guard<std::mutex> lock(mutex_);
    output_ += text;
}

std::string Job::output_since(std::size_t offset, std::size_t* total) const {
    std::lock_guard<std::mutex> lock(mutex_);
    if (total) *total = output_.size();
    if (offset >= output_.size()) return "";
    return output_.substr(offset);
}

void Job::set_progress_text(const std::string& text) {
    std::lock_guard<std::mutex> lock(mutex_);
    progress_text_ = text;
}

std::string Job::progress_text() const {
    std::lock_guard<std::mutex> lock(mutex_);
    return progress_text_;
}

}  // namespace macqc
