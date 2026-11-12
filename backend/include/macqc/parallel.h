// Minimal thread helper: runs body(thread_index) on `threads` std::threads and
// rethrows the first exception after all have finished.
#pragma once

#include <exception>
#include <mutex>
#include <thread>
#include <vector>

#include "macqc/job.h"

namespace macqc {

template <class Body>
void run_parallel(int threads, Body body) {
    if (threads <= 1) {
        body(0);
        return;
    }
    std::vector<std::thread> pool;
    std::exception_ptr error;
    std::mutex mutex;
    Job* job = current_job();  // workers see the caller's job (cancellation, progress)
    for (int t = 0; t < threads; ++t)
        pool.emplace_back([&, t, job] {
            JobScope scope(job);
            try {
                body(t);
            } catch (...) {
                std::lock_guard<std::mutex> lock(mutex);
                if (!error) error = std::current_exception();
            }
        });
    for (auto& th : pool) th.join();
    if (error) std::rethrow_exception(error);
}

}  // namespace macqc
