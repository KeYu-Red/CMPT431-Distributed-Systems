#ifndef UTILS_H
#define UTILS_H

#include "cxxopts.h"
#include "get_time.h"
#include <atomic>
#include <condition_variable>
#include <iostream>
#include <limits.h>
#include <mutex>

#define intV int32_t
#define uintV int32_t

#define intE int32_t
#define uintE int32_t

#define DEFAULT_NUMBER_OF_WORKERS "1"
#define DEFAULT_MAX_ITER "10"
#define TIME_PRECISION 5
#define VAL_PRECISION 14
#define THREAD_LOGS 0

struct CustomBarrier {
  int num_of_workers_;
  int current_waiting_;
  int barrier_call_;
  std::mutex my_mutex_;
  std::condition_variable my_cv_;

  CustomBarrier(int t_num_of_workers)
      : num_of_workers_(t_num_of_workers), current_waiting_(0),
        barrier_call_(0) {}

  void wait() {
    std::unique_lock<std::mutex> u_lock(my_mutex_);
    int c = barrier_call_;
    current_waiting_++;
    if (current_waiting_ == num_of_workers_) {
      current_waiting_ = 0;
      // unlock and send signal to wake up
      barrier_call_++;
      u_lock.unlock();
      my_cv_.notify_all();
      return;
    }
    my_cv_.wait(u_lock, [&] { return (c != barrier_call_); });
    //  Condition has been reached. return
  }
};

#endif