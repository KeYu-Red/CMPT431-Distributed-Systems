#include <iostream>
#include <thread>
#include <chrono>
#include <future>
#include <cstdlib>
#include <vector>
#include "common/cxxopts.h"
#include "common/utils.h"
#include "common/get_time.h"

#define MOD_VALUE 10000

#ifdef ONE_LOCK_QUEUE
#include "queues/one_lock_queue.h"
OneLockQueue<int> my_queue;
#elif defined TWO_LOCK_QUEUE
#include "queues/two_lock_queue.h"
TwoLockQueue<int> my_queue;
#elif defined NON_BLOCKING_QUEUE
#include "queues/non_blocking_queue.h"
NonBlockingQueue<int> my_queue;
#else
#include "queues/non_blocking_queue.h"
NonBlockingQueue<int> my_queue;
#endif

struct Producer
{
    int thread_id_;
    std::atomic<int> &stop_signal_;
    long *produced_;
    double *total_time_;

    Producer(int t_thread_id,
                     std::atomic<int> &t_stop_signal) : thread_id_(t_thread_id),
                                                        stop_signal_(t_stop_signal)
    {
        produced_ = new long;
        *produced_ = 0;
        total_time_ = new double;
        *total_time_ = 0;
    }

    void del(){
        delete produced_;
        delete total_time_;
    }

    void operator()()
    {
        timer t;
        t.start();
        long local_produced = 0;
        int value = thread_id_;
        bool flag = false;
        while (stop_signal_.load(std::memory_order_relaxed) == 0)
        {
            value = (value + thread_id_) % MOD_VALUE;
            // Enqueue
            my_queue.enqueue(value);
            local_produced++;
        }
        *produced_ = local_produced;
        *total_time_ = t.stop();
    }

    void print()
    {
        std::cout << thread_id_ << ", " << *produced_ << ", " << *total_time_ << std::endl;
    }
};

struct Consumer
{
    int thread_id_;
    std::atomic<int> &stop_signal_;
    long *consumed_;
    long *failed_dequeues_;
    double *total_time_;

    Consumer(int t_thread_id,
                     std::atomic<int> &t_stop_signal) : thread_id_(t_thread_id),
                                                        stop_signal_(t_stop_signal)
    {
        consumed_ = new long;
        *consumed_ = 0;
        failed_dequeues_ = new long;
        *failed_dequeues_ = 0;
        total_time_ = new double;
        *total_time_ = 0;
    }

    void del(){
        delete consumed_;
        delete failed_dequeues_;
        delete total_time_;
    }

    void operator()()
    {
        timer t;
        t.start();
        long local_consumed = 0;
        long local_failed_dequeues = 0;
        int value = thread_id_;
        bool flag = false;
        while (stop_signal_.load(std::memory_order_relaxed) == 0)
        {
            // Dequeue
            flag = my_queue.dequeue(&value);
            if (flag)
                local_consumed++;
            else
                local_failed_dequeues++;
        }
        *consumed_ = local_consumed;
        *failed_dequeues_ = local_failed_dequeues;
        *total_time_ = t.stop();
    }

    void print()
    {
        std::cout << thread_id_ << ", " << *consumed_ << ", " << *failed_dequeues_ << ", " << *total_time_ << std::endl;
    }
};

int main(int argc, char *argv[])
{
    cxxopts::Options options("Queue", "Test");
    options.add_options("custom", {
                                      {"n_producers", "Number of producers", cxxopts::value<uint>()->default_value(DEFAULT_NUMBER_OF_PRODUCERS)},
                                      {"n_consumers", "Number of producers", cxxopts::value<uint>()->default_value(DEFAULT_NUMBER_OF_CONSUMERS)},
                                      {"init_allocator", "Number of nodes to pre-allocate", cxxopts::value<long>()->default_value(DEFAULT_INIT_ALLOCATOR)},
                                      {"seconds", "Number of seconds to run the experiment", cxxopts::value<uint>()->default_value(DEFAULT_SECONDS)},
                                  });

    auto cl_options = options.parse(argc, argv);
    uint n_producers = cl_options["n_producers"].as<uint>();
    uint n_consumers = cl_options["n_consumers"].as<uint>();
    uint seconds = cl_options["seconds"].as<uint>();
    long init_allocator = cl_options["init_allocator"].as<long>();
    // serialTest(5);
    // exit(1);

    std::cout << "Pre-allocate with " << init_allocator << " elements\n";
    my_queue.initQueue(init_allocator);
    std::cout << "n_producers = " << n_producers << "\n";
    std::cout << "n_consumers = " << n_consumers << "\n";

    // Setup Producers -------------------------------------------
    std::vector<std::thread> producers;
    std::vector<Producer> producer_data;
    std::atomic<int> stop_signal;
    stop_signal.store(0);
    for (int i = 0; i < n_producers; ++i)
    {
        producer_data.push_back(Producer(i, std::ref(stop_signal)));
        producers.push_back(std::thread(producer_data[i]));
    }

    // Setup Consumers -------------------------------------------
    std::vector<std::thread> consumers;
    std::vector<Consumer> consumer_data;
    for (int i = 0; i < n_consumers; ++i)
    {
        consumer_data.push_back(Consumer(i, std::ref(stop_signal)));
        consumers.push_back(std::thread(consumer_data[i]));
    }

    // Wait for specified time and send stop signal ----------
    std::this_thread::sleep_for(std::chrono::seconds(seconds));
    stop_signal.store(1);

    // Wait for producers and consumers ----------------------------
    for (int i = 0; i < n_producers; ++i)
    {
        producers[i].join();
    }
    for (int i = 0; i < n_consumers; ++i)
    {
        consumers[i].join();
    }

    // Producer stats ----------------------------
    std::cout << "Producer data\n";
    uint64_t total_produced = 0;
    std::cout << "thread_id, produced, total_time\n";
    for (int i = 0; i < n_producers; ++i)
    {
        total_produced += *(producer_data[i].produced_);
        producer_data[i].print();
        producer_data[i].del();
    }
    // Consumer stats ----------------------------
    std::cout << "Consumer data\n";
    uint64_t total_consumed = 0;
    uint64_t total_failed = 0;
    std::cout << "thread_id, consumed, failed, total_time\n";
    for (int i = 0; i < n_consumers; ++i)
    {
        total_consumed += *(consumer_data[i].consumed_);
        total_failed += *(consumer_data[i].failed_dequeues_);
        consumer_data[i].print();
        consumer_data[i].del();
    }

    std::cout << "Total produced = " << total_produced << "\n";
    std::cout << "Total consumed = " << total_consumed << "\n";
    std::cout << "Total failed = " << total_failed << "\n";
    std::cout << "Total throughput = " << (total_produced+total_consumed)/seconds << "\n";
    // free
    my_queue.cleanup();

    return 0;
}
