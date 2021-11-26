#include <iostream>
#include <thread>
#include <chrono>
#include <cstdlib>
#include <vector>
#include <fstream>
#include <atomic>
#include "common/cxxopts.h"
#include "common/utils.h"
#include "common/get_time.h"

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

struct ProducerWithInputs
{
    int thread_id_;
    std::atomic<int> *producers_finished_;
    long start_;
    long end_;
    int *inputs_;
    double *total_time_;

    ProducerWithInputs(int t_thread_id,
                       std::atomic<int> *t_producers_finished,
                       long t_start,
                       long t_end,
                       int *t_inputs) : thread_id_(t_thread_id),
                                        producers_finished_(t_producers_finished),
                                        start_(t_start),
                                        end_(t_end),
                                        inputs_(t_inputs)
    {
        total_time_ = new double;
        *total_time_ = 0.0;
    }

    void del()
    {
        delete total_time_;
    }
    long getProcessedElements()
    {
        return end_ - start_;
    }
    void operator()()
    {
        timer t;
        t.start();
        // std::cout << "Start producing..\n";
        for (long i = start_; i < end_; i++)
        {
            my_queue.enqueue(inputs_[i]);
        }

        *total_time_ = t.stop();
        producers_finished_->fetch_add(1);
    }

    void print()
    {
        std::cout << thread_id_ << ", " << end_ - start_ << ", " << *total_time_ << "\n";
    }
};

struct ConsumerWithOutput
{
    int thread_id_;
    std::atomic<int> *producers_finished_;
    int n_producers_;
    int *my_output_;
    long *processed_elements_;
    double *total_time_;

    ConsumerWithOutput(int t_thread_id,
                       std::atomic<int> *t_producers_finished,
                       int t_n_producers,
                       int *t_output) : thread_id_(t_thread_id),
                                        producers_finished_(t_producers_finished),
                                        n_producers_(t_n_producers),
                                        my_output_(t_output)
    {
        processed_elements_ = new long;
        *processed_elements_ = 0;
        total_time_ = new double;
        *total_time_ = 0.0;
    }

    void del()
    {
        delete processed_elements_;
        delete total_time_;
    }
    long getProcessedElements()
    {
        return *processed_elements_;
    }
    int *getOutput()
    {
        return my_output_;
    }
    void operator()()
    {
        timer t;
        t.start();
        int value = 0;
        long count = 0;
        // std::cout << "Start consuming..\n";
        while (true)
        {
            if (my_queue.dequeue(&value))
            {
                my_output_[count] = value;
                count++;
            }
            else if (producers_finished_->load() == n_producers_)
            {
                if (my_queue.dequeue(&value))
                {
                    my_output_[count] = value;
                    count++;
                }
                else
                {
                    break;
                }
            }
        }
        *processed_elements_ = count;
        *total_time_ = t.stop();
    }

    void print()
    {
        std::cout << thread_id_ << ", " << *processed_elements_ << ", " << *total_time_ << "\n";
    }
};

int *readInputFromFile(std::string &inputFilePath, long &n_inputs)
{
    std::ifstream input_file(inputFilePath, std::ifstream::in | std::ios::binary);
    input_file.seekg(0, std::ios::end);
    long size = input_file.tellg();
    input_file.seekg(0);
    // cout << "Size : " << size << "\n";
    char *contents = new char[size];
    input_file.read(contents, size);
    input_file.close();
    n_inputs = (long)size / (long)sizeof(int);
    // cout << "m : " << m << "\n";
    int *inputs = (int *)contents;
    return inputs;
}

int main(int argc, char *argv[])
{
    cxxopts::Options options("Queue", "Test");
    options.add_options("custom", {
                                      {"n_producers", "Number of n_producers", cxxopts::value<uint>()->default_value(DEFAULT_NUMBER_OF_PRODUCERS)},
                                      {"n_consumers", "Number of n_consumers", cxxopts::value<uint>()->default_value(DEFAULT_NUMBER_OF_CONSUMERS)},
                                      {"init_allocator", "Number of nodes to pre-allocate", cxxopts::value<long>()->default_value(DEFAULT_INIT_ALLOCATOR)},
                                      {"input_file", "Input file path for random numbers file", cxxopts::value<std::string>()->default_value("/scratch/assignment5/inputs/rand_10M")},
                                  });

    auto cl_options = options.parse(argc, argv);
    uint n_producers = cl_options["n_producers"].as<uint>();
    uint n_consumers = cl_options["n_consumers"].as<uint>();
    long init_allocator = cl_options["init_allocator"].as<long>();
    std::string input_file_path = cl_options["input_file"].as<std::string>();
    long n;
    int *input_values = readInputFromFile(input_file_path, n);

    std::cout << "Number of elements in test file : " << n << "\n";
    long pre_alloc = std::max(n+1000, init_allocator);
    std::cout << "Pre-allocate with " << pre_alloc << " elements\n";
    my_queue.initQueue(pre_alloc);
    std::cout << "n_producers = " << n_producers << "\n";
    std::cout << "n_consumers = " << n_consumers << "\n";

    // Setup producers -------------------------------------------
    std::thread producers[n_producers];
    std::vector<ProducerWithInputs> producer_data;
    std::atomic<int> producers_finished;
    producers_finished.store(0);
    
    // Setup consumers -------------------------------------------
    std::thread consumers[n_consumers];
    std::vector<ConsumerWithOutput> consumer_data;
    // Create output container for each consumer. For worst case allocate each consumer with n.
    int **output_values;
    output_values = new int*[n_consumers];
    for (int i = 0; i < n_consumers; i++)
    {
        output_values[i] = new int[n];
    }

    // Start producers -------------------------------------------
    long start, end = 0;
    std::cout << "Creating producers and consumers\n";
    timer full_timer;
    full_timer.start();
    for (int i = 0; i < n_producers; ++i)
    {
        start = end;
        end = (i == n_producers - 1) ? n : end + n / n_producers;
        producer_data.push_back(ProducerWithInputs(i, &producers_finished, start, end, input_values));
        producers[i] = std::thread(producer_data[i]);
    }
        
    // Start consumers --------------------------------------------
    // std::cout << "Creating consumers\n";
    for (int i = 0; i < n_consumers; ++i)
    {
        consumer_data.push_back(ConsumerWithOutput(i, &producers_finished, n_producers, output_values[i]));
        consumers[i] = std::thread(consumer_data[i]);
    }

    // Wait for producers and consumers ----------------------------
    for (int i = 0; i < n_producers; ++i)
    {
        producers[i].join();
    }
    for (int i = 0; i < n_consumers; ++i)
    {
        consumers[i].join();
    }
    double total_time = full_timer.stop();

    // Producer and consumer stats ----------------------------
    uint64_t total_produced = 0;
    uint64_t total_consumed = 0;
    std::cout << "Producer data\n";
    for (int i = 0; i < n_producers; ++i)
    {
        total_produced += producer_data[i].getProcessedElements();
        producer_data[i].print();
        producer_data[i].del();
    }
    std::cout << "Consumer data\n";
    for (int i = 0; i < n_consumers; ++i)
    {
        total_consumed += consumer_data[i].getProcessedElements();
        consumer_data[i].print();
    }

    // Verify correctness -------------------------------
    // Validation 1 -------------------------------------
    std::cout << "Verifying correctness\n";
    bool verificationSuccess = true;
    if (total_consumed != total_produced)
    {
        verificationSuccess = false;
        std::cout << "ERROR (Consumer stats) : Cumulative number of elements dequeued (" << total_consumed << ") is not equal to number of elements enqueued(" << n << ")\n";
    }

    // Validation 2 -------------------------------------
    if(verificationSuccess){
        // Check if the enqueued values and dequeud values are the same
        int *combined_output = new int[n];
        long curr_index = 0;
        for (int i = 0; i < n_consumers; i++)
        {
            long curr_processed = consumer_data[i].getProcessedElements();
            int *curr_output = consumer_data[i].getOutput();
            for (long j = 0; j < curr_processed; j++)
            {
                combined_output[curr_index] = curr_output[j];
                curr_index++;
            }
        }
        if(checkEqual<int>(input_values, combined_output, n) == false){
            verificationSuccess = false;
            std::cout << "Not equal : Values not inserted\n";
        }
        else{
            std::cout << "Equal : All inserted elements are present in the queue\n";
        }
        delete[] combined_output;
    }

    // Overall stats
    std::cout << "Total time taken = " << total_time << "\n";
    std::cout << "Total produced = " << total_produced << "\n";
    std::cout << "Total consumed = " << total_consumed << "\n";
    if(verificationSuccess){
        std::cout << "Verification successful\n";
    }
    else{
        std::cout << "Verification failed\n";
    }

    // free 
    for(int i = 0; i < n_consumers; i++)
    {
        consumer_data[i].del();
        delete[] output_values[i];
    }
    delete[] input_values;
    delete[] output_values;
    my_queue.cleanup();

    return 0;
}
