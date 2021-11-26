#include "core/graph.h"
#include "core/utils.h"
#include <iomanip>
#include <iostream>
#include <stdlib.h>
#include <thread>
#include <vector>
#include <atomic>

#ifdef USE_INT
#define INIT_PAGE_RANK 100000
#define EPSILON 1000
#define PAGE_RANK(x) (15000 + (5 * x) / 6)
#define CHANGE_IN_PAGE_RANK(x, y) std::abs(x - y)
typedef int64_t PageRankType;
#else
#define INIT_PAGE_RANK 1.0
#define EPSILON 0.01
#define DAMPING 0.85
#define PAGE_RANK(x) (1 - DAMPING + DAMPING * x)
#define CHANGE_IN_PAGE_RANK(x, y) std::fabs(x - y)
typedef float PageRankType;
#endif

std::mutex first_mutex;

static std::atomic<PageRankType>& operator+= (std::atomic<PageRankType>& atomicFloat, PageRankType increment)
{
    PageRankType oldValue;
    PageRankType newValue;
    
    do
    {
        oldValue = atomicFloat.load (std::memory_order_relaxed);
        newValue = oldValue + increment;
    } while (! atomicFloat.compare_exchange_weak (oldValue, newValue,
                                                  std::memory_order_release,
                                                  std::memory_order_relaxed));

    return atomicFloat;
}

void singleThread(uintV start, uintV end, PageRankType *pr_curr,  std::atomic<PageRankType> *pr_next, Graph &g, CustomBarrier *barrier, int thread_id, int max_iters)
{
    timer t1;
    double time_taken = 0.0;
    t1.start();

  for(int iter = 0; iter < max_iters; iter++)
  {
    for (uintV u = start; u < end; u++) {
      uintE out_degree = g.vertices_[u].getOutDegree();
      
      for (uintE i = 0; i < out_degree; i++) {
        uintV v = g.vertices_[u].getOutNeighbor(i);
        // Atomic Operation
        pr_next[v] += pr_curr[u] / out_degree;
      }
    }

    barrier->wait();
    for (uintV v = start; v <end ; v++) {
        PageRankType temp = pr_next[v];
        pr_curr[v] = PAGE_RANK(temp);
        pr_next[v] = (PageRankType)0;
      }
    barrier->wait();
  }
  time_taken = t1.stop();
  
  printf("%d, %lf\n", thread_id, time_taken);
    
}

void pageRankSerial(Graph &g, int max_iters, uint n_threads) {
  uintV n = g.n_;
  uintV node_per_thread = n/n_threads;

  PageRankType *pr_curr = new PageRankType[n];
  std::atomic<PageRankType> *pr_next = new std::atomic<PageRankType>[n];

  for (uintV i = 0; i < n; i++) {
    pr_curr[i] = INIT_PAGE_RANK;
    pr_next[i] = (PageRankType)0;
  }

  // Push based pagerank
  timer t1;
  double time_taken = 0.0;
  // Create threads and distribute the work across T threads
  // -------------------------------------------------------------------
  std::vector<std::thread> thread_vector;
  CustomBarrier my_barrier(n_threads);
  std::cout << "thread_id, time_taken\n";
  // for each vertex 'u', process all its outNeighbors 'v'
  for(int i=0; i<n_threads; i++)
  {
      if(i==n_threads-1)
      {
          thread_vector.push_back(std::thread(singleThread, node_per_thread*i, n, pr_curr, pr_next, std::ref(g), &my_barrier, i, max_iters));
      }else{
          thread_vector.push_back(std::thread(singleThread, node_per_thread*i, node_per_thread*(i+1), pr_curr, pr_next, std::ref(g), &my_barrier, i, max_iters));
      }
  }

  t1.start();
  for(int i=0; i<n_threads; i++)
  {
      thread_vector[i].join();
  }
    //-----------------------//

  time_taken = t1.stop();

  PageRankType sum_of_page_ranks = 0;
  for (uintV u = 0; u < n; u++) {
    sum_of_page_ranks += pr_curr[u];
  }
  std::cout << "Sum of page rank : " << sum_of_page_ranks << "\n";
  std::cout << "Time taken (in seconds) : " << time_taken << "\n";
  delete[] pr_curr;
  delete[] pr_next;
}

int main(int argc, char *argv[]) {
  cxxopts::Options options(
      "page_rank_push",
      "Calculate page_rank using serial and parallel execution");
  options.add_options(
      "",
      {
          {"nWorkers", "Number of workers",
           cxxopts::value<uint>()->default_value(DEFAULT_NUMBER_OF_WORKERS)},
          {"nIterations", "Maximum number of iterations",
           cxxopts::value<uint>()->default_value(DEFAULT_MAX_ITER)},
          {"inputFile", "Input graph file path",
           cxxopts::value<std::string>()->default_value(
               "/scratch/input_graphs/roadNet-CA")},
      });

  auto cl_options = options.parse(argc, argv);
  uint n_workers = cl_options["nWorkers"].as<uint>();
  uint max_iterations = cl_options["nIterations"].as<uint>();
  std::string input_file_path = cl_options["inputFile"].as<std::string>();

#ifdef USE_INT
  std::cout << "Using INT\n";
#else
  std::cout << "Using FLOAT\n";
#endif
  std::cout << std::fixed;
  std::cout << "Number of workers : " << n_workers << "\n";

  Graph g;
  std::cout << "Reading graph\n";
  g.readGraphFromBinary<int>(input_file_path);
  std::cout << "Created graph\n";
  pageRankSerial(g, max_iterations, n_workers);

  return 0;
}
