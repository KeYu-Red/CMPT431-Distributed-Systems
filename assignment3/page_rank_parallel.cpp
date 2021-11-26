#include "core/graph.h"
#include "core/utils.h"
#include <iomanip>
#include <iostream>
#include <stdlib.h>
#include <thread>
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
void singleThread(uintV start, uintV end, PageRankType *pr_curr,  std::atomic<PageRankType> *pr_next, Graph &g, CustomBarrier *barrier, int thread_id, int max_iters);
void pageRankStrategyOne(Graph &g, int max_iters, uint n_threads);
void pageRankStrategyTwo(Graph &g, int max_iters, uint n_threads);

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

void pageRankSerial(Graph &g, int max_iters) {
  uintV n = g.n_;

  PageRankType *pr_curr = new PageRankType[n];
  PageRankType *pr_next = new PageRankType[n];

  for (uintV i = 0; i < n; i++) {
    pr_curr[i] = INIT_PAGE_RANK;
    pr_next[i] = 0.0;
  }

  // Push based pagerank
  timer t1;
  double time_taken = 0.0;
  // Create threads and distribute the work across T threads
  // -------------------------------------------------------------------
  t1.start();
  for (int iter = 0; iter < max_iters; iter++) {
    // for each vertex 'u', process all its outNeighbors 'v'
    for (uintV u = 0; u < n; u++) {
      uintE out_degree = g.vertices_[u].getOutDegree();
      for (uintE i = 0; i < out_degree; i++) {
        uintV v = g.vertices_[u].getOutNeighbor(i);
        pr_next[v] += (pr_curr[u] / out_degree);
      }
    }
    for (uintV v = 0; v < n; v++) {
      pr_next[v] = PAGE_RANK(pr_next[v]);

      // reset pr_curr for the next iteration
      pr_curr[v] = pr_next[v];
      pr_next[v] = 0.0;
    }
  }
  time_taken = t1.stop();
  // -------------------------------------------------------------------

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
          {"strategy", "Strategy to be used",
           cxxopts::value<uint>()->default_value(DEFAULT_STRATEGY)},
          {"inputFile", "Input graph file path",
           cxxopts::value<std::string>()->default_value(
               "/scratch/input_graphs/roadNet-CA")},
      });

  auto cl_options = options.parse(argc, argv);
  uint n_workers = cl_options["nWorkers"].as<uint>();
  uint strategy = cl_options["strategy"].as<uint>();
  uint max_iterations = cl_options["nIterations"].as<uint>();
  std::string input_file_path = cl_options["inputFile"].as<std::string>();

#ifdef USE_INT
  std::cout << "Using INT\n";
#else
  std::cout << "Using FLOAT\n";
#endif
  std::cout << std::fixed;
  std::cout << "Number of workers : " << n_workers << "\n";
  std::cout << "Task decomposition strategy : " << strategy << "\n";
  std::cout << "Iterations : " << max_iterations << "\n";

  Graph g;
  std::cout << "Reading graph\n";
  g.readGraphFromBinary<int>(input_file_path);
  std::cout << "Created graph\n";
  switch (strategy) {
  case 0:
    std::cout << "\nSerial\n";
    pageRankSerial(g, max_iterations);
    break;
  case 1:
    std::cout << "\nVertex-based work partitioning\n";
    pageRankStrategyOne(g, max_iterations, n_workers);
    break;
  case 2:
    std::cout << "\nEdge-based work partitioning\n";
    pageRankStrategyTwo(g, max_iterations, n_workers);
    break;
  default:
    break;
  }

  return 0;
}


void singleThread(uintV start, uintV end, PageRankType *pr_curr,  std::atomic<PageRankType> *pr_next, Graph &g, CustomBarrier *barrier, int thread_id, int max_iters)
{
  timer t1, t_barrier1_time, t_barrier2_time;
  double time_taken = 0.0;
  double barrier1_time = 0.0;
  double barrier2_time = 0.0;
  int num_edges = 0;
  int vertices_processed = 0;
  t1.start();
  for(int iter = 0; iter < max_iters; iter++)
  {
    timer t_barrier1_time, t_barrier2_time;
    double barrier1_time_current = 0.0;
    double barrier2_time_current = 0.0;
    for (uintV u = start; u < end; u++) {
      uintE out_degree = g.vertices_[u].getOutDegree();
      num_edges += out_degree;
      for (uintE i = 0; i < out_degree; i++) {
        uintV v = g.vertices_[u].getOutNeighbor(i);
        // Atomic Operation
        pr_next[v] += pr_curr[u] / out_degree;
      }
    }
    
    // barrier 1 time:
    t_barrier1_time.start();
    barrier->wait();
    barrier1_time_current = t_barrier1_time.stop();
    barrier1_time+=barrier1_time_current;

    for (uintV v = start; v <end ; v++) {
        vertices_processed ++;
        PageRankType temp = pr_next[v];
        pr_curr[v] = PAGE_RANK(temp);
        pr_next[v] = (PageRankType)0;
      }
    // barrier 2 time:
    t_barrier2_time.start();
    barrier->wait();
    barrier2_time_current = t_barrier2_time.stop();
    barrier2_time+=barrier2_time_current;
  }
  time_taken = t1.stop();
  //thread_id, num_vertices, num_edges, barrier1_time, barrier2_time, total_time
  printf("%d, %d, %d, %lf, %lf, %lf\n", thread_id, vertices_processed, num_edges, barrier1_time, barrier2_time, time_taken);
    
}


/*******************************************Strategy 1 *********************************************************************/
void pageRankStrategyOne(Graph &g, int max_iters, uint n_threads) {
  uintV n = g.n_;
  uintV node_per_thread = n/n_threads;

  PageRankType *pr_curr = new PageRankType[n];
  std::atomic<PageRankType> *pr_next = new std::atomic<PageRankType>[n];

  for (uintV i = 0; i < n; i++) {
    pr_curr[i] = INIT_PAGE_RANK;
    pr_next[i] = (PageRankType)0;
  }

  // Partition Time

  double partitioning_time = 0.0;
  timer t0;
  t0.start();
  int *start_vertex = new int[n_threads];
  int *end_vertex = new int[n_threads];

  for(int i=0; i<n_threads; i++)
  {
    start_vertex[i] = i*node_per_thread;
    if(i==n_threads-1)
      end_vertex[i] = n;
    else
      end_vertex[i] = node_per_thread*(i+1);
    
  }
  partitioning_time = t0.stop();

  // Push based pagerank
  timer t1;
  double time_taken = 0.0;
  t1.start();
  // Create threads and distribute the work across T threads
  // -------------------------------------------------------------------
  std::vector<std::thread> thread_vector;
  CustomBarrier my_barrier(n_threads);
  std::cout << "thread_id, num_vertices, num_edges, barrier1_time, barrier2_time, total_time\n";
  // for each vertex 'u', process all its outNeighbors 'v'
  for(int i=0; i<n_threads; i++)
  {
    thread_vector.push_back(std::thread(singleThread, start_vertex[i], end_vertex[i], pr_curr, pr_next, std::ref(g), &my_barrier, i, max_iters));
  }

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
  std::cout << "Partitioning time (in seconds) : " << std::setprecision(8)<< partitioning_time<< "\n";
  std::cout << "Time taken (in seconds) : " << time_taken << "\n";
  delete[] pr_curr;
  delete[] pr_next;
}

/*******************************************Strategy 2 *********************************************************************/
void pageRankStrategyTwo(Graph &g, int max_iters, uint n_threads) 
{
  uintV n = g.n_;
  uintE m = g.m_;
  uintV node_per_thread = n/n_threads;

  PageRankType *pr_curr = new PageRankType[n];
  std::atomic<PageRankType> *pr_next = new std::atomic<PageRankType>[n];

  for (uintV i = 0; i < n; i++) {
    pr_curr[i] = INIT_PAGE_RANK;
    pr_next[i] = (PageRankType)0;
  }

  // Partition Time

  double partitioning_time = 0.0;
  timer t0;
  t0.start();
  uintV single_edges = m / n_threads;
  int *start_vertex = new int[n_threads];
  int *end_vertex = new int[n_threads];
  
  uintV current_vertex = -1;
  uintV current_edges = 0;
  for(int i=0; i<n_threads; i++)
  {
    if(i != n_threads-1)
    {
      while(current_edges <= single_edges)
      {
        current_vertex++;
        current_edges += g.vertices_[current_vertex].getOutDegree();
      }
      current_edges = 0;
      end_vertex[i] = current_vertex;
      current_vertex--;
    }else
    {
      end_vertex[i] = n;
    }
  }
  start_vertex[0] = 0;
  for(int i=1; i<n_threads; i++)
  {
    start_vertex[i] = end_vertex[i-1];
  }
  partitioning_time = t0.stop();

  // Push based pagerank
  timer t1;
  double time_taken = 0.0;
  t1.start();
  // Create threads and distribute the work across T threads
  // -------------------------------------------------------------------
  std::vector<std::thread> thread_vector;
  CustomBarrier my_barrier(n_threads);
  std::cout << "thread_id, num_vertices, num_edges, barrier1_time, barrier2_time, total_time\n";
  // for each vertex 'u', process all its outNeighbors 'v'
  for(int i=0; i<n_threads; i++)
  {
    thread_vector.push_back(std::thread(singleThread, start_vertex[i], end_vertex[i], pr_curr, pr_next, std::ref(g), &my_barrier, i, max_iters));
  }
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
  std::cout << "Partitioning time (in seconds) : " << std::setprecision(8)<< partitioning_time<< "\n";
  std::cout << "Time taken (in seconds) : " << time_taken << "\n";
  delete[] pr_curr;
  delete[] pr_next;
}