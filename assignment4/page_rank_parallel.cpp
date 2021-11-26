#include "core/graph.h"
#include "core/utils.h"
#include <iomanip>
#include <iostream>
#include <stdlib.h>
#include <thread>
#include <atomic>

const uintV INVALID_VERTEX = -1;
std::atomic<uintV> current_vertex(0);

//uintV current_vertex = 0;
std::mutex getNextVertexLock;

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
uintV getNextVertexToBeProcessed(const int max, const int granularity);
void singleThreadFuntion(PageRankType *pr_curr, std::atomic<PageRankType> *pr_next, Graph &g, CustomBarrier *barrier, int const thread_id, int const max_iters, const uintV max, const int n_threads, const int granularity);
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
//pageRankSerial(g, max_iterations, n_workers, granularity);
void pageRankParallel(Graph &g, uint max_iters, int n_threads, uint granularity) {
  uintV n = g.n_;
  PageRankType *pr_curr = new PageRankType[n];
  std::atomic<PageRankType> *pr_next = new std::atomic<PageRankType>[n];

  for (uintV i = 0; i < n; i++) {
    pr_curr[i] = INIT_PAGE_RANK;
    pr_next[i] = (PageRankType)0;
  }

  // Push based pagerank
  timer t1;
  double time_taken = 0.0;
  t1.start();
  // Create threads and distribute the work across T threads
  // -------------------------------------------------------------------
  std::vector<std::thread> thread_vector;
  CustomBarrier my_barrier(n_threads);
  std::cout << "thread_id, num_vertices, num_edges, barrier1_time, barrier2_time, getNextVertex_time, total_time\n";
  for(int i=0; i<n_threads; i++)
  {
    thread_vector.push_back(std::thread(singleThreadFuntion, pr_curr, pr_next, std::ref(g), &my_barrier, i, max_iters, n, n_threads, granularity));
  }

  for(int i=0; i<n_threads; i++)
  {
      thread_vector[i].join();
  }
  // -------------------------------------------------------------------

  PageRankType sum_of_page_ranks = 0;
  for (uintV u = 0; u < n; u++) {
    sum_of_page_ranks += pr_curr[u];
  }
  time_taken = t1.stop();
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
          {"granularity", "Granularity to be used", cxxopts::value<uint>()->default_value(DEFAULT_GRANULARITY)},
          {"inputFile", "Input graph file path",
           cxxopts::value<std::string>()->default_value(
               "/scratch/input_graphs/roadNet-CA")},
      });

  auto cl_options = options.parse(argc, argv);
  uint n_workers = cl_options["nWorkers"].as<uint>();
  uint max_iterations = cl_options["nIterations"].as<uint>();
  uint granularity = cl_options["granularity"].as<uint>();    
  std::string input_file_path = cl_options["inputFile"].as<std::string>();

#ifdef USE_INT
  std::cout << "Using INT\n";
#else
  std::cout << "Using FLOAT\n";
#endif
  std::cout << std::fixed;
  std::cout << "Number of workers : " << n_workers << "\n";
  std::cout << "Granularity : " << granularity << "\n";
  std::cout << "Iterations : " << max_iterations << "\n";

  Graph g;
  std::cout << "Reading graph\n";
  g.readGraphFromBinary<int>(input_file_path);
  std::cout << "Created graph\n";

  pageRankParallel(g, max_iterations, n_workers, granularity);
  return 0;
}

uintV getNextVertexToBeProcessed(const int max, const int granularity)
{
  // getNextVertexLock.lock();
  // uintV current = INVALID_VERTEX;
  // if(current_vertex < max)
  // {
  //   current = (uintV)current_vertex;
  //   current_vertex += granularity;
  // }
  // getNextVertexLock.unlock();
  // return current;
  
  uintV current = current_vertex.fetch_add(granularity);
  return (current < max) ?current: INVALID_VERTEX;

}

void singleThreadFuntion(PageRankType *pr_curr, std::atomic<PageRankType> *pr_next, Graph &g, CustomBarrier *barrier, int const thread_id, int const max_iters, const uintV max, const int n_threads, const int granularity)
{
  timer t1;
  double time_taken = 0.0;
  double barrier1_time = 0.0;
  double barrier2_time = 0.0;
  double getNextVertex_time = 0.0;
  int num_edges = 0;
  int vertices_processed = 0;
  t1.start();
  for(int iter = 0; iter < max_iters; iter++)
  {
  
    timer t_barrier1_time, t_barrier2_time;
    double barrier1_time_current = 0.0;
    double barrier2_time_current = 0.0;
    std::vector<uintV> processed_vectex;
    while(true) {
      timer t_getNextVertex_time;
      double getNextVertex_time_current = 0.0;
      t_getNextVertex_time.start();
      uintV u = getNextVertexToBeProcessed(max, granularity);
      getNextVertex_time_current = t_getNextVertex_time.stop();
      getNextVertex_time+=getNextVertex_time_current;
      if(u == INVALID_VERTEX) break;
      for(int j=0; j<granularity; j++)
      {
        uintE out_degree = g.vertices_[u].getOutDegree();
        num_edges += out_degree;
        for (uintE i = 0; i < out_degree; i++) {
          uintV v = g.vertices_[u].getOutNeighbor(i);
          // Atomic Operation
          pr_next[v] += pr_curr[u] / out_degree; 
        }
        u++;
        if(u>=max) break;
      }
      
    }
    // barrier 1 time:
    t_barrier1_time.start();
    barrier->wait();
    barrier1_time_current = t_barrier1_time.stop();
    barrier1_time+=barrier1_time_current;

    if(thread_id == 0)
      current_vertex = 0;
    barrier->wait();

    while(true) {
        timer t_getNextVertex_time;
        double getNextVertex_time_current = 0.0;
        t_getNextVertex_time.start();
        uintV v = getNextVertexToBeProcessed(max, granularity);
        getNextVertex_time_current = t_getNextVertex_time.stop();
        getNextVertex_time+=getNextVertex_time_current;
        if(v == INVALID_VERTEX) break;
        for(int j=0; j<granularity; j++)
        {
          vertices_processed++;
          PageRankType temp = pr_next[v];
          pr_curr[v] = PAGE_RANK(temp);
          pr_next[v] = (PageRankType)0;
          v++;
          if(v>=max) break;
        }
      }

    // barrier 2 time:
    t_barrier2_time.start();
    barrier->wait();
    barrier2_time_current = t_barrier2_time.stop();
    barrier2_time+=barrier2_time_current;

    if(thread_id == 0)
      current_vertex = 0;
    barrier->wait();
  }
  time_taken = t1.stop();
  //thread_id, num_vertices, num_edges, barrier1_time, barrier2_time, getNextVertex_time, total_time
  printf("%d, %d, %d, %lf, %lf, %lf, %lf\n", thread_id, vertices_processed, num_edges, barrier1_time, barrier2_time, getNextVertex_time, time_taken);
    
}
