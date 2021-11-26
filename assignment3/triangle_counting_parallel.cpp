#include "core/graph.h"
#include "core/utils.h"
#include <iomanip>
#include <iostream>
#include <stdlib.h>
#include <thread>

std::atomic<long> total_triangles(0);

void singleThreadTriangleCountStrategy(Graph &g, uintV startIndex, uintV endIndex, int thread_id);
void triangleCountStrategyOne(Graph &g, uint n_threads);
void triangleCountStrategyTwo(Graph &g, uint n_threads);
void triangleCountStrategyTwo(Graph &g, uint n_threads);


uintV countTriangles(uintV *array1, uintE len1, uintV *array2, uintE len2,
                     uintV u, uintV v) {
  uintE i = 0, j = 0; // indexes for array1 and array2
  uintV count = 0;

  if (u == v)
    return count;

  while ((i < len1) && (j < len2)) {
    if (array1[i] == array2[j]) {
      if ((array1[i] != u) && (array1[i] != v)) {
        count++;
      } else {
        // triangle with self-referential edge -> ignore
      }
      i++;
      j++;
    } else if (array1[i] < array2[j]) {
      i++;
    } else {
      j++;
    }
  }
  return count;
}

void triangleCountSerial(Graph &g) {
  uintV n = g.n_;
  long triangle_count = 0;
  double time_taken = 0.0;
  timer t1;
  t1.start();
  for (uintV u = 0; u < n; u++) {
    uintE out_degree = g.vertices_[u].getOutDegree();
    for (uintE i = 0; i < out_degree; i++) {
      uintV v = g.vertices_[u].getOutNeighbor(i);
      triangle_count += countTriangles(g.vertices_[u].getInNeighbors(),
                                       g.vertices_[u].getInDegree(),
                                       g.vertices_[v].getOutNeighbors(),
                                       g.vertices_[v].getOutDegree(), u, v);
    }
  }
  time_taken = t1.stop();
  std::cout << "Number of triangles : " << triangle_count << "\n";
  std::cout << "Number of unique triangles : " << triangle_count / 3 << "\n";
  std::cout << "Time taken (in seconds) : " << std::setprecision(TIME_PRECISION)
            << time_taken << "\n";
}





int main(int argc, char *argv[]) {
  cxxopts::Options options(
      "triangle_counting_serial",
      "Count the number of triangles using serial and parallel execution");
  options.add_options(
      "custom",
      {
          {"nWorkers", "Number of workers",
           cxxopts::value<uint>()->default_value(DEFAULT_NUMBER_OF_WORKERS)},
          {"strategy", "Strategy to be used",
           cxxopts::value<uint>()->default_value(DEFAULT_STRATEGY)},
          {"inputFile", "Input graph file path",
           cxxopts::value<std::string>()->default_value(
               "/scratch/input_graphs/roadNet-CA")},
      });

  auto cl_options = options.parse(argc, argv);
  uint n_workers = cl_options["nWorkers"].as<uint>();
  uint strategy = cl_options["strategy"].as<uint>();
  std::string input_file_path = cl_options["inputFile"].as<std::string>();
  std::cout << std::fixed;
  std::cout << "Number of workers : " << n_workers << "\n";
  std::cout << "Task decomposition strategy : " << strategy << "\n";

  Graph g;
  std::cout << "Reading graph\n";
  g.readGraphFromBinary<int>(input_file_path);
  std::cout << "Created graph\n";

  switch (strategy) {
  case 0:
    std::cout << "\nSerial\n";
    triangleCountSerial(g);
    break;
  case 1:
    std::cout << "\nVertex-based work partitioning\n";
    triangleCountStrategyOne(g, n_workers);
    break;
  case 2:
    std::cout << "\nEdge-based work partitioning\n";
    triangleCountStrategyTwo(g, n_workers);
    break;
  default:
    break;
  }

  return 0;
}



void singleThreadTriangleCountStrategy(Graph &g, uintV startIndex, uintV endIndex, int thread_id)
{

  double time_taken = 0.0;
  timer t1;
  t1.start();
  long triangle_count = 0;
  uintE num_edges = 0;
  for (uintV u = startIndex; u < endIndex; u++) {
    // For each outNeighbor v, find the intersection of inNeighbor(u) and
    // outNeighbor(v)
    uintE out_degree = g.vertices_[u].getOutDegree();
    num_edges += out_degree;
    for (uintE i = 0; i < out_degree; i++) {
      uintV v = g.vertices_[u].getOutNeighbor(i);
      triangle_count += countTriangles(g.vertices_[u].getInNeighbors(),
                                       g.vertices_[u].getInDegree(),
                                       g.vertices_[v].getOutNeighbors(),
                                       g.vertices_[v].getOutDegree(), u, v);
      
    }
  }
  // Aotmic add
  total_triangles += triangle_count;
  time_taken = t1.stop();
  printf("%d, %d, %d, %ld, %lf\n", thread_id, endIndex-startIndex ,num_edges,  triangle_count, time_taken);
}

/*******************************************Strategy 1 *********************************************************************/
void triangleCountStrategyOne(Graph &g, uint n_threads) {
  uintV n = g.n_;

  double partitioning_time = 0.0;
  timer t0;
  t0.start();
  uintV single_vertex = n / n_threads;
  int *start_vertex = new int[n_threads];
  int *end_vertex = new int[n_threads];

  for(int i=0; i<n_threads; i++)
  {
    start_vertex[i] = i*single_vertex;
    if(i==n_threads-1)
      end_vertex[i] = n;
    else
      end_vertex[i] = (i+1)*single_vertex;
  }
  partitioning_time = t0.stop();

  double time_taken = 0.0;
  std::vector<std::thread> threads_vec;
  timer t1;
  t1.start();
  printf( "thread_id, num_vertices, num_edges, triangle_count, time_taken\n");
  for(int i=0; i<n_threads; i++)
  {
    threads_vec.push_back(std::thread(singleThreadTriangleCountStrategy, std::ref(g), start_vertex[i], end_vertex[i], i));
  }
  for(int i=0; i<n_threads; i++)
  {
    threads_vec[i].join();
  }
    

  time_taken = t1.stop();
  // -------------------------------------------------------------------
  // Print the overall statistics
  std::cout << "Number of triangles : " << total_triangles << "\n";
  std::cout << "Number of unique triangles : " << total_triangles / 3 << "\n";
  //partitioning_time
  std::cout << "Partitioning time (in seconds) : " << std::setprecision(8)<< partitioning_time<< "\n";
  std::cout << "Time taken (in seconds) : " << std::setprecision(TIME_PRECISION)
            << time_taken << "\n";
}

/*******************************************Strategy 2 *********************************************************************/
void triangleCountStrategyTwo(Graph &g, uint n_threads) 
{
  uintE n = g.n_;
  uintE m = g.m_;
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

    double time_taken = 0.0;
    std::vector<std::thread> threads_vec;
    timer t1;
    t1.start();
    printf( "thread_id, num_vertices, num_edges, triangle_count, time_taken\n");
    for(int i=0; i<n_threads; i++)
    {
      threads_vec.push_back(std::thread(singleThreadTriangleCountStrategy, std::ref(g), start_vertex[i], end_vertex[i], i));
    }
    for(int i=0; i<n_threads; i++)
    {
      threads_vec[i].join();
    }
      

    time_taken = t1.stop();
    // -------------------------------------------------------------------
    // Print the overall statistics
    std::cout << "Number of triangles : " << total_triangles << "\n";
    std::cout << "Number of unique triangles : " << total_triangles / 3 << "\n";
    //partitioning_time
    std::cout << "Partitioning time (in seconds) : " << std::setprecision(8)<< partitioning_time<< "\n";
    std::cout << "Time taken (in seconds) : " << std::setprecision(TIME_PRECISION)
              << time_taken << "\n";
}