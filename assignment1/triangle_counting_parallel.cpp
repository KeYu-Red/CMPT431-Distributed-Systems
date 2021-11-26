#include "core/graph.h"
#include "core/utils.h"
#include <future>
#include <iomanip>
#include <iostream>
#include <stdlib.h>
#include <thread>

std::atomic<long> total_triangles(0);
std::atomic<int> current_thread_id(0);


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

void singleThreadTriangleCount(Graph &g, uintV startIndex, uintV endIndex)
{

  double time_taken = 0.0;
  timer t1;
  t1.start();
  long triangle_count = 0;
  for (uintV u = startIndex; u < endIndex; u++) {
    // For each outNeighbor v, find the intersection of inNeighbor(u) and
    // outNeighbor(v)
    uintE out_degree = g.vertices_[u].getOutDegree();
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
//  std::cout << "thread_id, triangle_count, time_taken\n";
  //std::cout << current_thread_id++<<", "<< triangle_count<<", "<< time_taken << "\n";
  printf("%d, %ld, %lf\n", int(current_thread_id++), triangle_count, time_taken);
}

void triangleCountSerial(Graph &g, uint n_threads) {
  uintV n = g.n_;
  uintV single_vertex = n / n_threads;
  double time_taken = 0.0;
  std::vector<std::thread> threads_vec;
  timer t1;
  t1.start();
  printf( "thread_id, triangle_count, time_taken\n");
  for(int i=0; i<n_threads; i++)
  {
    if(i == n_threads-1)
      threads_vec.push_back(std::thread(singleThreadTriangleCount, std::ref(g), i*single_vertex, n));
    else
      threads_vec.push_back(std::thread(singleThreadTriangleCount, std::ref(g), i*single_vertex, (i+1)*single_vertex));
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
          {"inputFile", "Input graph file path",
           cxxopts::value<std::string>()->default_value(
               "/scratch/assignment1/input_graphs/roadNet-CA")},
      });

  auto cl_options = options.parse(argc, argv);
  uint n_workers = cl_options["nWorkers"].as<uint>();
  std::string input_file_path = cl_options["inputFile"].as<std::string>();
  std::cout << std::fixed;
  std::cout << "Number of workers : " << n_workers << "\n";

  Graph g;
  std::cout << "Reading graph\n";
  g.readGraphFromBinary<int>(input_file_path);
  std::cout << "Created graph\n";

  triangleCountSerial(g, n_workers);

  return 0;
}