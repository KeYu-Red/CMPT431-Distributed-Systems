#include "core/graph.h"
#include "core/utils.h"
#include <iomanip>
#include <iostream>
#include <stdlib.h>
#include <thread>

std::atomic<long> total_triangles(0);
std::atomic<uintV> current_vertex(0);
const uintV INVALID_VERTEX = -1;

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

uintV getNextVertexToBeProcessed(const int max)
{
  if(current_vertex >= max) return INVALID_VERTEX;
  uintV current = (uintV)current_vertex++;
  return current;
}

void singleThreadTriangleCount(Graph& g, const int max, const int thread_id)
{

  double time_taken = 0.0;
  timer t1;
  t1.start();
  long triangle_count = 0;
  uintE num_edges = 0;
  uintV num_vertex = 0;
  while(true)
  {
    uintV u = getNextVertexToBeProcessed(max);
    if(u == INVALID_VERTEX) break;
    num_vertex++;
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
  printf("%d, %d, %d, %ld, %lf\n", thread_id, num_vertex ,num_edges,  triangle_count, time_taken);
}

void triangleCountParallel(Graph &g, int n_threads) {
  uintV n = g.n_;
  double time_taken = 0.0;
  printf( "thread_id, num_vertices, num_edges, triangle_count, time_taken\n");
  std::vector<std::thread> threads_vec;
  timer t1;
  t1.start();
  for(int i=0; i<n_threads; i++)
  {
    threads_vec.push_back(std::thread(singleThreadTriangleCount, std::ref(g), n, i));
  }
  for(int i=0; i<n_threads; i++)
  {
    threads_vec[i].join();
  }

  time_taken = t1.stop();
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
               "/scratch/input_graphs/roadNet-CA")},
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

  triangleCountParallel(g, n_workers);
  return 0;
}
