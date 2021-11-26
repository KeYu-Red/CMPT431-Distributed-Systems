#include <iostream>
#include <cstdio>
#include "core/utils.h"
#include "core/graph.h"

#include <mpi.h>

const int ROOT_PROCESS = 0;

void triangleCountParallelStrategyOne(Graph &g, int world_rank, int world_size);
void triangleCountParallelStrategyTwo(Graph &g, int world_rank, int world_size);
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

void triangleCountSerial(Graph &g)
{
    uintV n = g.n_;
    long triangle_count = 0;
    double time_taken;
    timer t1;
    t1.start();
    for (uintV u = 0; u < n; u++)
    {
        uintE out_degree = g.vertices_[u].getOutDegree();
        for (uintE i = 0; i < out_degree; i++)
        {
            uintV v = g.vertices_[u].getOutNeighbor(i);
            triangle_count += countTriangles(g.vertices_[u].getInNeighbors(),
                                             g.vertices_[u].getInDegree(),
                                             g.vertices_[v].getOutNeighbors(),
                                             g.vertices_[v].getOutDegree(),
                                             u,
                                             v);
        }
    }

    // For every thread, print out the following statistics:
    // rank, edges, triangle_count, communication_time
    // 0, 17248443, 144441858, 0.000074
    // 1, 17248443, 152103585, 0.000020
    // 2, 17248443, 225182666, 0.000034
    // 3, 17248444, 185596640, 0.000022

    time_taken = t1.stop();

    // Print out overall statistics
    std::printf("Number of triangles : %ld\n", triangle_count);
    std::printf("Number of unique triangles : %ld\n", triangle_count / 3);
    std::printf("Time taken (in seconds) : %f\n", time_taken);
}

int main(int argc, char *argv[])
{
    cxxopts::Options options("triangle_counting_serial", "Count the number of triangles using serial and parallel execution");
    options.add_options("custom", {
                                      {"strategy", "Strategy to be used", cxxopts::value<uint>()->default_value(DEFAULT_STRATEGY)},
                                      {"inputFile", "Input graph file path", cxxopts::value<std::string>()->default_value("/scratch/input_graphs/roadNet-CA")},
                                  });

    auto cl_options = options.parse(argc, argv);
    uint strategy = cl_options["strategy"].as<uint>();
    std::string input_file_path = cl_options["inputFile"].as<std::string>();
    
    // Initialize the MPI environment
    MPI_Init(nullptr, nullptr);

    // Get the number of processes
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    // Get the rank of the process
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);


    // Get the world size and print it out here
    if(world_rank == ROOT_PROCESS)
    {
      std::printf("World size : %d\n", world_size);
      std::printf("Communication strategy : %d\n", strategy);
      std::printf("rank, edges, triangle_count, communication_time\n");
    }
    

    Graph g;
    g.readGraphFromBinary<int>(input_file_path);

    switch(strategy)
    {
      case 0:
        if(world_rank == ROOT_PROCESS)
          triangleCountSerial(g);
        break;
      case 1:
        triangleCountParallelStrategyOne(g, world_rank, world_size);
        break;
      case 2:
        triangleCountParallelStrategyTwo(g, world_rank, world_size);
        break;
    }

    MPI_Finalize();
    
    return 0;
}



void triangleCountParallelStrategyOne(Graph &g, int world_rank, int world_size)
{
    uintV n = g.n_;
    uintE m = g.m_;
    uintV single_edges = m/world_size;
    long total_triangle_count = 0;
    long local_count = 0;

    // Define the time of ROOT_PROCESS is the total time
    double total_time_taken;
    timer total_timer;
    if(world_rank == ROOT_PROCESS)
    {
      total_timer.start();
    }

    // Define the start, end vertex of the process
    // Using edge decomposition strategy 
    uintV startIndex =0 , endIndex = 0;
    for(int i = ROOT_PROCESS; i < world_size; i++) {
        startIndex = endIndex;
        if(i == world_size-1)
        {
          endIndex = n;
          break;
        }
        long process_edges = 0;
        while(endIndex < g.n_) {
            process_edges += g.vertices_[endIndex].getOutDegree();
            endIndex ++;
            if(process_edges >= single_edges) {
                break;
            }
        }
        if(i == world_rank) {
            break;
        }
    }

    // Calculate the triangles
    uintE local_edges = 0;
    for (uintV u = startIndex; u < endIndex; u++)
    {
        uintE out_degree = g.vertices_[u].getOutDegree();
        local_edges += out_degree;
        for (uintE i = 0; i < out_degree; i++)
        {
            uintV v = g.vertices_[u].getOutNeighbor(i);
            local_count += countTriangles(g.vertices_[u].getInNeighbors(),
                                             g.vertices_[u].getInDegree(),
                                             g.vertices_[v].getOutNeighbors(),
                                             g.vertices_[v].getOutDegree(),
                                             u,
                                             v);
        }
    }

    //std::printf("Finish Counting, current count=%ld, current process=%d\n ", local_count, world_rank);
    double communicate_time_taken;
    timer communicate_timer; 
    communicate_timer.start();
    if(world_rank != ROOT_PROCESS)
    {
      // int MPI_Gather(const void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf, int recvcount, MPI_Datatype recvtype, int root, MPI_Comm comm)
      MPI_Gather(&local_count, 1, MPI_LONG, nullptr, 0, MPI_LONG, ROOT_PROCESS, MPI_COMM_WORLD);
    }
    else{
      long *buffer = new long[world_size];
      MPI_Gather(&local_count, 1, MPI_LONG, buffer, 1, MPI_LONG, ROOT_PROCESS, MPI_COMM_WORLD);
      for (int i = 0; i < world_size; i++)
      {
        total_triangle_count+= buffer[i];
      }
      delete[] buffer;
 
    }
    communicate_time_taken = communicate_timer.stop();

    std::printf("%d, %u, %ld ,%lf\n",world_rank, local_edges, local_count, communicate_time_taken );

    if(world_rank == ROOT_PROCESS)
    {
      total_time_taken = total_timer.stop();
    // Print out overall statistics
    std::printf("Number of triangles : %ld\n", total_triangle_count);
    std::printf("Number of unique triangles : %ld\n", total_triangle_count / 3);
    std::printf("Time taken (in seconds) : %f\n", total_time_taken);
    }
}

void triangleCountParallelStrategyTwo(Graph &g, int world_rank, int world_size)
{
    uintV n = g.n_;
    uintE m = g.m_;
    uintV single_edges = m/world_size;
    long total_triangle_count = 0;
    long local_count = 0;

    // Define the time of ROOT_PROCESS is the total time
    double total_time_taken;
    timer total_timer;
    if(world_rank == ROOT_PROCESS)
    {
      total_timer.start();
    }

    // Define the start, end vertex of the process
    // Using edge decomposition strategy 
    uintV startIndex =0 , endIndex = 0;
    for(int i = ROOT_PROCESS; i < world_size; i++) {
        startIndex = endIndex;
        if(i == world_size-1)
        {
          endIndex = n;
          break;
        }
        long process_edges = 0;
        while(endIndex < g.n_) {
            process_edges += g.vertices_[endIndex].getOutDegree();
            endIndex ++;
            if(process_edges >= single_edges) {
                break;
            }
        }
        if(i == world_rank) {
            break;
        }
    }

    // Calculate the triangles
    uintE local_edges = 0;
    for (uintV u = startIndex; u < endIndex; u++)
    {
        uintE out_degree = g.vertices_[u].getOutDegree();
        local_edges += out_degree;
        for (uintE i = 0; i < out_degree; i++)
        {
            uintV v = g.vertices_[u].getOutNeighbor(i);
            local_count += countTriangles(g.vertices_[u].getInNeighbors(),
                                             g.vertices_[u].getInDegree(),
                                             g.vertices_[v].getOutNeighbors(),
                                             g.vertices_[v].getOutDegree(),
                                             u,
                                             v);
        }
    }

    //std::printf("Finish Counting, current count=%ld, current process=%d\n ", local_count, world_rank);
    double communicate_time_taken;
    timer communicate_timer; 
    communicate_timer.start();
    if(world_rank != ROOT_PROCESS)
    {
      // int MPI_Reduce(const void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype, MPI_Op op, int root, MPI_Comm comm)
      MPI_Reduce(&local_count, nullptr, 1, MPI_LONG, MPI_SUM, ROOT_PROCESS, MPI_COMM_WORLD);
    }
    else{
      MPI_Reduce(&local_count, &total_triangle_count, 1, MPI_LONG, MPI_SUM, ROOT_PROCESS, MPI_COMM_WORLD);
    }
    communicate_time_taken = communicate_timer.stop();

    std::printf("%d, %u, %ld ,%lf\n",world_rank, local_edges, local_count, communicate_time_taken );

    if(world_rank == ROOT_PROCESS)
    {
      total_time_taken = total_timer.stop();
    // Print out overall statistics
    std::printf("Number of triangles : %ld\n", total_triangle_count);
    std::printf("Number of unique triangles : %ld\n", total_triangle_count / 3);
    std::printf("Time taken (in seconds) : %f\n", total_time_taken);
    }
}