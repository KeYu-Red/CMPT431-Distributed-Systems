#include <iostream>
#include <cstdio>
#include "core/utils.h"
#include "core/graph.h"
#include <mpi.h>

#ifdef USE_INT
#define INIT_PAGE_RANK 100000
#define EPSILON 1000
#define PAGE_RANK(x) (15000 + (5 * x) / 6)
#define CHANGE_IN_PAGE_RANK(x, y) std::abs(x - y)
#define PR_FMT "%ld"
#define ROOT_PROCESS 0
#define PAGERANK_MPI_TYPE MPI_LONG
typedef int64_t PageRankType;
#else
#define INIT_PAGE_RANK 1.0
#define EPSILON 0.01
#define DAMPING 0.85
#define PAGE_RANK(x) (1 - DAMPING + DAMPING * x)
#define CHANGE_IN_PAGE_RANK(x, y) std::fabs(x - y)
#define PR_FMT "%f"
#define PAGERANK_MPI_TYPE MPI_FLOAT
#define ROOT_PROCESS 0
typedef float PageRankType;
#endif

void pageRankParallelStrategyOne(Graph &g, int max_iters, int world_size, int world_rank);
void pageRankParallelStrategyTwo(Graph &g, int max_iters, int world_size, int world_rank);

void pageRankSerial(Graph &g, int max_iters)
{
    uintV n = g.n_;
    double time_taken;
    timer t1;
    PageRankType *pr_curr = new PageRankType[n];
    PageRankType *pr_next = new PageRankType[n];

    t1.start();
    for (uintV i = 0; i < n; i++)
    {
        pr_curr[i] = INIT_PAGE_RANK;
        pr_next[i] = 0.0;
    }

    // Push based pagerank
    // -------------------------------------------------------------------
    for (int iter = 0; iter < max_iters; iter++)
    {
        // for each vertex 'u', process all its outNeighbors 'v'
        for (uintV u = 0; u < n; u++)
        {
            uintE out_degree = g.vertices_[u].getOutDegree();
            for (uintE i = 0; i < out_degree; i++)
            {
                uintV v = g.vertices_[u].getOutNeighbor(i);
                pr_next[v] += (pr_curr[u] / out_degree);
            }
        }
        for (uintV v = 0; v < n; v++)
        {
            pr_next[v] = PAGE_RANK(pr_next[v]);

            // reset pr_curr for the next iteration
            pr_curr[v] = pr_next[v];
            pr_next[v] = 0.0;
        }
    }
    // -------------------------------------------------------------------

    // For every thread, print the following statistics:
    // rank, num_edges, communication_time
    // 0, 344968860, 1.297778
    // 1, 344968860, 1.247763
    // 2, 344968860, 0.956243
    // 3, 344968880, 0.467028

    PageRankType sum_of_page_ranks = 0;
    for (uintV u = 0; u < n; u++)
    {
        sum_of_page_ranks += pr_curr[u];
    }
    time_taken = t1.stop();
    std::printf("Sum of page rank : " PR_FMT "\n", sum_of_page_ranks);
    std::printf("Time taken (in seconds) : %f\n", time_taken);
    delete[] pr_curr;
    delete[] pr_next;
}

int main(int argc, char *argv[])
{
    cxxopts::Options options("page_rank_push", "Calculate page_rank using serial and parallel execution");
    options.add_options("", {
                                {"nIterations", "Maximum number of iterations", cxxopts::value<uint>()->default_value(DEFAULT_MAX_ITER)},
                                {"strategy", "Strategy to be used", cxxopts::value<uint>()->default_value(DEFAULT_STRATEGY)},
                                {"inputFile", "Input graph file path", cxxopts::value<std::string>()->default_value("/scratch/input_graphs/roadNet-CA")},
                            });

    auto cl_options = options.parse(argc, argv);
    uint strategy = cl_options["strategy"].as<uint>();
    uint max_iterations = cl_options["nIterations"].as<uint>();
    std::string input_file_path = cl_options["inputFile"].as<std::string>();

    // Initialize the MPI environment
    MPI_Init(nullptr, nullptr);

    // Get the number of processes
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    // Get the rank of the process
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    if(world_rank == ROOT_PROCESS)
    {

#ifdef USE_INT
    std::printf("Using INT\n");
#else
    std::printf("Using FLOAT\n");
#endif
    // Get the world size and print it out here
    std::printf("World size : %d\n", world_size);
    std::printf("Communication strategy : %d\n", strategy);
    std::printf("Iterations : %d\n", max_iterations);
    std::printf("rank, num_edges, communication_time\n");
    }

    Graph g;
    g.readGraphFromBinary<int>(input_file_path);

    switch(strategy)
    {
        case 0:
            if(world_rank == ROOT_PROCESS)
                pageRankSerial(g, max_iterations);
            break;
        case 1:
            pageRankParallelStrategyOne(g, max_iterations, world_size, world_rank);
            break;
        case 2:
            pageRankParallelStrategyTwo(g, max_iterations, world_size, world_rank);
            break;
    }

    MPI_Finalize();
    return 0;
}

void pageRankParallelStrategyOne(Graph &g, int max_iters, int world_size, int world_rank)
{
    uintV n = g.n_;
    uintE m = g.m_;
    // Define the time of ROOT_PROCESS is the total time
    double total_time_taken;
    timer total_timer;
    if(world_rank == ROOT_PROCESS)
    {
      total_timer.start();
    }
    // Define the communication timer
    double communicate_time_taken;
    timer communicate_timer;

    // Define the start/end of the vertex for each process
    // Using edge decomposition strategy 
    uintV startIndex = 0 , endIndex = 0;
    uintE single_edges = m / world_size;
    int *startIndexPtr = new int[world_size];
    int *endIndexPtr = new int[world_size];
    int *sendCountsArray = new int[world_size];
    
    uintV current_vertex = -1;
    uintV current_edges = 0;
    for(int i=0; i<world_size; i++)
    {
        if(i != world_size-1)
        {
        while(current_edges <= single_edges)
        {
            current_vertex++;
            current_edges += g.vertices_[current_vertex].getOutDegree();
        }
        current_edges = 0;
        endIndexPtr[i] = current_vertex;
        current_vertex--;
        }else
        {
        endIndexPtr[i] = n;
        }
    }
    startIndexPtr[0] = 0;
    sendCountsArray[ROOT_PROCESS]=0;
    for(int i=1; i<world_size; i++)
    {
        startIndexPtr[i] = endIndexPtr[i-1];
        sendCountsArray[i]=endIndexPtr[i]-startIndexPtr[i];
    }
    startIndex = startIndexPtr[world_rank];
    endIndex = endIndexPtr[world_rank];

    PageRankType *pr_curr = new PageRankType[n];
    PageRankType *pr_next = new PageRankType[n];
    PageRankType *pr_temp = new PageRankType[n];

    for (uintV i = 0; i < n; i++)
    {
        pr_curr[i] = INIT_PAGE_RANK;
        pr_next[i] = (PageRankType)0;
        pr_temp[i] = (PageRankType)0;
    }
    uintE num_edges = 0;

    // Push based pagerank
    // -------------------------------------------------------------------
    for (int iter = 0; iter < max_iters; iter++)
    {
        // for each vertex 'u', process all its outNeighbors 'v'
        for (uintV u = startIndex ; u < endIndex; u++)
        {
            uintE out_degree = g.vertices_[u].getOutDegree();
            num_edges += out_degree;
            for (uintE i = 0; i < out_degree; i++)
            {
                uintV v = g.vertices_[u].getOutNeighbor(i);
                pr_next[v] += (pr_curr[u] / out_degree);
            }
        }

        // Communication Section:
        communicate_timer.start();
/*
sendbuf
[in] address of send buffer (choice, significant only at root)
sendcounts
[in] integer array (of length group size) specifying the number of elements to send to each processor
displs
[in] integer array (of length group size). Entry i specifies the displacement (relative to sendbuf from which to take the outgoing data to process i
sendtype
[in] data type of send buffer elements (handle)
recvbuf
[out] address of receive buffer (choice)
recvcount
[in] number of elements in receive buffer (integer)
recvtype
[in] data type of receive buffer elements (handle)
root
[in] rank of sending process (integer)
comm
[in] communicator (handle)
*/
        // int MPI_Reduce(const void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype, MPI_Op op, int root, MPI_Comm comm)
        if(world_rank != ROOT_PROCESS)
        {
            // int MPI_Reduce(const void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype, MPI_Op op, int root, MPI_Comm comm)
            MPI_Reduce(pr_next, nullptr, n, PAGERANK_MPI_TYPE, MPI_SUM, ROOT_PROCESS, MPI_COMM_WORLD);
            // int MPI_Scatterv(const void *sendbuf, const int sendcounts[], const int displs[], MPI_Datatype sendtype, void *recvbuf, int recvcount, MPI_Datatype recvtype, int root, MPI_Comm comm)
            MPI_Scatterv(nullptr, sendCountsArray, startIndexPtr, PAGERANK_MPI_TYPE, pr_temp, endIndex - startIndex, PAGERANK_MPI_TYPE, ROOT_PROCESS, MPI_COMM_WORLD);
            for (int i = startIndex; i < endIndex; i++)
             {
                 pr_next[i] = pr_temp[i - startIndex];
             }
        }else{
            // int MPI_Reduce(const void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype, MPI_Op op, int root, MPI_Comm comm)
            MPI_Reduce(pr_next, pr_temp, n, PAGERANK_MPI_TYPE, MPI_SUM, 0, MPI_COMM_WORLD);
            // int MPI_Scatterv(const void *sendbuf, const int sendcounts[], const int displs[], MPI_Datatype sendtype, void *recvbuf, int recvcount, MPI_Datatype recvtype, int root, MPI_Comm comm)
            MPI_Scatterv(pr_temp, sendCountsArray, startIndexPtr, PAGERANK_MPI_TYPE, nullptr, 0, PAGERANK_MPI_TYPE, ROOT_PROCESS, MPI_COMM_WORLD);
            for (int i = startIndex; i < endIndex; i++)
            {
                pr_next[i] = pr_temp[i - startIndex];
            }
            
        }
        
        communicate_time_taken += communicate_timer.stop();
        for (uintV v = 0; v < n; v++)
        {
            if(v >= startIndex && v < endIndex)
            {
                pr_next[v] = PAGE_RANK(pr_next[v]);
                pr_curr[v] = pr_next[v];
            }
            // reset pr_curr for the next iteration
            pr_next[v] = 0.0;
        }
    }
    std::printf("%d, %u, %lf\n",world_rank, num_edges, communicate_time_taken);

    PageRankType local_sum_of_page_ranks = 0, total_page_ranks=0;

    for (uintV u = startIndex; u < endIndex; u++)
    {
        local_sum_of_page_ranks += pr_curr[u];
    }
    //std::printf("process=%d, Sum of page rank : %ld\n", world_rank, local_count);

    // Get and sum all the processes
    if(world_rank != ROOT_PROCESS)
    {
        MPI_Reduce(&local_sum_of_page_ranks, nullptr, 1, PAGERANK_MPI_TYPE, MPI_SUM, 0, MPI_COMM_WORLD); 
    }else{
        MPI_Reduce(&local_sum_of_page_ranks, &total_page_ranks, 1, PAGERANK_MPI_TYPE, MPI_SUM, 0, MPI_COMM_WORLD);
    }

    if(world_rank == ROOT_PROCESS)
    {
        total_time_taken=total_timer.stop();
        std::printf("Sum of page rank : " PR_FMT "\n", total_page_ranks);
        std::printf("Time taken (in seconds) : %f\n", total_time_taken);
    }

    //time_taken = t1.stop();
    //std::printf("Time taken (in seconds) : %f\n", time_taken);
    delete[] pr_curr;
    delete[] pr_next;
    delete[] startIndexPtr;
    delete[] endIndexPtr;
    delete[] sendCountsArray;
    delete[] pr_temp;
}


void pageRankParallelStrategyTwo(Graph &g, int max_iters, int world_size, int world_rank)
{
    uintV n = g.n_;
    uintE m = g.m_;
    // Define the time of ROOT_PROCESS is the total time
    double total_time_taken;
    timer total_timer;
    if(world_rank == ROOT_PROCESS)
    {
      total_timer.start();
    }
    // Define the communication timer
    double communicate_time_taken;
    timer communicate_timer;

    // Define the start/end of the vertex for each process
    // Using edge decomposition strategy 
    uintV startIndex = 0 , endIndex = 0;
    uintE single_edges = m / world_size;
    int *startIndexPtr = new int[world_size];
    int *endIndexPtr = new int[world_size];
    
    uintV current_vertex = -1;
    uintV current_edges = 0;
    for(int i=0; i<world_size; i++)
    {
        if(i != world_size-1)
        {
        while(current_edges <= single_edges)
        {
            current_vertex++;
            current_edges += g.vertices_[current_vertex].getOutDegree();
        }
        current_edges = 0;
        endIndexPtr[i] = current_vertex;
        current_vertex--;
        }else
        {
        endIndexPtr[i] = n;
        }
    }
    startIndexPtr[0] = 0;
    for(int i=1; i<world_size; i++)
    {
        startIndexPtr[i] = endIndexPtr[i-1];
    }
    startIndex = startIndexPtr[world_rank];
    endIndex = endIndexPtr[world_rank];

    PageRankType *pr_curr = new PageRankType[n];
    PageRankType *pr_next = new PageRankType[n];
    PageRankType *pr_temp = new PageRankType[n];;

    for (uintV i = 0; i < n; i++)
    {
        pr_curr[i] = INIT_PAGE_RANK;
        pr_next[i] = (PageRankType)0;
        pr_temp[i] = (PageRankType)0;
    }
    uintE num_edges = 0;

    // Push based pagerank
    // -------------------------------------------------------------------
    for (int iter = 0; iter < max_iters; iter++)
    {
        // for each vertex 'u', process all its outNeighbors 'v'
        for (uintV u = startIndex ; u < endIndex; u++)
        {
            uintE out_degree = g.vertices_[u].getOutDegree();
            num_edges += out_degree;
            for (uintE i = 0; i < out_degree; i++)
            {
                uintV v = g.vertices_[u].getOutNeighbor(i);
                pr_next[v] += (pr_curr[u] / out_degree);
            }
        }

        // Communication Section:
        communicate_timer.start();


        for (int i = ROOT_PROCESS; i < world_size; i++)
        {
            uintV reduce_start = startIndexPtr[i];
            uintV reduce_end = endIndexPtr[i];
            if(world_rank==ROOT_PROCESS)
            {
                // int MPI_Reduce(const void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype, MPI_Op op, int root, MPI_Comm comm)
                MPI_Reduce(pr_next + reduce_start, pr_temp, reduce_end - reduce_start, PAGERANK_MPI_TYPE, MPI_SUM, i, MPI_COMM_WORLD);
                if (i == 0)
                {
                    for (int i = startIndex; i < endIndex; i++)
                    {
                        pr_next[i] = pr_temp[i - startIndex];
                    }
                }
            }else{
                if(i==world_rank)
                {
                    MPI_Reduce(pr_next + startIndex, pr_temp, endIndex-startIndex, PAGERANK_MPI_TYPE, MPI_SUM, world_rank, MPI_COMM_WORLD);
                    for (int i = startIndex; i < endIndex; i++)
                    {
                        pr_next[i] = pr_temp[i - startIndex];
                    }
                }else{
                    MPI_Reduce(pr_next + reduce_start, nullptr, reduce_end - reduce_start, PAGERANK_MPI_TYPE, MPI_SUM, i, MPI_COMM_WORLD); 
                }
                    
            }
        }
        communicate_time_taken += communicate_timer.stop();
        for (uintV v = 0; v < n; v++)
        {
            if(v >= startIndex && v < endIndex)
            {
                pr_next[v] = PAGE_RANK(pr_next[v]);
                pr_curr[v] = pr_next[v];
            }
            // reset pr_curr for the next iteration
            pr_next[v] = 0.0;
        }

    
    }
    std::printf("%d, %u, %lf\n",world_rank, num_edges, communicate_time_taken);

    PageRankType local_sum_of_page_ranks = 0, total_page_ranks=0;

    for (uintV u = startIndex; u < endIndex; u++)
    {
        local_sum_of_page_ranks += pr_curr[u];
    }
    //std::printf("process=%d, Sum of page rank : %ld\n", world_rank, local_count);

    // Get and sum all the processes
    if(world_rank != ROOT_PROCESS)
    {
        MPI_Reduce(&local_sum_of_page_ranks, nullptr, 1, PAGERANK_MPI_TYPE, MPI_SUM, 0, MPI_COMM_WORLD); 
    }else{
        MPI_Reduce(&local_sum_of_page_ranks, &total_page_ranks, 1, PAGERANK_MPI_TYPE, MPI_SUM, 0, MPI_COMM_WORLD);
    }

    if(world_rank == ROOT_PROCESS)
    {
        total_time_taken=total_timer.stop();
        std::printf("Sum of page rank : " PR_FMT "\n", total_page_ranks);
        std::printf("Time taken (in seconds) : %f\n", total_time_taken);
    }

    //time_taken = t1.stop();
    //std::printf("Time taken (in seconds) : %f\n", time_taken);
    delete[] pr_curr;
    delete[] pr_next;
    delete[] startIndexPtr;
    delete[] endIndexPtr;
    delete[] pr_temp;
}