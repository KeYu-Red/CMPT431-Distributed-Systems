#include "core/utils.h"
#include <iomanip>
#include <iostream>
#include <stdlib.h>
#include <vector>
#include <thread>

#define sqr(x) ((x) * (x))
#define DEFAULT_NUMBER_OF_POINTS "12345678"

std::atomic<int> current_thread_id(0);
std::atomic<int> total_circle_points(0);

uint c_const = (uint)RAND_MAX + (uint)1;
inline double get_random_coordinate(uint *random_seed) {
  return ((double)rand_r(random_seed)) / c_const;
}

uint get_points_in_circle(uint n, uint random_seed) {
  uint circle_count = 0;
  double x_coord, y_coord;
  for (uint i = 0; i < n; i++) {
    x_coord = (2.0 * get_random_coordinate(&random_seed)) - 1.0;
    y_coord = (2.0 * get_random_coordinate(&random_seed)) - 1.0;
    if ((sqr(x_coord) + sqr(y_coord)) <= 1.0)
      circle_count++;
  }
  return circle_count;
}

void singleThreadPiCalculation(uint n)
{
  timer serial_timer;
  double time_taken = 0.0;
  uint random_seed = 1;

  serial_timer.start();
  // Create threads and distribute the work across T threads
  
  // -------------------------------------------------------------------
  uint circle_points = get_points_in_circle(n, random_seed);
  //double pi_value = 4.0 * (double)circle_points / (double)n;
  // -------------------------------------------------------------------
  time_taken = serial_timer.stop();

//std::cout << "thread_id, points_generated, circle_points, time_taken\n";
  std::cout << current_thread_id++ << ", "<< n << ", " << circle_points << ", " << time_taken<< "\n";
  total_circle_points += circle_points;
}

void piCalculation(uint n, uint n_threads) {
  timer total_serial_timer;
  total_serial_timer.start();
  double time_taken = 0.0;
  uint random_seed = 1;

  std::vector<std::thread> vec_threads;
  uint single_n = n/n_threads;
  std::cout << "thread_id, points_generated, circle_points, time_taken\n";
  for(int i=0; i<n_threads; i++)
  {
    if(i==n_threads-1)
      single_n = n - single_n*i;
    vec_threads.push_back(std::thread(singleThreadPiCalculation, single_n));
  }
  for(int i=0; i<n_threads; i++)
    vec_threads[i].join();

  double pi_value = 4.0 * (double)total_circle_points / (double)n;
  time_taken = total_serial_timer.stop();
  // Print the overall statistics
  std::cout << "Total points generated : " << n << "\n";
  std::cout << "Total points in circle : " << total_circle_points << "\n";
  std::cout << "Result : " << std::setprecision(VAL_PRECISION) << pi_value
            << "\n";
  std::cout << "Time taken (in seconds) : " << std::setprecision(2)
            << time_taken << "\n";
}

int main(int argc, char *argv[]) {
  // Initialize command line arguments
  cxxopts::Options options("pi_calculation",
                           "Calculate pi using serial and parallel execution");
  options.add_options(
      "custom",
      {
          {"nPoints", "Number of points",
           cxxopts::value<uint>()->default_value(DEFAULT_NUMBER_OF_POINTS)},
          {"nWorkers", "Number of workers",
           cxxopts::value<uint>()->default_value(DEFAULT_NUMBER_OF_WORKERS)},
      });

  auto cl_options = options.parse(argc, argv);
  uint n_points = cl_options["nPoints"].as<uint>();
  uint n_workers = cl_options["nWorkers"].as<uint>();
  std::cout << std::fixed;
  std::cout << "Number of points : " << n_points << "\n";
  std::cout << "Number of workers : " << n_workers << "\n";

  piCalculation(n_points, n_workers);

  return 0;
}
