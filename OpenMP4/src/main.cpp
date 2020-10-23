#include <iostream>
#include <iomanip>
#include <fstream>
#include <omp.h>

double calc( uint32_t x_last, uint32_t num_threads) {
    double *arr_currResults, *arr_currFact;
    int size_num, curr_thread, first_val, last_val;
    double currFact, result;
    
    arr_currResults = ( double*)malloc( sizeof( double) * num_threads);
    arr_currFact = ( double*)malloc( sizeof( double) * num_threads);
    size_num = x_last % num_threads == 0 ? x_last / num_threads : x_last / num_threads + 1;
    x_last--;
    result = 1;
    currFact = 1;

    #pragma omp parallel private(curr_thread, first_val, last_val) num_threads(num_threads)
    {
        curr_thread = omp_get_thread_num();
        first_val = size_num * curr_thread + 1;
        last_val = size_num * (curr_thread + 1);
        
        arr_currFact[curr_thread] = 1;
        arr_currResults[curr_thread] = 0;
        
        if (last_val > x_last) {
            last_val = x_last;
        }
        
        for (int i = first_val; i <= last_val; i++) {
            arr_currFact[curr_thread] /= i;
            arr_currResults[curr_thread] += arr_currFact[curr_thread];
        }
    }
    
    for (int i = 0; i < num_threads; i++) {
        result += currFact * arr_currResults[i];
        currFact *= arr_currFact[i];
    }

    return result;
}


int main(int argc, char** argv)
{
  // Check arguments
  if (argc != 3)
  {
    std::cout << "[Error] Usage <inputfile> <output file>\n";
    return 1;
  }

  // Prepare input file
  std::ifstream input(argv[1]);
  if (!input.is_open())
  {
    std::cout << "[Error] Can't open " << argv[1] << " for write\n";
    return 1;
  }

  // Prepare output file
  std::ofstream output(argv[2]);
  if (!output.is_open())
  {
    std::cout << "[Error] Can't open " << argv[2] << " for read\n";
    input.close();
    return 1;
  }

// Read arguments from input
  uint32_t x_last = 0, num_threads = 0;
  input >> x_last >> num_threads;

  // Calculation
  double res = calc(x_last, num_threads);

  // Write result
  output << std::setprecision(15) << res << std::endl;
  // Prepare to exit
  output.close();
  input.close();
  return 0;
}
