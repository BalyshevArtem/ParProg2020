#include <iostream>
#include <iomanip>
#include <fstream>
#include <omp.h>
#include <cmath>

double func(double x)
{
  return sin(x);
}


double kehan_alg( double *values, int num_steps)
{
    double sum, c;
    double t, y;
    
    sum = 0.0;
    c = 0.0;

    for ( int i = 0; i < num_steps; i++) {
        y = values[i] - c;
        t = sum + y;
        c = (t - sum) - y;
        sum = t;
    }

    return sum;
}

double calc( double x0, double x1, double dx, uint32_t num_threads)
{
    int num_steps;
    double *values;
    double sum;
    
    /* Считаем кол-во шагов, необходимое для подсчета интеграла */
    num_steps = ( x1 - x0) / dx;
    
    /* Выделяем массив для хранения значения sin на каждом шаге */
    values = ( double*)calloc( num_steps, sizeof( double) );
    
    /* Считаем параллельно значение в num_threads потоках значение синуса в различных точках*/
    #pragma omp parallel for num_threads(num_threads)
    for ( int i = 0; i <= num_steps; i++) {
        values[i] =  func( x0 + i * dx);
    }
    
    for ( int i = 0; i < num_steps; i++) {
        values[i] = ( ( values[i] + values[i + 1]) * dx) / 2;
    }

    /* Алгоритм для более точного подсчета суммы */
    sum = kehan_alg( values, num_steps);

    free(values);
    return sum;
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
  double x0 = 0.0, x1 =0.0, dx = 0.0;
  uint32_t num_threads = 0;
  input >> x0 >> x1 >> dx >> num_threads;

  // Calculation
  double res = calc(x0, x1, dx, num_threads);

  // Write result
  output << std::setprecision(13) << res << std::endl;
  // Prepare to exit
  output.close();
  input.close();
  return 0;
}
