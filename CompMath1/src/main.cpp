

#include <iostream>
#include <iomanip>
#include <fstream>
#include <mpi.h>
#include <unistd.h>
#include <cmath>

double acceleration(double t)
{
  return sin(t);
}

void iterativeCalcRes(double* local_arrPoints, double a_i, double dt, double v_i, uint32_t num_points, double local_t0) {
    local_arrPoints[0] = a_i;
    local_arrPoints[1] = a_i + dt * v_i;
    for (uint32_t i = 2; i < num_points; i++)
    {
        local_arrPoints[i] = dt * dt * acceleration (local_t0 + (i - 1) * dt)
                    + 2 * local_arrPoints[i - 1] - local_arrPoints[i - 2];
    }
}

void calc (double* trace, uint32_t traceSize, double t0, double dt, double y0, double y1, int rank, int size)
{
    
    MPI_Status status;
    double v0, a_i, v_i, b_i, u_i, a_prev, v_prev, tau;
    uint32_t first_point, last_point, num_points;
    double* local_arrPoints;

    MPI_Bcast (&traceSize, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast (&t0, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast (&dt, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    
    first_point = traceSize * rank / size;
    last_point = traceSize * (rank + 1) / size;
    num_points = last_point - first_point;
    local_arrPoints = (double *)calloc(num_points, sizeof (*local_arrPoints));
    v0 = 0.0;
    v_i = 0.0;
    b_i = 0.0;
    u_i = 0.0;
    a_prev = 0.0;
    v_prev = 0.0;
    tau = dt * traceSize / size;
    t0 = t0 + tau * rank;
    
    if (rank == 0 && size > 0) {
        a_i = y0;
    } else {
        a_i = 0.0;
    }
    
    //first step
    iterativeCalcRes(local_arrPoints, a_i, dt, v_i, num_points, t0);
    
    b_i = local_arrPoints[num_points - 1];
    u_i = (local_arrPoints[num_points - 1] - local_arrPoints[num_points - 2]) / dt;
    
    if (size > 1) {
        
        if (rank != size - 1 && size > 1) {
            MPI_Send (&u_i, 1, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD);
            MPI_Send (&b_i, 1, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD);
        }
        
        if (rank != 0 && size > 1) {
            MPI_Recv (&v_prev, 1, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD, &status);
            MPI_Recv (&a_prev, 1, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD, &status);
            
            v_i = v_prev;
            u_i += v_prev;
            
            a_i = a_prev;
            b_i += a_prev + v_prev * tau;
        }

        if (rank == size - 1 && size > 1) {
            MPI_Send (&b_i, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
        }
        
        if (rank == 0 && size > 1) {
            MPI_Recv (&a_prev, 1, MPI_DOUBLE, size - 1, 0, MPI_COMM_WORLD, &status);
            v0 = (y1 - a_prev) / (dt * traceSize);
        }
        
        MPI_Bcast (&v0, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        a_i += v0 * tau * rank;
        v_i += v0;
    }
    
    if (size == 1) {
        v0 = (y1 - b_i) / (dt * traceSize);
        v_i = v0;
    }

    // The final step
    if (rank == 0) {
        iterativeCalcRes(trace, a_i, dt, v_i, num_points, t0);
        
        for (int i = 1; i < size; i++)
        {
            uint32_t first = traceSize * i / size;
            uint32_t last = traceSize * (i + 1) / size;
            uint32_t len = last - first;
            MPI_Recv (trace + first, len, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, &status);
        }
    } else {
        iterativeCalcRes(local_arrPoints, a_i, dt, v_i, num_points, t0);
        MPI_Send (local_arrPoints, num_points, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
    }
     
    free (local_arrPoints);
}

int main(int argc, char** argv)
{
  int rank = 0, size = 0, status = 0;
  uint32_t traceSize = 0;
  double t0 = 0, t1 = 0, dt = 0, y0 = 0, y1 = 0;
  double* trace = 0;

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  if (rank == 0)
  {
    // Check arguments
    if (argc != 3)
    {
      std::cout << "[Error] Usage <inputfile> <output file>\n";
      status = 1;
      MPI_Bcast(&status, 1, MPI_INT, 0, MPI_COMM_WORLD);
      return 1;
    }

    // Prepare input file
    std::ifstream input(argv[1]);
    if (!input.is_open())
    {
      std::cout << "[Error] Can't open " << argv[1] << " for write\n";
      status = 1;
      MPI_Bcast(&status, 1, MPI_INT, 0, MPI_COMM_WORLD);
      return 1;
    }

    // Read arguments from input
    input >> t0 >> t1 >> dt >> y0 >> y1;
    MPI_Bcast(&status, 1, MPI_INT, 0, MPI_COMM_WORLD);
    traceSize = (t1 - t0)/dt;
    trace = new double[traceSize];

    input.close();
  } else {
    MPI_Bcast(&status, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if (status != 0)
    {
      return 1;
    }
  }

  calc(trace, traceSize, t0, dt, y0, y1, rank, size);

  if (rank == 0)
  {
    // Prepare output file
    std::ofstream output(argv[2]);
    if (!output.is_open())
    {
      std::cout << "[Error] Can't open " << argv[2] << " for read\n";
      delete trace;
      return 1;
    }

    for (uint32_t i = 0; i < traceSize; i++)
    {
      output << " " << trace[i];
    }
    output << std::endl;
    output.close();
    delete trace;
  }

  MPI_Finalize();
  return 0;
}
