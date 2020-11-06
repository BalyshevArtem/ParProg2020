#include <iostream>
#include <iomanip>
#include <fstream>
#include <mpi.h>
#include <unistd.h>
#include <cmath>

void mainCalc(double* arr, uint32_t begin, uint32_t end, uint32_t xSize) {
    for (uint32_t i = 0; i < (end - begin) * xSize; i++) {
        arr[i] = sin(0.00001 * arr[i]);
    }
}

void mainProcess(double* arr, uint32_t yLocSize, uint32_t ySize, uint32_t xSize, int size) {
    uint32_t begin, end;
    MPI_Status status;
    
    for (int i = 1; i < size; i++) {
        begin = i * yLocSize;
        end = (i + 1) * yLocSize;
        if (end > ySize) {
            end = ySize;
        }
        if (begin > ySize) {
            begin = ySize;
        }
        MPI_Send(&(arr[begin * xSize]), xSize * (end - begin), MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
    }
        
    begin = 0;
    if (yLocSize < ySize) {
        end = yLocSize;
    } else {
        end = ySize;
    }
        
    mainCalc(arr, begin, end, xSize);
        
    for (int i = 1; i < size; i++) {
        begin = i * yLocSize;
        end = (i + 1) * yLocSize;
        if (end > ySize)
            end = ySize;
        if (begin > ySize)
        begin = ySize;
        MPI_Recv (&(arr[begin * xSize]), xSize * (end - begin), MPI_DOUBLE, i, 0, MPI_COMM_WORLD, &status);
    }
}

void slaveProcesses(int rank, uint32_t yLocSize, uint32_t ySize, uint32_t xSize) {
    uint32_t begin, end;
    MPI_Status status;
    double *newArr;
    
    begin = rank * yLocSize;
    end = (rank + 1) * yLocSize;
    if (end > ySize)
        end = ySize;
    if (begin > ySize)
        begin = ySize;
    newArr = (double*) malloc (sizeof(double) * xSize * (end - begin));
           
    MPI_Recv(newArr, xSize * (end - begin), MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &status);
    mainCalc(newArr, begin, end, xSize);
    MPI_Send(newArr, xSize * (end - begin), MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
    free(newArr);
}

void calc(double* arr, uint32_t ySize, uint32_t xSize, int rank, int size)
{
    uint32_t yLocSize;
    
    MPI_Bcast(&ySize, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&xSize, 1, MPI_INT, 0, MPI_COMM_WORLD);

    yLocSize = ySize / size + 1;
    
    if (rank == 0 && size > 0) {
        mainProcess(arr, yLocSize, ySize, xSize, size);
    } else {
        slaveProcesses(rank, yLocSize, ySize, xSize);
    }
}


int main(int argc, char** argv)
{
  int rank = 0, size = 0, buf = 0;
  uint32_t ySize = 0, xSize = 0;
  double* arr = 0;

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  if (rank == 0)
  {
    // Check arguments
    if (argc != 3)
    {
      std::cout << "[Error] Usage <inputfile> <output file>\n";
      buf = 1;
      MPI_Bcast(&buf, 1, MPI_INT, 0, MPI_COMM_WORLD);
      return 1;
    }

    // Prepare input file
    std::ifstream input(argv[1]);
    if (!input.is_open())
    {
      std::cout << "[Error] Can't open " << argv[1] << " for write\n";
      buf = 1;
      MPI_Bcast(&buf, 1, MPI_INT, 0, MPI_COMM_WORLD);
      return 1;
    }

    // Read arguments from input
    input >> ySize >> xSize;
    MPI_Bcast(&buf, 1, MPI_INT, 0, MPI_COMM_WORLD);

    arr = new double[ySize * xSize];

    for (uint32_t y = 0; y < ySize; y++)
    {
     for (uint32_t x = 0; x < xSize; x++)
      {
        input >> arr[y*xSize + x];
      }
    }
    input.close();
  } else {
    MPI_Bcast(&buf, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if (buf != 0)
    {
      return 1;
    }
  }

  calc(arr, ySize, xSize, rank, size);

  if (rank == 0)
  {
    // Prepare output file
    std::ofstream output(argv[2]);
    if (!output.is_open())
    {
      std::cout << "[Error] Can't open " << argv[2] << " for read\n";
      delete arr;
      return 1;
    }
    for (uint32_t y = 0; y < ySize; y++)
    {
      for (uint32_t x = 0; x < xSize; x++)
      {
        output << " " << arr[y*xSize + x];
      }
      output << std::endl;
    }
    output.close();
    delete arr;
  }

  MPI_Finalize();
  return 0;
}
