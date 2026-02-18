#include <mpi.h>
#include <stdio.h>

int main(int argc, char** argv) {
  printf("Before MPI_Init\n");
  int status = MPI_Init(&argc, &argv);
  printf("After MPI_Init, status = %d\n", status);
  MPI_Finalize();
  printf("After MPI_Finalize\n");
  return 0;
}
