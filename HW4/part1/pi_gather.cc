#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <sys/types.h>
#include <unistd.h>
#include <cmath>

const int MAX = std::sqrt(RAND_MAX);

int main(int argc, char **argv)
{
    // --- DON'T TOUCH ---
    MPI_Init(&argc, &argv);
    double start_time = MPI_Wtime();
    double pi_result;
    double *pi_results;
    long long int tosses = atoi(argv[1]);
    int world_rank, world_size;
    // ---

    // TODO: MPI init
	MPI_Comm_size(MPI_COMM_WORLD, &world_size);
	MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    // TODO: use MPI_Gather
	MPI_Scatter(nullptr, 0, MPI_DOUBLE, &pi_result, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	pi_result = 0;
	unsigned int seed = world_rank;//time(NULL);
	for ( long long int toss = tosses / world_size; toss > 0; toss--) {
		int r = rand_r(&seed);
		double x = (double(r % MAX) / MAX - 0.5);
		double y = (double(r / MAX) / MAX - 0.5);
		
		double distance_squared = x * x + y * y;
		if ( distance_squared <= 0.250000)
			pi_result++;
	}
    if (world_rank == 0)
	{
		pi_results = (double *)malloc(sizeof(double) * world_size);
	}
	MPI_Gather(&pi_result, 1, MPI_DOUBLE, pi_results, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	
    if (world_rank == 0)
    {
        // TODO: PI result
		pi_result = 0;
		for ( int i = 0; i < world_size; i++)
		{
			pi_result += pi_results[i];
		}
		pi_result = pi_result * 4.0 / tosses;

        // --- DON'T TOUCH ---
        double end_time = MPI_Wtime();
        printf("%lf\n", pi_result);
        printf("MPI running time: %lf Seconds\n", end_time - start_time);
        // ---
    }
    
    MPI_Finalize();
    return 0;
}
