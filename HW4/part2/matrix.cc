#include <mpi.h>
#include <cstdio>
#include <iostream>
#include <stdio.h>

void construct_matrices(int *n_ptr, int *m_ptr, int *l_ptr, int **a_mat_ptr, int **b_mat_ptr)
{
    int world_rank, world_size;
	MPI_Comm_size(MPI_COMM_WORLD, &world_size);
	MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
	world_size -= 1;
	
	if (world_rank == 0)
	{
		std::cin >> *(n_ptr) >> *(m_ptr) >> *(l_ptr);
		
		int sizeA = (*n_ptr) * (*m_ptr);
		int sizeB = (*m_ptr) * (*l_ptr);
		int n = *n_ptr, m = *m_ptr, l = *l_ptr;
		*a_mat_ptr = new int[n*m]();
		*b_mat_ptr = new int[m*l]();
		//input shape
		
		for (int i = 0; i < sizeA; i++)
		{
			std::cin >> (*a_mat_ptr)[i];
		}
		for (int i = 0; i < sizeB; i++)
		{
			std::cin >> (*b_mat_ptr)[i];
		}
		
		/*for (int n = 0; n < *(n_ptr); n++)
		{
			a_mat_ptr[n] = new int[(*m_ptr)]();
			for (int m = 0; m < *(m_ptr); m++)
			{
				std::cin >> a_mat_ptr[n][m];
			}
		}
		for (int m = 0; m < *(m_ptr); m++)
		{
			b_mat_ptr[m] = new int[(*l_ptr)]();
			for (int l = 0; l < *(l_ptr); l++)
			{
				std::cin >> b_mat_ptr[m][l];
			}	
		}*/
	}
}

// Just matrix multiplication (your should output the result in this function)
// 
// n:     row number of matrix a
// m:     col number of matrix a / row number of matrix b
// l:     col number of matrix b
// a_mat: a continuous memory placing n * m elements of int
// b_mat: a continuous memory placing m * l elements of int
void matrix_multiply(const int n, const int m, const int l, const int *a_mat, const int *b_mat)
{
    int world_rank, world_size;
	MPI_Comm_size(MPI_COMM_WORLD, &world_size);
	MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
	
	int N, M, L;
	
	//SEND N, M, L
	if (world_rank == 0)
	{
		N = n;
		M = m;
		L = l;
		for (int i = 1; i < world_size; i++)
		{
			MPI_Send(&n, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
			MPI_Send(&m, 1, MPI_INT, i, 1, MPI_COMM_WORLD);
			MPI_Send(&l, 1, MPI_INT, i, 2, MPI_COMM_WORLD);
		}
	}
	else
	{
		MPI_Recv(&N, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		MPI_Recv(&M, 1, MPI_INT, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		MPI_Recv(&L, 1, MPI_INT, 0, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	}
	//printf("rank=%d => N=%d,M=%d,L=%d\n",world_rank,N,M,L);
	if (world_rank == 0)
	{
		//SEND DATA
		int stepN = (world_size-1);
		
		const int *a_first = a_mat;
		const int *b_first = b_mat;
		
		for (int in = 0; in < N; )
		{
			for (int i = 1; i < world_size && in < N; i++)
			{
				//SEND A
				//tag == row of N
				//rank 0 -> i+1
				MPI_Send(a_first, M, MPI_INT, i, in, MPI_COMM_WORLD);
				in++;
				a_first += M;
			}
		}
		for (int i = 1; i < world_size; i++)
		{
			//SEND B
			MPI_Send(b_mat, M*L, MPI_INT, i, 0, MPI_COMM_WORLD);
		}
		
		bool enterFirst = false;
		//RECV RESULT
		for (int in = 0; in < N; )
		{
			for (int i = 1; i < world_size && in < N; i++)
			{
				if (enterFirst)
				{
					//std::cout << "\n";
				}
				enterFirst = true;
				int *recv = new int[L]();
				MPI_Recv(recv, L, MPI_INT, i, in, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				//printf("RECV ROW %d FROM %d Data=(%d)\n", in, i, recv[0]);
				for (int il = 0; il < L; il++)
				{
					std::cout << recv[il] << " ";
				}
					std::cout << "\n";
				in++;
				delete[] recv;
			}
		}
	}
	else
	{
		int stepN = (world_size-1);
		
		int dataASize = (N + stepN - 1) / stepN * M;
		int dataCSize = (N + stepN - 1) / stepN * L;
		
		int *a_datas, *b_datas, *c_datas;
		a_datas = new int[dataASize];
		b_datas = new int[M*L];
		c_datas = new int[dataCSize]();
		
		//WORKER FOR RECV A AND B
		for (int in = world_rank - 1, idxA = 0; in < N; in += stepN, idxA += M)
		{
			MPI_Recv(&(a_datas[idxA]), M, MPI_INT, 0, in, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		}
		MPI_Recv(b_datas, M*L, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		
		
		
		//WORKER FOR CALCULATE
		int idxC = 0;
		int idxCFirst = 0;
		int idxAFirst = 0;
		for (int in = world_rank - 1; in < N; in += stepN)
		{
			for (int il = 0; il < L; il++)
			{
				int idxA = idxAFirst;
				for (int im = 0; im < M; im++)
				{
					int idxB = im * L + il;
					c_datas[idxC] += (a_datas[idxA] * b_datas[idxB]);
					//printf("rank=%d=====>%d<-%d*%d\n", world_rank, idxC, idxA, idxB);
					idxA++;
				}
				idxC++;
			}
			//SEND ROW RESULT
			//printf("SEND ROW %d FROM %d Data=(%d)\n", in, world_rank, c_datas[idxCFirst]);
			MPI_Send(&(c_datas[idxCFirst]), L, MPI_INT, 0, in, MPI_COMM_WORLD);
			idxCFirst += L;
			idxAFirst += M;
		}
		

		
	}
}
// Remember to release your allocated memory
void destruct_matrices(int *a_mat, int *b_mat)
{
    /*
    int world_rank, world_size;
	MPI_Comm_size(MPI_COMM_WORLD, &world_size);
	MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
	world_size -= 1;
	if (world_rank == world_size)
	{
		if (a_mat != nullptr)
		{	
			delete[] a_mat;
			a_mat = nullptr;
		}
		if (b_mat != nullptr)
		{	
			delete[] b_mat;
			b_mat = nullptr;
		}
	}*/
}