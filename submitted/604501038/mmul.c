#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

// Please modify this function


void mmul2(float* A, float* B, float* C, int n, int ni, int nj);

void mmul(float *A, float *B, float *C, int n)
{
	int tmp_pnum, tmp_pid;
	int i, j, k;
	int bsize;
	MPI_Comm_size(MPI_COMM_WORLD, &tmp_pnum);
	MPI_Comm_rank(MPI_COMM_WORLD, &tmp_pid);
	bsize = n/tmp_pnum;

	float* tmp;
	if (tmp_pid != 0) {
		// A = (float*)malloc( sizeof(float) * n * n );
		B = (float*)malloc( sizeof(float) * n * n );
	}
	tmp = (float*)malloc( sizeof(float) * n * bsize );	

	float* localA = (float*)malloc( sizeof(float) * n * bsize);

	MPI_Bcast(B, n*n, MPI_FLOAT, 0, MPI_COMM_WORLD);
	// MPI_Scatter(A, n*n/tmp_pnum, MPI_FLOAT, &A[from*n], n*n/tmp_pnum, MPI_FLOAT, 0, MPI_COMM_WORLD);

	MPI_Scatter(A, bsize*n, MPI_FLOAT, localA, bsize*n, MPI_FLOAT, 0, MPI_COMM_WORLD);
	
	mmul2(localA, B, tmp, bsize, n, n);

	// MPI_Gather(&C[from*n], n*n/tmp_pnum, MPI_FLOAT, C, n*n/tmp_pnum, MPI_FLOAT, 0, MPI_COMM_WORLD);
	for (i=0; i<bsize; i++){
		MPI_Gather(tmp, bsize*n, MPI_FLOAT, C, bsize*n, MPI_FLOAT, 0, MPI_COMM_WORLD);
	}

	if (tmp_pid != 0) {
		free(B);		
	}
	free(localA);
	free(tmp);
}



void mmul2(float* A, float* B, float* C, int ni, int nk, int nj)
{
	int i,j,k;
	int bk =64;
	int bi =bk;
	int bj =bk;


	for (i=0; i<ni/bi; i++) {
		for (j=0; j<nj/bj; j++) {
			float lc[bi][bj];
			int ii;
			int jj;
			
			for (ii=0; ii<bi; ii++) {
				for (jj=0; jj<bj; jj++) {
					lc[ii][jj] = 0;
				}
			}
			
			for (k=0; k<nk; k++) {
				for (ii=0; ii<bi; ii++) {
					for (jj=0; jj<bj; jj++) {
						lc[ii][jj] += A[nk*(i*bi+ii)+k]*B[k*nj+j*bj+jj];
					}
				}
			}
			
			for (ii=0; ii<bi; ii++) {
				for (jj=0; jj<bj; jj++) {
					C[nj*(i*bi+ii)+j*bj+jj] = lc[ii][jj];
				}
			}
			
		}
	}
	
}



