#include <stdio.h>
#include <math.h>
#include <stdlib.h>

// number of threads per block
//const int numThreadsPerBlock = 256;
const int numThreadsPerBlock = 1024;

//// device to use in case there is more than one
//static int selectedDevice = 0;

__global__ void kern_Dvxv(const int n, double* v1, double* v2, double* v3)
{
	int i = blockIdx.x * blockDim.x + threadIdx.x;
	while (i < n)
	{
		v3[i] = v1[i] * v2[i];
		i += blockDim.x * gridDim.x;
	}
}


extern "C" {

	//element-wise vector multiplication double precision
	void cuda_Dvxv(const int* n, double* v1, double* v2, double* v3)
	{
		int numBlocks = (*n+(numThreadsPerBlock-1)) / numThreadsPerBlock;

		kern_Dvxv<<<numBlocks,numThreadsPerBlock>>>(*n,v1,v2,v3);

		return;
	}

}

