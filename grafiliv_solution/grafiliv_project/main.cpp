#include <stdio.h>
#include "kernel.cuh"
#include "render.h"

//CUDA includes
#include "cuda_runtime.h"
#include "device_launch_parameters.h"

int main(int argc, char **argv)
{
	const int arraySize = 5;
	const int a[arraySize] = { 1, 2, 3, 4, 5 };
	const int b[arraySize] = { 10, 20, 30, 40, 50 };
	int c[arraySize] = { 0 };

	// Add vectors in parallel.
	addWithCuda(c, a, b, arraySize);

	printf("{1,2,3,4,5} + {10,20,30,40,50} = {%d,%d,%d,%d,%d}\n",
		c[0], c[1], c[2], c[3], c[4]);

	glInit(argc, argv);

	cudaDeviceReset();
	return 0;
}