__kernel void addmat(__global unsigned *A,
		     __global float *C)
{

	int col = get_global_id(0); //j
	int row = get_global_id(1); //i
	int N = get_global_size(0);

	//C[i][j]
	C[row*N+col] = A[row*N+col]+B[row*N+col];

}
