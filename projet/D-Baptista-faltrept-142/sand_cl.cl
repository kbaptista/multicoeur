__kernel void sand(	__global unsigned *ocean, 
					__global unsigned *output_buffer,
					__global int * res_buffer)
{

	int row = get_global_id(0)+1; //i
	int col = get_global_id(1)+1; //j

	output_buffer[col*SIZE+row] = 	ocean[col*SIZE+row]%4 +
									ocean[col*SIZE+row-1]/4 +
									ocean[col*SIZE+row+1]/4 +
									ocean[(col-1)*SIZE+row]/4 +
									ocean[(col+1)*SIZE+row]/4 ;

	__private tmp = output_buffer[col*SIZE+row] - ocean[col*SIZE+row];
	// no need to protect res_buffer since we dont care of the value in it
	// we just want to know if there are new values
	res_buffer += (tmp < 0) ? -tmp : tmp ; 
}
