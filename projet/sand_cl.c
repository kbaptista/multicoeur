
#include <stdio.h>
#include <string.h>
#include <stdbool.h>
#include <sys/time.h>
#ifdef __APPLE__
#  include <OpenCL/opencl.h>
#else
#  include <CL/opencl.h>
#endif
#include "util.h"

#define MAX_PLATFORMS 3
#define MAX_DEVICES   5

#define ITERATIONS    100

// Application-specific data
//
#define KERNEL_NAME  "sand"


unsigned SIZE = 64;
unsigned TILE = 16;
unsigned TILE2 = 1;

unsigned *ocean;

cl_mem ocean_buffer;
cl_mem output_buffer;

static bool is_end = false ;

unsigned get (unsigned x, unsigned y)
{
  return ocean[x*SIZE+y];
}

static void print()
{
  for (int x = 0; x < SIZE; x++)
  {
    printf("[%3d]",x );
    for (int y = 0; y < SIZE; y++) 
    {
      printf(" %3d ",ocean[x*SIZE+y]);
    }
    printf("\n");
  }
    printf("     ");
  for (int i = 0; i < SIZE; ++i)
  {
    printf("[%3d]",i );
  }
    printf("\n");
}

void sand_init_homogeneous()
{
  for (int x = 0; x < SIZE; x++)
  {
    for (int y = 0; y < SIZE; y++)
    {
      if(x > 0 && x < SIZE-1 && y > 0 && y < SIZE-1)
        ocean[x*SIZE+y] = 5;
      else
        
        ocean[x*SIZE+y] = 0;
    }
  }
}

/** 
  * initialize matrices with no one grain except in the middle case which have 100000 grain.
*/
void sand_init_center()
{
  int center = 100000;
  for (int x = 0; x < SIZE; x++)
  {
    for (int y = 0; y < SIZE; y++)
    {
      ocean[x*SIZE+y] = 0;
    }
  }
  ocean[SIZE*SIZE/2+SIZE/2] = center;
}
static void alloc_buffers_and_user_data(cl_context context)
{
  // CPU side
  ocean = malloc(SIZE * SIZE * sizeof(unsigned));

  sand_init_center();

  // Allocate buffers inside device memory
  //
  ocean_buffer = clCreateBuffer(context,  CL_MEM_READ_WRITE,  sizeof(unsigned) * SIZE * SIZE, NULL, NULL);
  if (!ocean_buffer)
    error("Failed to allocate input buffer A");
  output_buffer = clCreateBuffer(context,  CL_MEM_READ_WRITE,  sizeof(unsigned) * SIZE * SIZE, NULL, NULL);
  if (!output_buffer)
    error("Failed to allocate input buffer A");
}

static void check_output_data(void)
{
  //TODO vérification temporaire
  unsigned correct = 0;
  int s = (SIZE-1)*(SIZE-1);
  for(int x = 1; x < s ; x++)
    if (((int(*))output_buffer)[x]<4)
      correct++;

  printf("\tComputed '%d/%d' correct values!\n", correct, (SIZE-1)*(SIZE-1));
}

static void free_buffers_and_user_data(void)
{
  free(ocean);

  clReleaseMemObject(ocean_buffer);
  clReleaseMemObject(output_buffer);
}

static void send_input(cl_command_queue queue)
{
  cl_int err;

  err = clEnqueueWriteBuffer(queue, ocean_buffer, CL_TRUE, 0,
			     sizeof(unsigned) * SIZE * SIZE, ocean, 0, NULL, NULL);
  check(err, "Failed to write to input array A");

  err = clEnqueueWriteBuffer(queue, output_buffer, CL_TRUE, 0,
           sizeof(unsigned) * SIZE * SIZE, ocean, 0, NULL, NULL);
}

static void retrieve_output(cl_command_queue queue)
{
  cl_int err;

  int s = SIZE * SIZE;
  err = clEnqueueReadBuffer(queue, output_buffer, CL_TRUE, 0,
			    sizeof(unsigned) * s, ocean, 0, NULL, NULL );  
  check(err, "Failed to read output array C");
}

int main(int argc, char** argv)
{
  cl_platform_id pf[MAX_PLATFORMS];
  cl_uint nb_platforms = 0;
  cl_int err;                            // error code returned from api calls
  cl_device_type device_type = CL_DEVICE_TYPE_ALL;

  // Filter args
  //
  argv++;
  while (argc > 1)
  {
    if(!strcmp(*argv, "-g") || !strcmp(*argv, "--gpu-only"))
    {
      if(device_type != CL_DEVICE_TYPE_ALL)
        error("--gpu-only and --cpu-only can not be specified at the same time\n");
      device_type = CL_DEVICE_TYPE_GPU;
    } 
    else if(!strcmp(*argv, "-c") || !strcmp(*argv, "--cpu-only"))
    {
      if(device_type != CL_DEVICE_TYPE_ALL)
        error("--gpu-only and --cpu-only can not be specified at the same time\n");
      device_type = CL_DEVICE_TYPE_CPU;
    }
    else if(!strcmp(*argv, "-s") || !strcmp(*argv, "--size"))
    {
      unsigned i;
      int r;
      char c;

      r = sscanf(argv[1], "%u%[mMkK]", &SIZE, &c);

      if (r == 2)
      {
        if (c == 'k' || c == 'K')
          SIZE *= 1024;
        else if (c == 'm' || c == 'M')
          SIZE *= 1024 * 1024;
      }

      argc--;
      argv++;
    }
    else
      break;
    argc--;
    argv++;
  }

  if(argc > 1)
    TILE = atoi(*argv);

  if(argc > 2)
    TILE2 = atoi(argv[1]);

  // Get list of OpenCL platforms detected
  //
  err = clGetPlatformIDs(3, pf, &nb_platforms);
  check(err, "Failed to get platform IDs");

  printf("%d OpenCL platforms detected\n", nb_platforms);

  // For each platform do
  //
  for (cl_int p = 0; p < nb_platforms; p++)
  {
    cl_uint num;
    int platform_valid = 1;
    char name[1024], vendor[1024];
    cl_device_id devices[MAX_DEVICES];
    cl_uint nb_devices = 0;
    cl_context context;                 // compute context
    cl_program program;                 // compute program
    cl_kernel kernel;

    err = clGetPlatformInfo(pf[p], CL_PLATFORM_NAME, 1024, name, NULL);
    check(err, "Failed to get Platform Info");

    err = clGetPlatformInfo(pf[p], CL_PLATFORM_VENDOR, 1024, vendor, NULL);
    check(err, "Failed to get Platform Info");

    printf("Platform %d: %s - %s\n", p, name, vendor);

    // Get list of devices
    //
    err = clGetDeviceIDs(pf[p], device_type, MAX_DEVICES, devices, &nb_devices);
    printf("nb devices = %d\n", nb_devices);

    if(nb_devices == 0)
      continue;

    // Create compute context with "device_type" devices
    //
    context = clCreateContext (0, nb_devices, devices, NULL, NULL, &err);
    check(err, "Failed to create compute context");


    // Load program source into memory
    //
    const char  *opencl_prog;
    opencl_prog = file_load(KERNEL_FILE);

    // Attach program source to context
    //
    program = clCreateProgramWithSource(context, 1, &opencl_prog, NULL, &err);
    check(err, "Failed to create program");

    // Compile program
    //
    {
      char flags[1024];

      sprintf (flags,
         "-cl-mad-enable -cl-fast-relaxed-math -DSIZE=%d -DTILE=%d -DTYPE=%s",
         SIZE, TILE, "float");

      err = clBuildProgram (program, 0, NULL, flags, NULL, NULL);
      if(err != CL_SUCCESS) {
        size_t len;

        // Display compiler log
        //
        clGetProgramBuildInfo(program, devices[0], CL_PROGRAM_BUILD_LOG, 0, NULL, &len);
        {
          char buffer[len+1];

          fprintf(stderr, "--- Compiler log ---\n");
          clGetProgramBuildInfo(program, devices[0], CL_PROGRAM_BUILD_LOG, sizeof(buffer), buffer, NULL);
          fprintf(stderr, "%s\n", buffer);
          fprintf(stderr, "--------------------\n");
        }
        if(err != CL_SUCCESS)
          error("Failed to build program!\n");
      }
    }

    // Create the compute kernel in the program we wish to run
    //
    kernel = clCreateKernel(program, KERNEL_NAME, &err);
    check(err, "Failed to create compute kernel");

    // Allocate and initialize input data
    //
    alloc_buffers_and_user_data(context);

    // Iterate over devices
    //
    for(cl_int dev = 0; dev < nb_devices; dev++)
    {
      cl_command_queue queue;

      char name[1024];
      cl_device_type dtype;

      err = clGetDeviceInfo(devices[dev], CL_DEVICE_NAME, 1024, name, NULL);
      check(err, "Cannot get type of device");
      err = clGetDeviceInfo(devices[dev], CL_DEVICE_TYPE, sizeof(cl_device_type), &dtype, NULL);
      check(err, "Cannot get type of device");

      printf("\tDevice %d : %s [%s]\n", dev, (dtype == CL_DEVICE_TYPE_GPU) ? "GPU" : "CPU", name);

      // Create a command queue
      //
      queue = clCreateCommandQueue(context, devices[dev], CL_QUEUE_PROFILING_ENABLE, &err);
      check(err,"Failed to create command queue");

      // Write our data set into device buffer
      //
      send_input(queue);

      // Execute kernel
      //
      //for (int i = 0; i < 3; ++i)
      {
      	cl_event prof_event;
      	cl_ulong start, end;
      	struct timeval t1,t2;
      	double timeInMicroseconds;
      	size_t global[1] = { (SIZE-2) * (SIZE-2)};  // global domain size for our calculation
      	size_t local[1]  = { 2 };   // local domain size for our calculation

      	printf("\t%d Threads in workgroups of %d\n", global[0], local[0]);

      	// Set kernel arguments
      	//
      	err = 0;
      	err  = clSetKernelArg(kernel, 0, sizeof(cl_mem), &ocean_buffer);
        err |= clSetKernelArg(kernel, 1, sizeof(cl_mem), &output_buffer);
      	check(err, "Failed to set kernel arguments");

      	gettimeofday (&t1, NULL);

      	// for (unsigned iter = 0; iter < ITERATIONS; iter++)
        // {
      	  err = clEnqueueNDRangeKernel(queue, kernel, 1, NULL, global, local,
      				       0, NULL, &prof_event);
      	  check(err, "Failed to execute kernel");
      	// }

      	// Wait for the command commands to get serviced before reading back results
      	//
      	clFinish(queue);

      	gettimeofday (&t2,NULL);

      	// Check performance
      	//
      	timeInMicroseconds = (double)TIME_DIFF(t1, t2); // ITERATIONS;

      	printf("\tComputation performed in %lf µs over device #%d\n", timeInMicroseconds,dev);

      	clReleaseEvent(prof_event);
      }

      // Read back the results from the device to verify the output
      //
      retrieve_output(queue);
      

      // Validate computation
      //
      check_output_data();

      clReleaseCommandQueue(queue);
    }
      print(); 

    // Cleanup
    //
    free_buffers_and_user_data();

    clReleaseKernel(kernel);
    clReleaseProgram(program);
    clReleaseContext(context);
  }

  return 0;
}
