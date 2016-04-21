#include "parallel.h"
#include "parallel_task.h"
#include "seq.h"


#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int main (int argc, char **argv)
{
	printf("FORMAT 128 - ITERATION 10\n\n");

	printf("Sequentiel\n\n");
	seq(argc,argv,1);

	printf("Parallel_for\n\n");
	parallel(argc,argv,1);

	printf("parallel_task\n\n");	
	parallel_task(argc,argv,1);
	return 0;
}