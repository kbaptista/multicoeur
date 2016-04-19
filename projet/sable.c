

#include <stdio.h>
#include <stdlib.h>
#include <string.h>



int main (int argc, char **argv)
{

	unsigned function;
	unsigned configuration;

	if(argc<2)
	{
		printf("No argument, please use help argument to have indications.\n");
		sand_init_center();
	}
	else if (strcmp(argv[1],"sequentiel") || strcmp(argv[1],"sequentiel"))
	{
		function = 1;
	}
	else if (strcmp(argv[1],"parallel") || strcmp(argv[1],"parallel"))
	{
		function = 2;
	}
	if(strcmp(argv[1],"centered") || strcmp(argv[1],"centered"))
	{
		configuration = 1;
	}
	else if(strcmp(argv[1],"homogeneous") || strcmp(argv[1],"homogeneous"))
	{
		configuration = 2;
	}

	else if (atoi(argv[1])==2)
	{
	printf("Homogeneous Case\n");
	sand_init_homogeneous();
	}
	else
	{
	printf("usage : ./programme <number> where number is 1 or 2 \n -1 start the centered case \n -2 start the homogeneous case \n - no argument start the center case");
	return 1;
	}
}
