#include "display.h"
#include "treatment.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>


// ------------------------------------------------------------------------------
// -------------        Fonctions de traitement des Arguments       -------------
// ------------------------------------------------------------------------------

/** 
  * initialise matrices : 
  * 0 means centering configuration
  * 1 means homogeneous configuration
*/
void sand_init (int init)
{
  if(!init)
  {
    cb_sand_init_center();
  }
  else
  {
    cb_sand_init_homogeneous();
  }
}

char * man = "usage : ./sand <INITIALIZATION> <SIZE> <ALGORITHM> \n\n\t-INITIALIZATION can be :\n\t\t-homogeneous or h : it starts the homogeneous case ;\n\t\t-centered or c : it starts the centered case ;\n\t-SIZE can be :\n\t\t-128 ;\n\t\t-512 ;\n\t-ALGORITHM can be :\n\t\t-sequential or s : it runs the sequential method ;\n\t\t-parallel or p : it runs the parallel method ;\n\t\t-task or t : it runs the parallel task method ;\n";

/* Add your treatment choices here */
/*
void treatment_compare(int argc, char **argv, int init_compare_res, int size_compare_res){

  if(!strcmp(argv[3], "sequential") || !strcmp(argv[3], "s")){
    	printf("Treatment Sequential\n\n");
      sand_init(init_compare_res);
      cb_seq(argc, argv);
    }
    else if(!strcmp(argv[3], "for") || !strcmp(argv[3], "f")){
    	printf("Treatment Parallel For\n\n");
      sand_init(init_compare_res);
      cb_parallel(argc, argv);
    }
    else if(!strcmp(argv[3], "task") || !strcmp(argv[3], "t")){
    	printf("Treatment Parallel Task\n\n");
      sand_init(init_compare_res);
      cb_parallel_task(argc, argv);
    }
    
    else{
		printf("%s",man);
		return;
	}
}

void size_compare(int argc, char **argv, int init_compare_res){
	if (atoi(argv[2])==128)
	{
		treatment_compare(argc,argv,init_compare_res, 128);
	}
	else if (atoi(argv[2])==512)
	{
		treatment_compare(argc, argv, init_compare_res, 512);
	}
  	else{
  		printf("%s",man);
  		return;
  	}
}

void init_compare(int argc, char **argv){

  if(!strcmp(argv[1], "homogeneous") || !strcmp(argv[1], "h")){
  	printf("Case Homogeneous\n");
  	size_compare(argc, argv, 1);
  }
  else if(!strcmp(argv[1], "centered") || !strcmp(argv[1], "c")){
  	printf("Case Centered\n");
  	size_compare(argc, argv, 0);
  }
  else{
  	printf("%s",man);
  	return;
  }
}
*/
// ------------------------------------------------------------------------------
// -------------                         MAIN                       -------------
// ------------------------------------------------------------------------------

int main (int argc, char **argv)
{

	if(argc<2){
	  printf("%s", man);
	  /*
    sand_init(1);
		cb_seq(argc,argv);
		return 0;*/
	}
	
	treatment(argc,argv);

	//init_compare(argc, argv);

	return 0;
}