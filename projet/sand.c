
#include "parallel.h"
#include "parallel_task.h"
#include "seq.h"


#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int main (int argc, char **argv)
{

	//	parallel_task(argc,argv,1);
	seq(argc,argv,1);
	return 0;

	if(argc<2){
		seq(argc,argv,1);
		return 0;
	}



	if(!strcmp(argv[1], "homogeneous") || !strcmp(argv[1], "h")){

		if(!strcmp(argv[2], "sequential") || !strcmp(argv[2], "s")){
		  seq(argc,argv,1);
		}
		else if(!strcmp(argv[2], "parallel") || !strcmp(argv[2], "p")){
		  parallel(argc,argv,1);
		}
		else{
			  printf("usage : ./sand <INITIALIZATION> <ALGORITHM> \n\t-INITIALIZATION can be :\n\t\t-homogeneous, h : it starts the homogeneous case ;\n\t\t-centered, c : it starts the centered case ;\n\t-ALGORITHM can be :\n\t\t-sequential, s : it runs the sequential method ;\n\t\t-parallel, p : it runs the parallel method ;");
		}
	}
  
	else if(!strcmp(argv[1], "centered") || !strcmp(argv[1], "c")){
		if(!strcmp(argv[2], "sequential") || !strcmp(argv[2], "s")){
		  seq(argc,argv,0);
		}
		else if(!strcmp(argv[2], "parallel") || !strcmp(argv[2], "p")){
		  parallel(argc,argv,0);
		}
		else{
			  printf("usage : ./sand <INITIALIZATION> <ALGORITHM> \n\t-INITIALIZATION can be :\n\t\t-homogeneous, h : it starts the homogeneous case ;\n\t\t-centered, c : it starts the centered case ;\n\t-ALGORITHM can be :\n\t\t-sequential, s : it runs the sequential method ;\n\t\t-parallel, p : it runs the parallel method ;");
		  }
	}
  
  else{
		printf("usage : ./sand <INITIALIZATION> <ALGORITHM> \n\t-INITIALIZATION can be :\n\t\t-homogeneous, h : it starts the homogeneous case ;\n\t\t-centered, c : it starts the centered case ;\n\t-ALGORITHM can be :\n\t\t-sequential, s : it runs the sequential method ;\n\t\t-parallel, p : it runs the parallel method ;");
	}
	return 0;
}