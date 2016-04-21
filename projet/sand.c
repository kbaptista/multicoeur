
#include "parallel.h"
#include "seq.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int main (int argc, char **argv)
{
  char * man = "usage : ./sand <INITIALIZATION> <ALGORITHM> \n\t-INITIALIZATION can be :\n\t\t-homogeneous, h or 0 : it starts the homogeneous case ;\n\t\t-centered, c or 1 : it starts the centered case ;\n\t-ALGORITHM can be :\n\t\t-sequential, s or 0 : it runs the sequential method ;\n\t\t-parallel, p or 1 : it runs the parallel method ;";


	if(argc<3){
		seq(argc,argv,1);
		return 0;
	}

  if(!strcmp(argv[1], "homogeneous") || !strcmp(argv[1], "h") || atoi(argv[1])==0){
    if(!strcmp(argv[2], "sequential") || !strcmp(argv[2], "s") || atoi(argv[2])==0){
      seq(argc,argv,1);
    }
    else if(!strcmp(argv[2], "for") || !strcmp(argv[2], "f") || atoi(argv[2])==1){
      parallel(argc,argv,1);
    }
    else if(!strcmp(argv[2], "task") || !strcmp(argv[2], "t") || atoi(argv[2])==3){
      parallel_task(argc,argv,1);
    }
    
    else{
		  printf("%s",man);
	  }
  }
  
  else if(!strcmp(argv[1], "centered") || !strcmp(argv[1], "c") || atoi(argv[1])==1){
    if(!strcmp(argv[2], "sequential") || !strcmp(argv[2], "s") || atoi(argv[2])==0){
      seq(argc,argv,0);
    }
    else if(!strcmp(argv[2], "parallel") || !strcmp(argv[2], "p") || atoi(argv[2])==1){
      parallel(argc,argv,0);
    }
    else if(!strcmp(argv[2], "task") || !strcmp(argv[2], "t") || atoi(argv[2])==3){
      parallel_task(argc,argv,0);
    }
    else{
		  printf("%s",man);
    }
  }
  
  else{
		  printf("%s",man);
	}
	return 0;
}