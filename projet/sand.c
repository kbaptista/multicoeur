#include "display.h"
#include "treatment.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

char * man = "usage : ./sand <GUI> <INITIALIZATION> <SIZE> <ALGORITHM> \n\n\t-GUI can be :\n\t\t-0 to disable GUI\n\t\t-1 to enable GUI\n\t-INITIALIZATION can be :\n\t\t-homogeneous or h : it starts the homogeneous case ;\n\t\t-centered or c : it starts the centered case ;\n\t-SIZE can be :\n\t\t-128 ;\n\t\t-512 ;\n\t-ALGORITHM can be :\n\t\t-sequential or s : it runs the sequential method ;\n\t\t-parallel or p : it runs the parallel method ;\n\t\t-task or t : it runs the parallel task method ;\n";

// ------------------------------------------------------------------------------
// -------------                         MAIN                       -------------
// ------------------------------------------------------------------------------

int main (int argc, char **argv)
{
  
	if(argc<5){
	  printf("%s", man);
	}
	
	treatment(argc,argv);
	return 0;
}