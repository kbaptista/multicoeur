#include "display.h"
#include "treatment.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>


static char * man = {"usage : ./sand <GUI> <INITIALIZATION> <SIZE> <ALGORITHM>\n\n"
    "\t-GUI can be :\n"
    "\t\t- 0 to disable GUI\n"
    "\t\t- 1 to enable GUI\n"
    "\t-INITIALIZATION can be :\n"
    "\t\t- h : it starts the homogeneous case ;\n"
    "\t\t- c : it starts the centered case ;\n"
    "\t-SIZE can be :\n"
    "\t\t- 128 ;\n"
    "\t\t- 512 ;\n"
    "\t-ALGORITHM can be :\n"
    "\t\t- s : it runs the sequential expander method ;\n"
    "\t\t- S : it runs the sequential gatherer method ;\n"
    "\t\t- u : it runs the sequential unwrapped expander method ;\n"
    "\t\t- U : it runs the sequential unwrapped gatherer method ;\n"
    "\t\t- P : it runs the parallel gatherer method ;\n"
    "\t\t- I : it runs the parallel gatherer method which synchronise each p iterations ;\n"};

// ------------------------------------------------------------------------------
// -------------                         MAIN                       -------------
// ------------------------------------------------------------------------------

int main (int argc, char **argv)
{
  
	if(argc<5){
	  printf("%s", man);
	  return -1;
	}
	
	treatment(argc,argv);
	return 0;
}