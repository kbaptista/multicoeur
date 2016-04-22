#include "display.h"

#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define DIM 128
#define MAX_HEIGHT  4

#define TIME_DIFF(t1, t2) \
        ((t2.tv_sec - t1.tv_sec) * 1000000 + (t2.tv_usec - t1.tv_usec))

unsigned ocean[DIM][DIM][2];

// indice indiquant laquelle des deux tables est actuellement en lecture.
unsigned table; //init 0 implicite

//booléen qui change après le dernier traitement, indique donc si compute a encore du traitement à faire ou non.
static unsigned is_end = 1;

// ------------------------------------------------------------------------------
// -------------               Fonctions Utilitaires                -------------
// ------------------------------------------------------------------------------

/** 
  * Print matrices time after time on the shell.
  * Is adapted to print values until 999 without deformations.
*/
static void print(int table)
{
    for (int y = 0; y < DIM; y++)
    {
      for (int x = 0; x < DIM; x++) 
      {
        if(ocean[y][x][table] < 10)
        {
          printf("  %d ",ocean[y][x][table]);
        }
        else if (ocean[y][x][table] < 100)
        {
          printf(" %d ",ocean[y][x][table]);
        }
        else
        {
          printf("%d ",ocean[y][x][table]);
        }
      }
      printf("\n");
    }
    printf("\n");
}

/** 
  * initialise matrices with 5 grain in every case.
*/
static void sand_init_homogeneous()
{
//#pragma omp parallel for schedule(static) collapse(2)

  for (int y = 0; y < DIM; y++){
    for (int x = 0; x < DIM; x++){
      ocean[y][x][0] = 5;
      ocean[y][x][1] = 5;
    }
  }
}

/** 
  * initialise matrices with no one grain except in the middle case which have 100000 grain.
*/
static void sand_init_center()
{
  int center_value = 100000;
//#pragma omp parallel for schedule(static) collapse(2)
  for (int y = 0; y < DIM; y++){
    for (int x = 0; x < DIM; x++){
      ocean[y][x][0] = 0;
      ocean[y][x][1] = 0;
    }
  }
  ocean[DIM/2][DIM/2][0] = center_value;
  ocean[DIM/2][DIM/2][1] = center_value;
}

/** 
  * initialise matrices : 
  * 0 means centering configuration
  * 1 means homogeneous configuration
*/
void sand_init (int init)
{
  if(!init)
  {
    sand_init_center();
  }
  else
  {
    sand_init_homogeneous();
  }
}



// ------------------------------------------------------------------------------
// -------------        Fonctions de traitement des Arguments       -------------
// ------------------------------------------------------------------------------

char * man = "usage : ./sand <INITIALIZATION> <SIZE> <ALGORITHM> \n\n\t-INITIALIZATION can be :\n\t\t-homogeneous or h : it starts the homogeneous case ;\n\t\t-centered or c : it starts the centered case ;\n\t-SIZE can be :\n\t\t-128 ;\n\t\t-512 ;\n\t-ALGORITHM can be :\n\t\t-sequential or s : it runs the sequential method ;\n\t\t-parallel or p : it runs the parallel method ;\n\t\t-task or t : it runs the parallel task method ;\n";

/* Add your treatment choices here */
void treatment_compare(int argc, char **argv, int init_compare_res, int size_compare_res){

  if(!strcmp(argv[3], "sequential") || !strcmp(argv[3], "s")){
    	printf("Treatment Sequential\n\n");
      sand_init(init_compare_res);
      seq(argc, argv);
    }
    else if(!strcmp(argv[3], "for") || !strcmp(argv[3], "f")){
    	printf("Treatment Parallel For\n\n");
      sand_init(init_compare_res);
      parallel(argc, argv);
    }
    else if(!strcmp(argv[3], "task") || !strcmp(argv[3], "t")){
    	printf("Treatment Parallel Task\n\n");
      sand_init(init_compare_res);
      parallel_task(argc, argv);
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

// ------------------------------------------------------------------------------
// -------------                         MAIN                       -------------
// ------------------------------------------------------------------------------

int main (int argc, char **argv)
{

	if(argc<2){
    sand_init(1);
		seq(argc,argv);
		return 0;
	}

	init_compare(argc, argv);

	return 0;
}