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

unsigned table; //init 0 implicite

static unsigned is_end = 1;

// callback
unsigned get (unsigned x, unsigned y)
{
  return ocean[y][x][table];
}

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

static void sand_init_homogeneous()
{
#pragma omp parallel for schedule(static) collapse(2)

  for (int y = 0; y < DIM; y++){
    for (int x = 0; x < DIM; x++){
      ocean[y][x][0] = 5;
      ocean[y][x][1] = 5;
    }
  }
}

static void sand_init_center()
{
  int center_value = 100000;
#pragma omp parallel for schedule(static) collapse(2)
  for (int y = 0; y < DIM; y++){
    for (int x = 0; x < DIM; x++){
      ocean[y][x][0] = 0;
      ocean[y][x][1] = 0;
    }
  }
  ocean[DIM/2][DIM/2][0] = center_value;
  ocean[DIM/2][DIM/2][1] = center_value;
}

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

static void copy(int table){
  #pragma omp parallel for collapse(2)
  for (int y = 0; y < DIM; y++){
    for (int x = 0; x < DIM; x++){
      ocean[y][x][1-table] = ocean[y][x][table];
    }
  }
}

void compute_cell(int x, int y, int div4){
  ocean[x][y][1-table] -= div4*4;
    ocean[x-1][y][1-table] += div4;
    ocean[x+1][y][1-table] += div4;
    ocean[x][y-1][1-table] += div4;
    ocean[x][y+1][1-table] += div4;
    is_end = 1;
}

float *compute_seq(unsigned iterations){
  if(!is_end){
    return DYNAMIC_COLORING;
  }
  is_end = 0;

  for (unsigned i = 0; i < iterations; i++){
    for (int x = 1; x < DIM-1; x=x+2){
      for (int y = 1; y < DIM-1; y++){
        if(ocean[y][x][table] >= MAX_HEIGHT)
        {
        int div4 = ocean[x][y][table]/4;
            compute_cell(x,y,div4);
        }
      } 
    }
    for (int x = 2; x < DIM-1; x=x+2){
      for (int y = DIM-2; y > 0; y--){
        if(ocean[y][x][table] >= MAX_HEIGHT)
        {
          int div4 = ocean[x][y][table]/4;
            compute_cell(x,y,div4);
        }
      }
    }
    table = 1 - table;
    if(is_end){
      copy(table);
    }
  }
  return DYNAMIC_COLORING;
}

float *compute_parallel(unsigned iterations){
  if(!is_end){
    return DYNAMIC_COLORING;
  }
  is_end = 0;


  for (unsigned i = 0; i < iterations; i++){
    #pragma omp parallel for collapse(2)
    for (int x = 1; x < DIM-1; x++){
      for (int y = 1; y < DIM-1; y++){
        //printf("%d et %d\n",omp_get_thread_num(),omp_get_num_threads());
        if(ocean[y][x][table] >= MAX_HEIGHT)
        {
            int div4 = ocean[y][x][table]/4;
            compute_cell(x, y, div4);
        }
      } 
    }
    table = 1 - table;
    if(is_end){
      copy(table);
    }
  }
  return DYNAMIC_COLORING;
}


void treatment_task(int x, int y){
  int div4 = ocean[x][y][table]/4;

  #pragma omp critical
  {
    compute_cell(x, y, div4);
  }
}

float *compute_task(unsigned iterations)
{
  if(!is_end){
    return DYNAMIC_COLORING;
  }
  is_end = 0;


  for (unsigned i = 0; i < iterations; i++){
    for (int x = 1; x < DIM-1; x =x+2){
      for (int y = 1; y < DIM-1; y++){
        #pragma omp task
        if(ocean[x][y][table] >= MAX_HEIGHT)
        {
          //depend(in:ocean[x][y-1][table], ocean[x][y+1][table]) depend(out:ocean[x][y-1][table], ocean[x][y+1][table])
          treatment_task(x,y);
        }
      }
    }

    for (int x = 2; x < DIM-1; x = x+2){
      for (int y = DIM-2; y > 0 ; y--){
        #pragma omp task
        if(ocean[x][y][table] >= MAX_HEIGHT)
        {
          //depend(in:ocean[x][y-1][table], ocean[x][y+1][table]) depend(out:ocean[x][y-1][table], ocean[x][y+1][table])
          treatment_task(x,y);
        }
      }
    }

    #pragma omp taskwait
    table = 1 - table;
    if(is_end){
      copy(table);
    }
  }
  return DYNAMIC_COLORING;
}

int seq (int argc, char **argv)
{
  display_init (argc, argv,
                DIM,                // dimension ( = x = y) du tas
                MAX_HEIGHT,         // hauteur maximale du tas
                get,                // callback func
                compute_seq);       // callback func
}

int parallel (int argc, char **argv)
{
  omp_set_nested(1);
  display_init (argc, argv,
                DIM,                // dimension ( = x = y) du tas
                MAX_HEIGHT,         // hauteur maximale du tas
                get,                // callback func
                compute_parallel);  // callback func

}

int parallel_task (int argc, char **argv)
{
  #pragma omp parallel

  #pragma omp single
  {
    display_init (argc, argv,
                  DIM,              // dimension ( = x = y) du tas
                  MAX_HEIGHT,       // hauteur maximale du tas
                  get,              // callback func
                  compute_task);    // callback func
  }
}


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

int main (int argc, char **argv)
{

	if(argc<3){
    sand_init(1);
		seq(argc,argv);
		return 0;
	}

	init_compare(argc, argv);

	return 0;
}