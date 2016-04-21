
#include "display.h"
#include "parallel.h"

#include <math.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>


#define DIM 128
#define MAX_HEIGHT  4

#define TIME_DIFF(t1, t2) \
        ((t2.tv_sec - t1.tv_sec) * 1000000 + (t2.tv_usec - t1.tv_usec))

unsigned ocean[DIM][DIM][2];

unsigned table; //init 0 implicite

static unsigned is_end = 1;

// callback
unsigned get_parallel (unsigned x, unsigned y)
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



static void copy(int table){
  #pragma omp parallel for collapse(2)
  for (int y = 0; y < DIM; y++){
    for (int x = 0; x < DIM; x++){
      ocean[y][x][1-table] = ocean[y][x][table];
    }
  }
}

static float *compute_parallel(unsigned iterations){
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
            ocean[y][x][1-table] -= div4*4;
            ocean[y-1][x][1-table] += div4;
            ocean[y+1][x][1-table] += div4;
            ocean[y][x-1][1-table] += div4;
            ocean[y][x+1][1-table] += div4;
            is_end = 1;
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
int parallel (int argc, char **argv, int sand_init)
{
  DIM = size ;
  omp_set_nested(1);

  if(!sand_init)
  {
    sand_init_center();
  }else
  {
    sand_init_homogeneous();
  }

  display_init (argc, argv,
                DIM,              // dimension ( = x = y) du tas
                MAX_HEIGHT,       // hauteur maximale du tas
                get_parallel,              // callback func
                compute_parallel);         // callback func

}
