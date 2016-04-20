
#include "display.h"

#include <math.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>

#define DIM 128
#define MAX_HEIGHT  4

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
  for (int y = 0; y < DIM; y++){
    for (int x = 0; x < DIM; x++){
      ocean[y][x][1-table] = ocean[y][x][table];
    }
  }
}

void traitement(int x, int y){
  int mod4 = ocean[x][y][table]%4;
  int div4 = ocean[x][y][table]/4;
  #pragma omp critical
  {
    ocean[x][y][1-table] -= div4*4;
    ocean[x-1][y][1-table] += div4;
    ocean[x+1][y][1-table] += div4;
    ocean[x][y-1][1-table] += div4;
    ocean[x][y+1][1-table] += div4;
    is_end = 1;
  }
}

float *compute(unsigned iterations)
{
  if(!is_end){
    return DYNAMIC_COLORING;
  }
  is_end = 0;

  for (unsigned i = 0; i < iterations; i++){
    
    for (int x = 1; x < DIM-1; x =x+2){
      #pragma omp task
      for (int y = 1; y < DIM-1; y++){
        if(ocean[x][y][table] >= MAX_HEIGHT)
        {
          traitement(x,y);
        }
      }
    }

    for (int x = 2; x < DIM-1; x = x+2){
      #pragma omp task
      for (int y = DIM-2; y > 0 ; y--){
        if(ocean[x][y][table] >= MAX_HEIGHT)
        {
          traitement(x,y);
        }
      }
    }

    #pragma omp barrier
    table = 1 - table;
    if(is_end){
      copy(table);
    }
  }
  return DYNAMIC_COLORING;
}

int parallel_task (int argc, char **argv, int sand_init)
{
  #pragma omp parallel

  #pragma omp single
  {
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
                  get,              // callback func
                  compute);         // callback func
  }
}
