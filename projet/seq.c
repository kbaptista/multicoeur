
#include "display.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>

#define DIM 17
#define MAX_HEIGHT  4
#define TIME_DIFF(t1, t2) \
        ((t2.tv_sec - t1.tv_sec) * 1000000 + (t2.tv_usec - t1.tv_usec))

unsigned ocean[DIM][DIM][2];

unsigned table; //init 0 implicite

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
  int center_value = 10000 ;
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
float *compute(unsigned iterations){
  for (unsigned i = 0; i < iterations; i++){
    for (int y = 1; y < DIM-1; y++){
      for (int x = 1; x < DIM-1; x++){
        if(ocean[y][x][table] >= MAX_HEIGHT)
        {
            int mod4 = ocean[y][x][table]%4;
            int div4 = ocean[y][x][table]/4;
            ocean[y][x][1-table] -= div4*4;
            ocean[y-1][x][1-table] += div4;
            ocean[y+1][x][1-table] += div4;
            ocean[y][x-1][1-table] += div4;
            ocean[y][x+1][1-table] += div4;
        }
      } 
    }
    table = 1 - table;
    copy(table);
  }
  return DYNAMIC_COLORING;
}

int main (int argc, char **argv)
{
  unsigned long temps;
  struct timeval t1, t2;

  gettimeofday(&t1,NULL);

  //sand_init_homogeneous();
  sand_init_center();

  display_init (argc, argv,
                DIM,              // dimension ( = x = y) du tas
                MAX_HEIGHT,       // hauteur maximale du tas
                get,              // callback func
                compute);         // callback func
  
  gettimeofday(&t2,NULL);
   
  temps = TIME_DIFF(t1,t2);

  printf("time = %ld.%03ldms \n", temps/1000, temps%1000);
  
}
