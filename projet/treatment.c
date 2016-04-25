#include "display.h"
#include "treatment.h"

#include <omp.h>
#include <stdio.h>
#include <stdlib.h>


#define DIM 128
#define MAX_HEIGHT  4

#define TIME_DIFF(t1, t2) \
        ((t2.tv_sec - t1.tv_sec) * 1000000 + (t2.tv_usec - t1.tv_usec))

unsigned ocean[DIM*DIM];
// 
static unsigned is_end = 1;

// ------------------------------------------------------------------------------
// -------------               Fonctions Utilitaires                -------------
// ------------------------------------------------------------------------------

/** 
  * Print matrices line after line on the shell.
  * Is adapted to print values until 999 without deformations.
*/
static void print()
{
    for (int x = 0; x < DIM; x++)
    {
      for (int y = 0; y < DIM; y++) 
      {
        if(ocean[x*DIM+y]< 10)
        {
          printf("  %d ",ocean[x*DIM+y]);
        }
        else if (ocean[x*DIM+y] < 100)
        {
          printf(" %d ",ocean[x*DIM+y]);
        }
        else
        {
          printf("%d ",ocean[x*DIM+y]);
      }
      printf("\n");
    }
    printf("\n");
}

/** 
  * initialize matrices with 5 grain in every case.
*/
void sand_init_homogeneous()
{
//#pragma omp parallel for schedule(static) collapse(2)
    for (int x = 0; x < DIM; x++){
      for (int y = 0; y < DIM; y++){
      ocean[x*DIM+y] = 5;
    }
  }
}

/** 
  * initialize matrices with no one grain except in the middle case which have 100000 grain.
*/
void sand_init_center()
{
  int center_value = 100000;
//#pragma omp parallel for schedule(static) collapse(2)
  for (int x = 0; x < DIM; x++){
    for (int y = 0; y < DIM; y++){
      ocean[x*DIM+y] = 0;
    }
  }
  ocean[DIM*DIM/2+DIM/2] = center_value;
}

// ------------------------------------------------------------------------------
// -------------              Fonctions de traitement               -------------
// ------------------------------------------------------------------------------

    // ------------------------------------------------------------------------------
    // -------------                     checkboat                      -------------
    // ------------------------------------------------------------------------------


// callback
unsigned get (unsigned x, unsigned y)
{
  return ocean[x*DIM+y];
}

/** 
  * Applied modifications into 1-table for ocean[y][x] sandpile.
*/
void compute_cell(int x, int y, int div4){
  ocean[x*DIM+y]      -= div4*4;
  ocean[x*DIM+y+1]    += div4;
  ocean[x*DIM+y-1]    += div4;
  ocean[(x-1)*DIM+y]  += div4;
  ocean[(x+1)*DIM+y]  += div4;
  is_end = 1;
}

/** 
  * Compute fonction for sequentiel treatment
*/
static float *compute_seq(unsigned iterations){
  if(!is_end){
    return DYNAMIC_COLORING;
  }
  is_end = 0;

 
  for (unsigned i = 0; i < iterations; i++)
    {
        for (int x = 1; x < DIM-1; x++)
      {
         for (int y = 1; y < DIM-1; y++)
         {
            if(ocean[x*DIM+y] >= MAX_HEIGHT)
            {
              int div4 = ocean[x*DIM+y]/4;
              compute_cell(x,y,div4);
            }
         }
    	}
    }
  }
  return DYNAMIC_COLORING;
}

/** 
  * Compute fonction for sequentiel treatment which alternate at every line
*/
static float *compute_seq_alternative(unsigned iterations){
  if(!is_end){
    return DYNAMIC_COLORING;
  }
  is_end = 0;

  for (unsigned i = 0; i < iterations; i++){
    for (int x = 1; x < DIM-1; x=x+2){
      for (int y = 1; y < DIM-1; y++){
        if(ocean[x*DIM+y] >= MAX_HEIGHT)
        {
          int div4 = ocean[x*DIM+y]/4;
          compute_cell(x,y,div4);
        }
      } 
    }
    for (int x = 2; x < DIM-1; x=x+2){
      for (int y = DIM-2; y > 0; y--){
        if(ocean[x*DIM+y] >= MAX_HEIGHT)
        {
          int div4 = ocean[x*DIM+y]/4;
          compute_cell(x,y,div4);
        }
      }
    }
  }
  return DYNAMIC_COLORING;
}

/** 
  * Compute fonction for collapsed treatment
*/
float *compute_parallel(unsigned iterations){
  if(!is_end){
    return DYNAMIC_COLORING;
  }
  is_end = 0;

  for (unsigned i = 0; i < iterations; i++){
    #pragma omp parallel for collapse(2)
    for (int x = 1; x < DIM-1; x=x+2){
      for (int y = 1; y < DIM-1; y++){
        if(ocean[x*DIM+y] >= MAX_HEIGHT)
        {
          int div4 = ocean[x*DIM+y]/4;
            #pragma omp critical
          compute_cell(x,y,div4);
        }
      } 
    }
    #pragma omp parallel for collapse(2)
    for (int x = 2; x < DIM-1; x=x+2){
      for (int y = DIM-2; y > 0; y--){
        if(ocean[x*DIM+y] >= MAX_HEIGHT)
        {
          int div4 = ocean[x*DIM+y]/4;
            #pragma omp critical
          compute_cell(x,y,div4);
        }
      }
    }
    table = 1 - table;
  }
  return DYNAMIC_COLORING;
}

/** 
  * Adapted treatment to one case, with task specificity
*/
static void treatment_task(int x, int y){
  int div4 = ocean[x*DIM+y]/4;

  //#pragma omp critical
  {
    compute_cell(x,y,div4);
  }
}

/** 
  * Compute fonction for task treatment
*/
static float *compute_task(unsigned iterations)
{
  if(!is_end){
    return DYNAMIC_COLORING;
  }
  is_end = 0;

    // 2 passage (pour optimiser le cache ?)
  for (unsigned i = 0; i < iterations; i++){
    for (int x = 1; x < DIM-1; x =x+2){
      for (int y = 1; y < DIM-1; y++){
        #pragma omp task depend(in:ocean[x*DIM+y], ocean[(x-1)*DIM+y], ocean[(x+1)*DIM+y], ocean[x*DIM+y-1], ocean[x*DIM+y+1]) depend(out:ocean[x*DIM+y], ocean[(x-1)*DIM+y], ocean[(x+1)*DIM+y], ocean[x*DIM+y-1], ocean[x*DIM+y+1])
        if(ocean[x*DIM+y] >= MAX_HEIGHT)
        {
          treatment_task(x,y);
        }
      }
    }

    for (int x = 2; x < DIM-1; x = x+2){
      for (int y = DIM-2; y > 0 ; y--){
        #pragma omp task depend(in:ocean[x*DIM+y], ocean[(x-1)*DIM+y], ocean[(x+1)*DIM+y], ocean[x*DIM+y-1], ocean[x*DIM+y+1]) depend(out:ocean[x*DIM+y], ocean[(x-1)*DIM+y], ocean[(x+1)*DIM+y], ocean[x*DIM+y-1], ocean[x*DIM+y+1])
        if(ocean[x*DIM+y] >= MAX_HEIGHT)
        {
          treatment_task(x,y);
        }
      }
    }

    //#pragma omp taskwait
    table = 1 - table;
  }
  return DYNAMIC_COLORING;
}

    // ------------------------------------------------------------------------------
    // -------------                     doubleline                     -------------
    // ------------------------------------------------------------------------------


/** 
  * Compute fonction for sequentiel treatment which alternate at every line
*/
static float *compute_seq_alternative(unsigned iterations){
  if(!is_end){
    return DYNAMIC_COLORING;
  }
  is_end = 0;

  for (unsigned i = 0; i < iterations; i++){
    for (int x = 1; x < DIM-1; x=x+2){
      for (int y = 1; y < DIM-1; y++){
        if(ocean[x*DIM+y] >= MAX_HEIGHT)
        {
          int div4 = ocean[x*DIM+y]/4;
          compute_cell(x,y,div4);
        }
      } 
    }
    for (int x = 2; x < DIM-1; x=x+2){
      for (int y = DIM-2; y > 0; y--){
        if(ocean[x*DIM+y] >= MAX_HEIGHT)
        {
          int div4 = ocean[x*DIM+y]/4;
          compute_cell(x,y,div4);
        }
      }
    }
  }
  return DYNAMIC_COLORING;
}

/** 
  * Compute fonction for collapsed treatment
*/
float *compute_parallel(unsigned iterations){
  if(!is_end){
    return DYNAMIC_COLORING;
  }
  is_end = 0;

  for (unsigned i = 0; i < iterations; i++){
    #pragma omp parallel for collapse(2)
    for (int x = 1; x < DIM-1; x=x+2){
      for (int y = 1; y < DIM-1; y++){
        if(ocean[x*DIM+y] >= MAX_HEIGHT)
        {
          int div4 = ocean[x*DIM+y]/4;
            #pragma omp critical
          compute_cell(x,y,div4);
        }
      } 
    }
    #pragma omp parallel for collapse(2)
    for (int x = 2; x < DIM-1; x=x+2){
      for (int y = DIM-2; y > 0; y--){
        if(ocean[x*DIM+y] >= MAX_HEIGHT)
        {
          int div4 = ocean[x*DIM+y]/4;
            #pragma omp critical
          compute_cell(x,y,div4);
        }
      }
    }
    table = 1 - table;
  }
  return DYNAMIC_COLORING;
}

/** 
  * Adapted treatment to one case, with task specificity
*/
static void treatment_task(int x, int y){
  int div4 = ocean[x*DIM+y]/4;

  //#pragma omp critical
  {
    compute_cell(x,y,div4);
  }
}

/** 
  * Compute fonction for task treatment
*/
static float *compute_task(unsigned iterations)
{
  if(!is_end){
    return DYNAMIC_COLORING;
  }
  is_end = 0;

    // 2 passage (pour optimiser le cache ?)
  for (unsigned i = 0; i < iterations; i++){
    for (int x = 1; x < DIM-1; x =x+2){
      for (int y = 1; y < DIM-1; y++){
        #pragma omp task depend(in:ocean[x*DIM+y], ocean[(x-1)*DIM+y], ocean[(x+1)*DIM+y], ocean[x*DIM+y-1], ocean[x*DIM+y+1]) depend(out:ocean[x*DIM+y], ocean[(x-1)*DIM+y], ocean[(x+1)*DIM+y], ocean[x*DIM+y-1], ocean[x*DIM+y+1])
        if(ocean[x*DIM+y] >= MAX_HEIGHT)
        {
          treatment_task(x,y);
        }
      }
    }

    for (int x = 2; x < DIM-1; x = x+2){
      for (int y = DIM-2; y > 0 ; y--){
        #pragma omp task depend(in:ocean[x*DIM+y], ocean[(x-1)*DIM+y], ocean[(x+1)*DIM+y], ocean[x*DIM+y-1], ocean[x*DIM+y+1]) depend(out:ocean[x*DIM+y], ocean[(x-1)*DIM+y], ocean[(x+1)*DIM+y], ocean[x*DIM+y-1], ocean[x*DIM+y+1])
        if(ocean[x*DIM+y] >= MAX_HEIGHT)
        {
          treatment_task(x,y);
        }
      }
    }

    //#pragma omp taskwait
    table = 1 - table;
  }
  return DYNAMIC_COLORING;
}


// ------------------------------------------------------------------------------
// ------------- Fonction initiale dépendante de la méthode appelée -------------
// ------------------------------------------------------------------------------

int cb_seq (int argc, char **argv)
{
  display_init (argc, argv,
                DIM,                // dimension ( = x = y) du tas
                MAX_HEIGHT,         // hauteur maximale du tas
                cb_get,                // callback func
                cb_compute_seq);       // callback func
}

int cb_seq_alternance(int argc, char **argv)
{
  display_init (argc, argv,
              DIM,                // dimension ( = x = y) du tas
              MAX_HEIGHT,         // hauteur maximale du tas
              cb_get,                // callback func
              cb_compute_seq_alternance);       // callback func
}

int cb_parallel (int argc, char **argv)
{
  //omp_set_nested(1);
  display_init (argc, argv,
                DIM,                // dimension ( = x = y) du tas
                MAX_HEIGHT,         // hauteur maximale du tas
                cb_get,                // callback func
                cb_compute_parallel);  // callback func

}

int cb_parallel_task (int argc, char **argv)
{
  #pragma omp parallel

  #pragma omp single
  {
    display_init (argc, argv,
                  DIM,              // dimension ( = x = y) du tas
                  MAX_HEIGHT,       // hauteur maximale du tas
                  cb_get,              // callback func
                  cb_compute_task);    // callback func
  }
}