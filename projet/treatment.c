#include "display.h"
#include "treatment.h"

#include <omp.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>

#define MAX_HEIGHT  4

#define TIME_DIFF(t1, t2) \
        ((t2.tv_sec - t1.tv_sec) * 1000000 + (t2.tv_usec - t1.tv_usec))

unsigned DIM;
unsigned *ocean;

static bool is_end = false;

// ------------------------------------------------------------------------------
// -------------               Utility Functions                    -------------
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
    }
    printf("\n");
  }
}

/** 
  * initialize matrices with 5 squares in every square.
*/
void sand_init_homogeneous()
{
  for (int x = 0; x < DIM; x++)
  {
    for (int y = 0; y < DIM; y++)
    {
      if(x > 0 && x < DIM-1 && y > 0 && y < DIM-1)
        ocean[x*DIM+y] = 5;
      else
        
        ocean[x*DIM+y] = 0;
    }
  }
}

/** 
  * initialize matrices with no one grain except in the middle case which have 100000 grain.
*/
void sand_init_center()
{
  int center_value = 100000;
  for (int x = 0; x < DIM; x++)
  {
    for (int y = 0; y < DIM; y++)
    {
      ocean[x*DIM+y] = 0;
    }
  }
  ocean[DIM*DIM/2+DIM/2] = center_value;
}

// callback
unsigned get (unsigned x, unsigned y)
{
  return ocean[x*DIM+y];
}

// ------------------------------------------------------------------------------
// -------------              Treatments Functions               ----------------
// ------------------------------------------------------------------------------

/** 
  * Applied modifications into 1-table for ocean[x][y] sandpile.
*/
static inline void compute_cell(int x, int y, int div4)
{
// is the pragma atoms used when there is no pragma parallel before ?
#pragma omp atoms
  {
    ocean[x*DIM+y]      -= div4*4;
    ocean[x*DIM+y+1]    += div4;
    ocean[x*DIM+y-1]    += div4;
    ocean[(x-1)*DIM+y]  += div4;
    ocean[(x+1)*DIM+y]  += div4;
  }
  is_end = false;
}



/** 
  * Compute fonction for sequential treatment
*/
static inline float *compute_seq(unsigned iterations)
{
  if(is_end == true)
  {
    return DYNAMIC_COLORING;
  }
  is_end = true;
 
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
  return DYNAMIC_COLORING;
}

/** 
  * Compute fonction for sequential treatment which alternate at every line | Alternative Method to use the cache differently
*/
static inline float *compute_seq_alternative(unsigned iterations)
{
  if(is_end == true)
  {
    return DYNAMIC_COLORING;
  }
  is_end = true;

  for (unsigned i = 0; i < iterations; i++)
  {
    for (int x = 1; x < DIM-1; x=x+2)
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
    for (int x = 2; x < DIM-1; x=x+2)
    {
      for (int y = DIM-2; y > 0; y--)
      {
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

/*
TODO : ne peut pas fonctionner tel quel.
*/
bool init = false;
unsigned * tmp;
static inline float *compute_seq_see_stabilize(unsigned iterations)
{
  if(init == false){
    tmp = calloc(DIM,sizeof(unsigned));
    init = true;
  }

  if(is_end == true)
  {
    return DYNAMIC_COLORING;
    free(tmp);
  }
  is_end = true;
 
  for (unsigned i = 0; i < iterations; i++)
  {
    for (int x = 1; x < DIM-1; x++)
    {
      if(tmp[i] == true)
      {
        for (int y = 1; y < DIM-1; y++)
        {
          if(ocean[x*DIM+y] >= MAX_HEIGHT)
          {
            int div4 = ocean[x*DIM+y]/4;
            compute_cell(x,y,div4);
            if(ocean[(x+1)*DIM+y] >= MAX_HEIGHT)
              tmp[i+1] = true;
            if(ocean[(x-1)*DIM+y] >= MAX_HEIGHT)
              tmp[i-1] = true;
          }
        }
      }

    }
  }
  return DYNAMIC_COLORING;
}


/* methode nécessitant deux matrices : une écriture une lecture
int bol = 0;
unsigned * tmp;
static inline float *compute_seq(unsigned iterations)
{
  if(bol == 0){
    tmp = malloc(sizeof(unsigned)*DIM*DIM);
    bol = 1;
  }

  if(is_end != 1)
  {
    print();
    return DYNAMIC_COLORING;
  }
  is_end = 0;
 
  for (unsigned i = 0; i < iterations; i++)
  {
    for (int x = 1; x < DIM-1; x++)
    {
      for (int y = 1; y < DIM-1; y++)
      {
        int val = ocean[x*DIM+y]%4;
        //if(x > 0)
          val += ocean[(x-1)*DIM+y]/4;
        //if(x < DIM-1)
          val += ocean[(x+1)*DIM+y]/4;
        //if(y > 0)
          val += ocean[x*DIM+y-1]/4;
        //if(y < DIM-1)
          val += ocean[x*DIM+y+1]/4;
        if(val != ocean[x*DIM+y]){
          is_end = 1;
        }
        tmp[x*DIM+y] = val;
      }
    }

    //copy nécessaire pour fonctionner. aura nécessairement un temps horrible.

  }
  return DYNAMIC_COLORING;
}
*/

static inline float *compute_seq_doubleline(unsigned iterations)
{
  if(is_end == true)
  {
    return DYNAMIC_COLORING;
  }
  is_end = true;

  for (unsigned i = 0; i < iterations; i++)
  {
    for (int x = 1; x < DIM-1; x++)
    {
      for (int y = 1; y < DIM-5; y+=4)
      {
        if(ocean[x*DIM+y] >= MAX_HEIGHT)
        {
          int div4 = ocean[x*DIM+y]/4;
          compute_cell(x,y,div4);
        }
        if(ocean[x*DIM+y+1] >= MAX_HEIGHT)
        {
          int div4 = ocean[x*DIM+y+1]/4;
          compute_cell(x,y+1,div4);
        }
        if(ocean[x*DIM+y+2] >= MAX_HEIGHT)
        {
          int div4 = ocean[x*DIM+y+2]/4;
          compute_cell(x,y+2,div4);
        }
        if(ocean[x*DIM+y+3] >= MAX_HEIGHT)
        {
          int div4 = ocean[x*DIM+y+3]/4;
          compute_cell(x,y+3,div4);
        }
      } 
    }
  }
  return DYNAMIC_COLORING;
}

/** 
  * Compute fonction for collapsed treatment 
*/
static inline float *compute_parallel_for(unsigned iterations)
{
  if(is_end == true)
  {
    return DYNAMIC_COLORING;
  }
  is_end = true;

 
  for (unsigned i = 0; i < iterations; i++)
  {
#pragma omp parallel for //collapse(2)
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
  return DYNAMIC_COLORING;
}

/** 
  * Compute fonction for collapsed alternative treatment | Alternative Method to use the cache differently
*/
static inline float *compute_parallel_for_alternative(unsigned iterations)
{

  if(is_end == true)
  {
    return DYNAMIC_COLORING;
  }
  is_end = true;

  for (unsigned i = 0; i < iterations; i++)
  {
#pragma omp parallel for
    for (int x = 1; x < DIM-1; x=x+2)
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
#pragma omp parallel for
    for (int x = 2; x < DIM-1; x=x+2)
    {
      for (int y = DIM-2; y > 0; y--)
      {
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
//TODO
float *compute_parallel_multiple_lines(unsigned iterations)
{
  if(is_end == true)
  {
    return DYNAMIC_COLORING;
  }
  is_end = true;

  omp_set_num_threads(17);
  int nb_thread = omp_get_num_threads();
  int nb_lines = DIM / nb_thread;

  for (unsigned i = 0; i < iterations; i++)
  {
#pragma omp parallel
    {
      int thread_num = omp_get_thread_num();
      int my_first_line = nb_lines*thread_num;
      int my_last_line = my_first_line + nb_lines*DIM;
      for (int x = my_first_line+1; ((x < my_last_line)&&(x < DIM-1)); x++)
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

/*
      for(int j = my_first_line ; j < my_first_line+nb_lines*(DIM-1); j++)
      {
        if(ocean[j] >= MAX_HEIGHT)
        {
          compute_cell(j / DIM, J % DIM, ocean[j]/4);
        }
      }
*/

    }
  }
  return DYNAMIC_COLORING;
}

/** 
  * Adapted treatment to one case, with tasks specificities
*/
static inline void treatment_task(int x, int y)
{
  int div4 = ocean[x*DIM+y]/4;
  compute_cell(x,y,div4);
}

/** 
  * Compute fonction for task treatment
*/
static inline float *compute_parallel_task(unsigned iterations)
{
  if(is_end == true)
  {
    return DYNAMIC_COLORING;
  }
  is_end = true;

  for (unsigned i = 0; i < iterations; i++)
  {
    for (int x = 1; x < DIM-1; x++)
    {
      for (int y = 1; y < DIM-1; y++)
      {
#pragma omp task depend(in:ocean[x*DIM+y], ocean[(x-1)*DIM+y], ocean[(x+1)*DIM+y], ocean[x*DIM+y-1], ocean[x*DIM+y+1]) depend(out:ocean[x*DIM+y], ocean[(x-1)*DIM+y], ocean[(x+1)*DIM+y], ocean[x*DIM+y-1], ocean[x*DIM+y+1])
        if(ocean[x*DIM+y] >= MAX_HEIGHT)
        {
          treatment_task(x,y);
        }
      }
    }
  }
  return DYNAMIC_COLORING;
}

/** 
  * Compute fonction for task treatment | Alternative Method to move the cache differently
*/
static inline float *compute_parallel_task_alternative(unsigned iterations)
{
  if(is_end == true)
  {
    return DYNAMIC_COLORING;
  }
  is_end = true;

  for (unsigned i = 0; i < iterations; i++)
  {
    for (int x = 1; x < DIM-1; x =x+2)
    {
      for (int y = 1; y < DIM-1; y++)
      {
#pragma omp task depend(in:ocean[x*DIM+y], ocean[(x-1)*DIM+y], ocean[(x+1)*DIM+y], ocean[x*DIM+y-1], ocean[x*DIM+y+1]) depend(out:ocean[x*DIM+y], ocean[(x-1)*DIM+y], ocean[(x+1)*DIM+y], ocean[x*DIM+y-1], ocean[x*DIM+y+1])
        if(ocean[x*DIM+y] >= MAX_HEIGHT)
        {
          treatment_task(x,y);
        }
      }
    }

    for (int x = 2; x < DIM-1; x = x+2)
    {
      for (int y = DIM-2; y > 0 ; y--)
      {
#pragma omp task depend(in:ocean[x*DIM+y], ocean[(x-1)*DIM+y], ocean[(x+1)*DIM+y], ocean[x*DIM+y-1], ocean[x*DIM+y+1]) depend(out:ocean[x*DIM+y], ocean[(x-1)*DIM+y], ocean[(x+1)*DIM+y], ocean[x*DIM+y-1], ocean[x*DIM+y+1])
        if(ocean[x*DIM+y] >= MAX_HEIGHT)
        {
          treatment_task(x,y);
        }
      }
    }
  }
  return DYNAMIC_COLORING;
}

// --------------------------------------------------------------------------------------------------------
// ----------------------------- Initials functions to adapt at arguments  --------------------------------
// --------------------------------------------------------------------------------------------------------

int display(int argc, char ** argv)
{
  switch((int)argv[4][0])
  {
    case 115 : // ascii s
      display_init (argc, argv,
                DIM,                // dimension ( = x = y) du tas
                MAX_HEIGHT,         // hauteur maximale du tas
                get,                // callback func
                compute_seq);       // callback func
      break;
    case 97 : //ascii a
      display_init (argc, argv,
              DIM,                // dimension ( = x = y) du tas
              MAX_HEIGHT,         // hauteur maximale du tas
              get,                // callback func
              compute_seq_alternative); // callback func
      break;
    case 102 : //ascii f
      omp_set_nested(1);
      display_init (argc, argv,
                DIM,                // dimension ( = x = y) du tas
                MAX_HEIGHT,         // hauteur maximale du tas
                get,                // callback func
                compute_parallel_for);  // callback func
      break;
    case 116 : //ascii t
#pragma omp parallel

#pragma omp single
      {
        display_init (argc, argv,
                      DIM,              // dimension ( = x = y) du tas
                      MAX_HEIGHT,       // hauteur maximale du tas
                      get,              // callback func
                      compute_parallel_task);    // callback func
      }
      break;
    case 109 : //ascii m
      display_init (argc, argv,
                DIM,                // dimension ( = x = y) du tas
                MAX_HEIGHT,         // hauteur maximale du tas
                get,                // callback func
                compute_parallel_multiple_lines);  // callback func
      break;
    default :
      printf("Unrecognize Algorithm. Please, read our manual by using ./sand\n");
      break;
  }
}

void without_display(float * (*compute_func_t) (unsigned iterations))
{
  // TODO : prendre une décision pour la valeur de interation par défault.
  do{
    compute_func_t(10);
  }
  while(is_end != 0);
}

int performance(int argc, char ** argv)
{
  unsigned long temps;
  struct timeval t1, t2;

  gettimeofday(&t1,NULL);

  switch((int)argv[4][0])
  {
    case 115 : // ascii s
      printf("Sequential Algorithm\n");
      without_display(compute_seq);
      break;
    case 83 : // ascii S
      without_display(compute_seq_alternative);
      printf("Sequential Algorithm | Double Loop\n");
      break;
    case 108 : //ascii l
      without_display(compute_seq_doubleline);
      printf("Sequential Algorithm | 4 lines each loop turn\n");
      break;
    case 102 : //ascii f
      without_display(compute_parallel_for);
      printf("Parallel Algorithm\n");
      break;
    case 70 : //ascii F
      without_display(compute_parallel_for_alternative);
      printf("Parallel Algorithm | Double Loop\n");
      break;
    case 116 : //ascii t
      without_display(compute_parallel_task);
      printf("Parallel Task Algorithm\n");
      break;
    case 84 : //ascii T
      without_display(compute_parallel_task_alternative);
      printf("Parallel Task Algorithm | Double Loop\n");
      break;
    case 109 : //ascii m
      without_display(compute_parallel_multiple_lines);
      printf("Parallel Algorithm | Multiple Lines\n");
      break;
    default :
      printf("Unrecognize Algorithm. Please, read our manual by using ./sand\n");
      break;
  }
  gettimeofday(&t2,NULL);
   
  temps = TIME_DIFF(t1,t2);
  printf("Algorithm time = %ld.%03ldms \n", temps/1000, temps%1000);
  print();
  free(ocean);
}


void treatment(int argc, char ** argv)
{
  DIM = strtol(argv[3],NULL,10);
  if(!DIM)
  {
    printf("Unrecognize size. Please, read our manual by using ./sand\n");
    return;
  }
  
  // pas d'autres allocations durant le traitement,
  // ne sera pas free par nous à la fin si on entre dans display car aurait nécessité de modifier display 
  // MAIS sera nettoyé sans perte à la fermture du processus en tout les cas.
  ocean = malloc(sizeof(unsigned)*DIM*DIM); 

   //initialisation de la matrice
  switch((int)argv[2][0])
  {
    case 99 : //ascii c
      sand_init_center();
      printf("Centered Case\n");
      break;
    case 104: //ascii h
      sand_init_homogeneous();
      printf("Homogeneous Case\n");
      break;
    default :
      printf("Unrecognize configuration. Please, read our manual by using ./sand\n");
      break;
  }

  // objectif
  if(!atoi(argv[1]))
  {
    performance(argc,argv);
  }
  else
  {
    display(argc,argv);  
  }
}
