#include "display.h"
#include "treatment.h"

#include <math.h>
#include <omp.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>

#define MAX_HEIGHT  4

#define TIME_DIFF(t1, t2) \
        ((t2.tv_sec - t1.tv_sec) * 1000000 + (t2.tv_usec - t1.tv_usec))

unsigned DIM;
unsigned *ocean;

static bool is_end = false;

// The comments may refer to ocean[x][y] in place of the real value which is ocean[x*DIM+y] for lisibility concerns.

// ------------------------------------------------------------------------------
// -------------               Utility Functions                    -------------
// ------------------------------------------------------------------------------

/** 
  * Print matrices line after line on the shell.
  * Rows and Columns are printed too.
*/
static void print()
{
  for (int x = 0; x < DIM; x++)
  {
    //print the line number, or the ordinate value
    printf("[%3d]",x );
    for (int y = 0; y < DIM; y++) 
    {
      printf(" %3d ",ocean[x*DIM+y]);
    }
    printf("\n");
  }
  printf("     ");

  //when printing ocean is done, we print the row number, or the abcissa value
  for (int i = 0; i < DIM; ++i)
  {
    printf("[%3d]",i );
  }
  printf("\n");
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
      //if the cell is intern
      if(x > 0 && x < DIM-1 && y > 0 && y < DIM-1)
        ocean[x*DIM+y] = 5;
      //if the cell is in the border, it simulates the void, the table's edges.
      else
        ocean[x*DIM+y] = 0;
    }
  }
}

/** 
  * initialize an empty matrice with 100 000 grains in its center.
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
  * Divide the [x][y] cell content and repart it in the neighbours.
  * div4 HAVE TO BE the results of ocean[x][y]/4
*/
static inline void compute_cell_expander(int x, int y, int div4)
{
  ocean[x*DIM+y]      -= div4*4;
  ocean[x*DIM+y+1]    += div4;
  ocean[x*DIM+y-1]    += div4;
  ocean[(x-1)*DIM+y]  += div4;
  ocean[(x+1)*DIM+y]  += div4;
  is_end = false;
}

/** 
  * Divide the [x][y] cell content and repart it in the neighbours.
*/
static inline unsigned compute_cell_gatherer(int x, int y, int table)
{
  unsigned val = ocean[table+x*DIM+y]%4 + ocean[table+(x-1)*DIM+y]/4 + ocean[table+(x+1)*DIM+y]/4 + ocean[table+x*DIM+y-1]/4 + ocean[table+x*DIM+y+1]/4;
  if (val >= MAX_HEIGHT ) is_end = false ;

  //printf("\n1:[%d] %d = %d + %d + %d + %d + %d \n",x,val,ocean[table+x*DIM+y]%4,ocean[table+(x-1)*DIM+y]/4,ocean[table+(x+1)*DIM+y]/4,ocean[table+x*DIM+y-1]/4,ocean[table+x*DIM+y+1]/4);
  return val;
}

/** 
  * Compute fonction for sequential treatment
*/
static inline float *compute_seq_expander(unsigned iterations)
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
          compute_cell_expander(x,y,div4);
        }
      }
    }
  }
  return DYNAMIC_COLORING;
}

/** 
  * Compute fonction for sequential treatment
*/
static inline float *compute_seq_gatherer(unsigned iterations)
{
  if(is_end == true)
  {
    return DYNAMIC_COLORING;
  }
  is_end = true;
 
  int table = 0;
  int table_2 = DIM*DIM;
 
  for (unsigned i = 0; i < iterations; i++)
  {
    for (int x = 1; x < DIM-1; x++)
    {
      for (int y = 1; y < DIM-1; y++)
      {
        ocean[table_2+x*DIM+y] = compute_cell_gatherer(x,y,table);
      }
    }
    //switch table after each iteration
    if (table==0)
    {
      table = DIM*DIM;
      table_2 = 0;
    }
    else
    {
      table = 0;
      table_2 = DIM*DIM;
    }
  }
  //we have to be sure that the last calculated table is in the first half of ocean
  if(table != 0)
  {
    for (int x = 1; x < DIM-1; x++)
    {
      for (int y = 1; y < DIM-1; y++)
      {
        ocean[table_2+x*DIM+y] = ocean[table+x*DIM+y];
      }
    }
  }
  return DYNAMIC_COLORING;
}



/*
  *Unwrapped sequential expander
*/
static inline float *compute_seq_multipleline_expander(unsigned iterations)
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
          compute_cell_expander(x,y,div4);
        }
        if(ocean[x*DIM+y+1] >= MAX_HEIGHT)
        {
          int div4 = ocean[x*DIM+y+1]/4;
          compute_cell_expander(x,y+1,div4);
        }
        if(ocean[x*DIM+y+2] >= MAX_HEIGHT)
        {
          int div4 = ocean[x*DIM+y+2]/4;
          compute_cell_expander(x,y+2,div4);
        }
        if(ocean[x*DIM+y+3] >= MAX_HEIGHT)
        {
          int div4 = ocean[x*DIM+y+3]/4;
          compute_cell_expander(x,y+3,div4);
        }
      } 
    }
  }
  return DYNAMIC_COLORING;
}

/*
  *Unwrapped sequential gatherer
*/
static inline float *compute_seq_multipleline_gatherer(unsigned iterations)
{
  if(is_end == true)
  {
    return DYNAMIC_COLORING;
  }
  is_end = true;

  int table = 0;
  int table_2 = DIM*DIM;

  for (unsigned i = 0; i < iterations; i++)
  {
    for (int x = 1; x < DIM-1; x++)
    {
      for (int y = 1; y < DIM-5; y+=4)
      {
        ocean[table_2+x*DIM+y]   = compute_cell_gatherer(x,y,table);
        ocean[table_2+x*DIM+y+1] = compute_cell_gatherer(x,y+1,table);
        ocean[table_2+x*DIM+y+2] = compute_cell_gatherer(x,y+2,table);
        ocean[table_2+x*DIM+y+3] = compute_cell_gatherer(x,y+3,table);

      }
      // boucle nécessaire pour finir le traitement de la matrice si DIM-1 n'est pas multiple de 5
      int res = (DIM-2)%4;
      if(res!= 0)
      {
        for (int y = res; y < DIM-1; y++)
        {
          ocean[table_2+x*DIM+y] = compute_cell_gatherer(x,y,table);
        }
      }
    }

    //switch table after each iteration
    if (table==0)
    {
      table = DIM*DIM;
      table_2 = 0;
    }
    else
    {
      table = 0;
      table_2 = DIM*DIM;
    }

  }
  //we have to be sure that the last calculated table is in the first half of ocean
  if(table != 0)
  {
    for (int x = 1; x < DIM-1; x++)
    {
      for (int y = 1; y < DIM-1; y++)
      {
        ocean[table_2+x*DIM+y] = ocean[table+x*DIM+y];
      }
    }
  }
  return DYNAMIC_COLORING;
}

static inline float *compute_parallel_for_gatherer2(unsigned iterations)
{
  if(is_end == true)
  {
    return DYNAMIC_COLORING;
  }
  is_end = true;

  for (unsigned i = 0; i < iterations; i++)
  {
#pragma omp parallel shared(is_end,ocean)
    {
    unsigned ocean_private[DIM*DIM];
    for(int j = 0 ; j < DIM*DIM ; j++)
      ocean_private[j] = 0;


#pragma omp for schedule(static,4)
        for (int x = 1; x < DIM-1; x++)
        {
          for (int y = 1; y < DIM-1; y++)
          {

            int val = ocean[x*DIM+y]%4;
            val += ocean[(x-1)*DIM+y]/4;
            val += ocean[(x+1)*DIM+y]/4;
            val += ocean[x*DIM+y-1]/4;
            val += ocean[x*DIM+y+1]/4;

            ocean_private[x*DIM+y] = val;

            if(ocean[x*DIM+y] >= MAX_HEIGHT)
              is_end = false;

          }
        }
#pragma omp for schedule(static,4)
      for (int x = 1; x < DIM-1; x++){
        for (int y = 1; y < DIM-1; y++){
          ocean[x*DIM+y] = ocean_private[x*DIM+y];
        }
      }
    }

  }
  return DYNAMIC_COLORING;
}

/** 
  * Compute fonction for collapsed treatment 
*/
static inline float *compute_parallel_for_gatherer(unsigned iterations)
{
  if(is_end == true)
  {
    return DYNAMIC_COLORING;
  }
  is_end = true;

  for (unsigned i = 0; i < iterations; i++)
  {
#pragma omp parallel shared(is_end,ocean)
    {
    unsigned ocean_private[DIM*DIM];
    for(int j = 0 ; j < DIM*DIM ; j++)
      ocean_private[j] = 0;


#pragma omp for schedule(static,4)
        for (int x = 1; x < DIM-1; x++)
        {
          for (int y = 1; y < DIM-1; y++)
          {
          int val = ocean[x*DIM+y]%4 + 
                    ocean[(x-1)*DIM+y]/4 +
                    ocean[(x+1)*DIM+y]/4 +
                    ocean[x*DIM+y-1]/4 +
                    ocean[x*DIM+y+1]/4;

            ocean_private[x*DIM+y] = val;

            if(ocean[x*DIM+y] >= MAX_HEIGHT)
              is_end = false;

          }
        }
#pragma omp for schedule(static,4)
      for (int x = 1; x < DIM-1; x++){
        for (int y = 1; y < DIM-1; y++){
          ocean[x*DIM+y] = ocean_private[x*DIM+y];
        }
      }
    }

  }
  return DYNAMIC_COLORING;
}

static inline float *compute_parallel_p_iteration(unsigned iterations)
{
  if(is_end == true)
  {
    return DYNAMIC_COLORING;
  }
  is_end = true;

  unsigned nb_lines = round((DIM-2) / omp_get_max_threads());

  //nb de lignes spécifiques au dernier thread : il prend la différence.
  unsigned last_thread_lines = nb_lines*omp_get_max_threads() != DIM-2 ? nb_lines+abs(DIM-2 - nb_lines*omp_get_max_threads()) : nb_lines;

#pragma omp parallel shared(is_end, iterations) firstprivate(nb_lines)
  {
    unsigned my_num = omp_get_thread_num();
    int table = 0;

    int begin = ((my_num)*nb_lines)-iterations; 

    unsigned ocean_private[DIM*((my_num==omp_get_max_threads()-1?last_thread_lines:nb_lines)+iterations*2)][2];

    for (int x = 0, final = my_num==omp_get_max_threads()-1?last_thread_lines+iterations*2:nb_lines+iterations*2 ; x < final ; x++)
    {
      for(int y = 0 ; y < DIM ; y++)
      {
        if(begin+x > 0 && begin+x < DIM-1)
        {
          ocean_private[x*DIM+y][table] = ocean[(begin+x/*+1*/)*DIM+y];
        }
        else{
          ocean_private[x*DIM+y][table] = 0;
        }
        ocean_private[x*DIM+y][1-table] = 0;
      }
    }

    if(my_num == 1){
      printf(" --- pre-iteration --- je traite %d\n",begin+iterations);
      for (int x = 0, final = my_num==omp_get_max_threads()-1?last_thread_lines+iterations*2:nb_lines+iterations*2; x < final ; x++)
      {
        for(int y = 0 ; y < DIM ; y++)
        {
          printf("%d ",ocean_private[x*DIM+y][table]);
          printf("%d - ",ocean_private[x*DIM+y][1-table]);
        }
        printf("[%d]\n",x);
      }
    }

    //fin de l'initialisation



    for (unsigned i = 1; i <= iterations; i++)
    {
      int taille = (iterations-i)*2+(my_num==omp_get_max_threads()-1?last_thread_lines:nb_lines);
      for (int x = i; x < i+taille; x++)
      {
        for (int y = 1; y < DIM-1; y++)
        {
          int val = ocean_private[x*DIM+y][table]%4 + 
                    ocean_private[(x-1)*DIM+y][table]/4 +
                    ocean_private[(x+1)*DIM+y][table]/4 +
                    ocean_private[x*DIM+y-1][table]/4 +
                    ocean_private[x*DIM+y+1][table]/4;

          ocean_private[x*DIM+y][1-table] = val;

          if(ocean_private[x*DIM+y][table] >= MAX_HEIGHT)
            is_end = false;
        }
        
      }

      
      table = 1-table;
    }
//tmp
    // if(my_num == 10){
    //   printf(" --- post_iterations ---\n");
    //   for (int x = 0, final = my_num==omp_get_max_threads()-1?last_thread_lines+iterations*2:nb_lines+iterations*2; x < final ; x++)
    //   {
    //     for(int y = 0 ; y < DIM ; y++)
    //     {
    //       printf("%d ",ocean_private[x*DIM+y][table]);
    //       printf("%d - ",ocean_private[x*DIM+y][1-table]);
    //     }
    //     printf("[%d]\n",x);
    //   }
    // }


    begin = my_num*nb_lines;
    for (int x = 0, final = my_num==omp_get_max_threads()-1?last_thread_lines:nb_lines; x < final ; x++)
    {
      for(int y = 1 ; y < DIM-1 ; y++)
      {
        if(begin+x+1 > 0 && begin+x < DIM-1)
        {
#pragma omp critical
          ocean[(begin+x+1)*DIM+y] = ocean_private[(x+iterations)*DIM+y][table];
          // printf("\n- %d - %d -\n",(begin+x+1), ocean_private[(x+iterations)*DIM+y][table]);
        }
      }
    }
  }
  //tmp
    print();
  // exit(0);

  return DYNAMIC_COLORING;
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
          int div4 = ocean[x*DIM+y]/4;
          compute_cell_expander(x,y,div4);
        }
      }
    }
  }
  return DYNAMIC_COLORING;
}

// --------------------------------------------------------------------------------------------------------
// ----------------------------- Initials functions to adapt at arguments  --------------------------------
// --------------------------------------------------------------------------------------------------------

void without_display(float * (*compute_func_t) (unsigned iterations))
{
  // TODO : prendre une décision pour la valeur de iterations par défault.
  do{
    compute_func_t(10);
  }
  while(is_end == false);
}

void expand_ocean(){

  ocean = realloc(ocean,sizeof(unsigned)*DIM*DIM*2);
  for (int x = 0; x < DIM; x++)
  {
    for (int y = 0; y < DIM; y++)
    {
      ocean[DIM*DIM+x*DIM+y] = 0;
    }
  }
}

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

void treatment(int argc, char ** argv)
{
  if(argc)
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
  unsigned long temps;
  struct timeval t1, t2;

  gettimeofday(&t1,NULL);

  // switch case to determine which algorithm the user want
  // the argument wanted is a char, convert in int with its ASCII code
  switch((int)argv[4][0])
  {
    //case sequential expander
    case 115 :  // ASCII for s
      printf("Sequential Expander Algorithm\n");
      //start the non-graphical version
      if (strtol(argv[1],NULL,10)==0) without_display(compute_seq_expander);
      //start the graphical version
      else display_init (argc, argv,
                      DIM,                    // dimensions ( = x = y)
                      MAX_HEIGHT,             // maximum number of grains in the square
                      get,                    // callback func
                      compute_seq_expander);  // callback func
      break;

    //case sequential gatherer
    case 83 :  // ASCII for S
      expand_ocean();
      printf("Sequential Gatherer Algorithm\n");
      //start the non-graphical version
      if (strtol(argv[1],NULL,10)==0) without_display(compute_seq_gatherer);
      //start the graphical version
      else display_init (argc, argv,
                      DIM,                    // dimensions ( = x = y)
                      MAX_HEIGHT,             // maximum number of grains in the square
                      get,                    // callback func
                      compute_seq_gatherer);  // callback func
      break;

    //case sequential unwrapped expander
    case 117 : //ascii u
      printf("Sequential expander Algorithm | 4 lines each loop turn\n");
      //start the non-graphical version
      if (strtol(argv[1],NULL,10)==0) without_display(compute_seq_multipleline_expander);
      //start the graphical version
      else display_init (argc, argv,
                      DIM,                    // dimensions ( = x = y)
                      MAX_HEIGHT,             // maximum number of grains in the square
                      get,                    // callback func
                      compute_seq_multipleline_expander);  // callback func
      break;

    //case sequential unwrapped gatherer
    case 85 : //ascii U
      expand_ocean();
      printf("Sequential gatherer Algorithm | 4 lines each loop turn\n");
      //start the non-graphical version
      if (strtol(argv[1],NULL,10)==0) without_display(compute_seq_multipleline_gatherer);
      //start the graphical version
      else display_init (argc, argv,
                      DIM,                                 // dimensions ( = x = y)
                      MAX_HEIGHT,                          // maximum number of grains in the square
                      get,                                 // callback func
                      compute_seq_multipleline_gatherer);  // callback func
      break;

    //case parallel for gatherer
    case 70 :  // ASCII for F
      printf("Parallel Gatherer Algorithm\n");
      if (strtol(argv[1],NULL,10)==0) without_display(compute_parallel_for_gatherer);
      //start the graphical version
      else display_init (argc, argv,
                      DIM,                                 // dimensions ( = x = y)
                      MAX_HEIGHT,                          // maximum number of grains in the square
                      get,                                 // callback func
                      compute_parallel_for_gatherer);      // callback func
      break;

    //case parallel gatherer which synchronise each p iterations
    case 105 : //ascii i
      printf("Parallel Gatherer P iterations Algorithm\n");
      //start the non-graphical version
      if (strtol(argv[1],NULL,10)==0) without_display(compute_parallel_p_iteration);
      //start the graphical version
      else display_init (argc, argv,
                      DIM,                                 // dimensions ( = x = y)
                      MAX_HEIGHT,                          // maximum number of grains in the square
                      get,                                 // callback func
                      compute_parallel_p_iteration);       // callback func
      break;
    
    //case parallel tasked gatherer
    case 116 : //ascii t
      printf("Parallel Task Expander Algorithm\n");
      //start the non-graphical version
      if (strtol(argv[1],NULL,10)==0) without_display(compute_parallel_task);
      //start the graphical version
      else display_init (argc, argv,
                      DIM,                                 // dimensions ( = x = y)
                      MAX_HEIGHT,                          // maximum number of grains in the square
                      get,                                 // callback func
                      compute_parallel_task);              // callback func
      break;

    default :
      printf("Unrecognize Algorithm. Please, read how to use this program\n\n");
      printf("%s\n", man);
      break;
  }
  gettimeofday(&t2,NULL);
   
  temps = TIME_DIFF(t1,t2);
  printf("Algorithm time = %ld.%03ldms \n", temps/1000, temps%1000);

  //print();
  free(ocean);
}
