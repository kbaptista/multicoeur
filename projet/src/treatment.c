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
  * Divide the [x][y] cell content then add to himself the quarter of neigbours squares.
*/
static inline unsigned compute_cell_gatherer(int x, int y, int table)
{
  unsigned val = ocean[table+x*DIM+y]%4 + ocean[table+(x-1)*DIM+y]/4 + ocean[table+(x+1)*DIM+y]/4 + ocean[table+x*DIM+y-1]/4 + ocean[table+x*DIM+y+1]/4;
  is_end &= (val==ocean[table+x*DIM+y]) ;
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
 
  int table = DIM*DIM;
  int index = 0;
 
  for (unsigned i = 0; i < iterations; i++)
  {
    for (int x = 1; x < DIM-1; x++)
    {
      for (int y = 1; y < DIM-1; y++)
      {
        ocean[table*(1-index)+x*DIM+y] = compute_cell_gatherer(x,y,table*index);
      }
    }
    //switch table after each iteration
    index = 1- index;
  }
  //we have to be sure that the last calculated table is in the first half of ocean
  if(table*(1-index) != 0)
  {
    for (int x = 1; x < DIM-1; x++)
    {
      for (int y = 1; y < DIM-1; y++)
      {
        ocean[table*(1-index)+x*DIM+y] = ocean[table*index+x*DIM+y];
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
      // boucle nécessaire pour finir le traitement de la matrice si DIM-1 n'est pas multiple de 4 
      int res = (DIM-2)%4;
      if(res!= 0)
      {
        for (int y = res; y < DIM-1; y++)
        {
          int div4 = ocean[x*DIM+y]/4;
          compute_cell_expander(x,y,div4);
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


void compute_parallel_p_iteration_inside(int iterations, int nb_lines, int my_lines, int my_num)
{
  int table = 0;
  int begin = my_num*nb_lines-iterations; 
  unsigned ocean_private[DIM*(my_lines+iterations*2)][2];
  unsigned taille = my_lines+iterations*2;

  for (int x = 0 ; x < taille ; x++)
  {
    for(int y = 0 ; y < DIM ; y++)
    {
      if(begin+x > 0 && begin+x < DIM-1)
      {
        ocean_private[x*DIM+y][table] = ocean[(begin+x+1)*DIM+y];
      }
      else{
        ocean_private[x*DIM+y][table] = 0;
      }
      ocean_private[x*DIM+y][1-table] = 0;
    }
  }

#pragma omp barrier
  //fin de l'initialisation

  for (unsigned i = 1; i <= iterations; i++)
  {
    taille = (iterations-i)*2+my_lines;
    //tmp
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

  begin = begin+iterations;
  for (int x = 0; x < my_lines ; x++)
  {
    for(int y = 1 ; y < DIM-1 ; y++)
    {
#pragma omp critical
      ocean[(begin+x+1)*DIM+y] = ocean_private[(iterations+x)*DIM+y][1-table];
    }
  }
}

/** 
  * Compute fonction for p iteration synchronization.
*/
static inline float *compute_parallel_p_iteration(unsigned iterations)
{
  if(is_end == true)
  {
    return DYNAMIC_COLORING;
  }
  is_end = true;

  unsigned nb_lines = round((DIM-2) / omp_get_max_threads());

  //nombre de lignes spécifiques au dernier thread : il prend/retire la différence.
  unsigned last_thread_lines = nb_lines*omp_get_max_threads() != DIM-2 ? nb_lines+ DIM-2 - nb_lines*omp_get_max_threads() : nb_lines;

#pragma omp parallel shared(is_end, iterations,nb_lines)
  {
    unsigned my_num = omp_get_thread_num();
    compute_parallel_p_iteration_inside(iterations,nb_lines,my_num==omp_get_max_threads()-1?last_thread_lines:nb_lines,my_num);
  }

  return DYNAMIC_COLORING;
}

/** 
  * Compute fonction with task treatment
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

/*
  function used when we don't use display to simulate all it calls.
*/
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



// void speedup_inside_init(int algo, int config, int iter)
// {
//   switch(algo)
//   {
//     case 0 : break;
//     case 1 : 
//     case 3 : expand_ocean(); break;
//     case 2 : 
//     case 4 : ocean = realloc(ocean,sizeof(unsigned)*DIM*DIM); break;
//     case 5 : break;
//     case 6 : break;
//   }

//   if(config == 0) sand_init_center();
//   else sand_init_homogeneous();
// }

// /*
//   calculate speed up for every algorithm on every size and configuration
// */
// void speedup(int iterations)
// {
//   int nb_configurations = 2;
//   int nb_algorithmes = 7;
//   int sizes[2]={128,512};
//   int nb_sizes = 2;

//   // 7 algo * 2 config * 2 size
//   double resultats[nb_algorithmes*nb_configurations*nb_sizes];
//   for(int i = 0 ; i < nb_algorithmes*nb_configurations*nb_sizes ; i++)
//   {
//     resultats[i] = 0;
//   }

//   unsigned long temps;
//   struct timeval t1, t2;



//   ocean = malloc(sizeof(unsigned)*sizes[0]*sizes[0]);

//   for(int algo = 0 ; algo < nb_algorithmes ; algo++)
//   {
//     for(int s = 0 ; s < nb_sizes ; s++)
//     {
//       DIM = sizes[s];

//       for(int c = 0 ; c < nb_configurations ; c++)
//       {
//         for(int iter = 0 ; iter < iterations ; iter++)
//         {

//           printf("");
//           speedup_inside_init(algo,c,iter);

//           gettimeofday(&t1,NULL);
          
//           switch(algo)
//           {
//             case 0 : without_display(compute_seq_expander); break;
//             case 1 : without_display(compute_seq_gatherer); break;
//             case 2 : without_display(compute_seq_multipleline_expander); break;
//             case 3 : without_display(compute_seq_multipleline_gatherer); break;
//             case 4 : without_display(compute_parallel_for_gatherer); break;
//             case 5 : without_display(compute_parallel_p_iteration); break;
//             case 6 : without_display(compute_parallel_task); break;
//           }

//           gettimeofday(&t2,NULL);
          
//           resultats[algo+s+c] += TIME_DIFF(t1,t2);
//         }
//         resultats[algo+s+c] = resultats[algo+s+c]/iterations;
//         printf("%f\n",resultats[algo+s+c]);
//       }
//     }
//     printf("algo : "+algo);
//   }


//   free(ocean);
//   exit(0);
// }

static char * man = {"usage : ./sand <GUI> <INITIALIZATION> <SIZE> <ALGORITHM>\n"
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
    "\t\t- F : it runs the parallel gatherer method ;\n"
    "\t\t- P : it runs the parallel gatherer method which synchronise each p iterations ;\n"
    "\t\t- t : it runs the parallel expander task method ;\n"
  };

void treatment(int argc, char ** argv)
{

  if(argc < 5) {
      printf("%s", man);
      exit(-1);
    }

  // //branchement vers speedup
  // if(strtol(argv[1],NULL,10)==2)
  //   speedup(strtol(argv[2],NULL,10));

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
    case 80 : //ascii P
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
      printf("Unrecognize Algorithm. Please, use an algorithm as you can see in the usage : ./sand\n\n");
      printf("%s\n", man);
      break;
  }
  gettimeofday(&t2,NULL);
   
  temps = TIME_DIFF(t1,t2);
  printf("Algorithm time = %ld.%03ldms \n", temps/1000, temps%1000);
  //Set the border at 0 (help to track differences between sequential expander and others algorithms)
  for (int i = 0; i < DIM; ++i)
  {
    ocean[i]=0;             //Top Line
    ocean[i+DIM*(DIM-1)]=0; //Bottom Line
    ocean[i*DIM]=0;         //left Line
    ocean[i*DIM+DIM-1]=0;   //Right Line
  }
  print();
  free(ocean);
}

int main (int argc, char **argv)
{ 
  treatment(argc,argv);
  return 0;
}