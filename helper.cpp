/* 
 * Utilities for the Aliev-Panfilov code
 * Scott B. Baden, UCSD
 * Nov 2, 2015
 */

#include <iostream>
#include <assert.h>
// Needed for memalign
#include <malloc.h>
#include <math.h>

#ifdef _MPI_
#include <mpi.h>
#endif
#include "cblock.h"
using namespace std;

extern control_block cb;

void printMat(const char mesg[], double *E, int m, int n);

//Used by each worker to compute the statistics for their partition only
void stats_block(double *E, int m, int n, double *_mx, double *sumSq, int row_skip){
    double mx = -1;
    double _sumSq = 0;
    int i, j;

    for (j=0; j<m; j++) {
      for(i=0;i<n;i++)
      {

        _sumSq += E[j*row_skip+i]*E[j*row_skip+i];
        double fe = fabs(E[j*row_skip+i]);
        if (fe > mx)
          mx = fe;

      }
    }
    *_mx = mx;
    *sumSq = _sumSq;
}

//
// Initialization
//
// We set the right half-plane of E_prev to 1.0, the left half plane to 0
// We set the bottom half-plane of R to 1.0, the top half plane to 0
// These coordinates are in world (global) coordinate and must
// be mapped to appropriate local indices when parallelizing the code
//
void init (double *E,double *E_prev,double *R,int m,int n){
  #ifndef _MPI_
    int i;

    for (i=0; i < (m+2)*(n+2); i++)
        E_prev[i] = R[i] = 0;

    for (i = (n+2); i < (m+1)*(n+2); i++) {
    int colIndex = i % (n+2);		// gives the base index (first row's) of the current index

        // Need to compute (n+1)/2 rather than n/2 to work with odd numbers
    if(colIndex == 0 || colIndex == (n+1) || colIndex < ((n+1)/2+1))
      continue;

        E_prev[i] = 1.0;
    }

    for (i = 0; i < (m+2)*(n+2); i++) {
    int rowIndex = i / (n+2);		// gives the current row number in 2D array representation
    int colIndex = i % (n+2);		// gives the base index (first row's) of the current index

        // Need to compute (m+1)/2 rather than m/2 to work with odd numbers
    if(colIndex == 0 || colIndex == (n+1) || rowIndex < ((m+1)/2+1))
      continue;

        R[i] = 1.0;
    }
  #endif
  #ifdef _MPI_
    int nprocs,myrank;
    MPI_Comm_size(MPI_COMM_WORLD,&nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD,&myrank);
    int x=cb.px;
    int y=cb.py;
    int idx=myrank/y;
    int idy=myrank-y*idx;
    //check if x,y divide m,n properly. and assign one extra row to each process with idx<m%x. assign one extra column to each process with idy<n%y;
    int mod_x=m%x;
    int mod_y=n%y;
    MPI_Request req_e[1000],req_r[1000];
    if(myrank==0)//initialize the whole matrix and send each portion to the corresponding process
    {
      
      int i;

      for (i=0; i < (m+2)*(n+2); i++)
          E_prev[i] = R[i] = 0;

      for (i = (n+2); i < (m+1)*(n+2); i++) {
      int colIndex = i % (n+2);   // gives the base index (first row's) of the current index

          // Need to compute (n+1)/2 rather than n/2 to work with odd numbers
      if(colIndex == 0 || colIndex == (n+1) || colIndex < ((n+1)/2+1))
        continue;

          E_prev[i] = 1.0;
      }

      for (i = 0; i < (m+2)*(n+2); i++) {
      int rowIndex = i / (n+2);   // gives the current row number in 2D array representation
      int colIndex = i % (n+2);   // gives the base index (first row's) of the current index

          // Need to compute (m+1)/2 rather than m/2 to work with odd numbers
      if(colIndex == 0 || colIndex == (n+1) || rowIndex < ((m+1)/2+1))
        continue;

          R[i] = 1.0;
      }
     
     
      MPI_Request dummy_req,another_dummy_req;
      for(int rank=1;rank<nprocs;rank++)  //send each process the correct portion 
      {   
          int i=rank/y;
          int j=rank-y*i;
          
          int added_x_prev=mod_x<i?mod_x:i;
          int added_y_prev=mod_y<j?mod_y:j;
          int to_add_x=added_x_prev<mod_x?1:0;
          int to_add_y=added_y_prev<mod_y?1:0;

          int row_start=i*(m/x)+added_x_prev;
          int row_end=row_start+(m/x)+to_add_x-1;

          int col_start=j*(n/y)+added_y_prev;
          int col_end=col_start+(n/y)+to_add_y-1;

          int innerBlockRowStartIndex = (n+2)*(row_start+1)+1;
          int innerBlockRowEndIndex = (n+2)*(row_end+1)+1;

          int row_count=row_end-row_start+1;
          
          int stride=n+2;
          int blocklength=col_end-col_start+1;
          MPI_Datatype final_result1,final_result2;
          MPI_Type_vector(row_count,blocklength,stride,MPI_DOUBLE,&final_result1);
          MPI_Type_vector(row_count,blocklength,stride,MPI_DOUBLE,&final_result2);

          MPI_Type_commit(&final_result1);
          MPI_Type_commit(&final_result2);

          MPI_Send((E_prev+innerBlockRowStartIndex+col_start), 1, final_result1, rank, 0, MPI_COMM_WORLD);
          MPI_Send((R+innerBlockRowStartIndex+col_start), 1, final_result2, rank, 1, MPI_COMM_WORLD);
      }
    }
    else //receive the correct portion from process 0
    {  
      int i=myrank/y;
      int j=myrank-y*i;
      
      int added_x_prev=mod_x<i?mod_x:i;
      int added_y_prev=mod_y<j?mod_y:j;
      int to_add_x=added_x_prev<mod_x?1:0;
      int to_add_y=added_y_prev<mod_y?1:0;

      int row_start=i*(m/x)+added_x_prev;
      int row_end=row_start+(m/x)+to_add_x-1;

      int col_start=j*(n/y)+added_y_prev;
      int col_end=col_start+(n/y)+to_add_y-1;

      int innerBlockRowStartIndex = (n+2)*(row_start+1)+1;
      int innerBlockRowEndIndex = (n+2)*(row_end+1)+1;

      int row_count=row_end-row_start+1;
      
      int stride=n+2;
      int blocklength=col_end-col_start+1;

      MPI_Datatype final_result1,final_result2;
      MPI_Type_vector(row_count,blocklength,stride,MPI_DOUBLE,&final_result1);
      MPI_Type_vector(row_count,blocklength,stride,MPI_DOUBLE,&final_result2);

      MPI_Type_commit(&final_result1);
      MPI_Type_commit(&final_result2);
      
      MPI_Recv((E_prev+innerBlockRowStartIndex+col_start), 1, final_result1, 0, 0, MPI_COMM_WORLD,MPI_STATUS_IGNORE);
      MPI_Recv((R+innerBlockRowStartIndex+col_start), 1, final_result2, 0, 1, MPI_COMM_WORLD,MPI_STATUS_IGNORE);
    }
    //wiat for all E and all R to be received
    
  #endif
    // We only print the meshes if they are small enough
#if 0
    printMat("E_prev",E_prev,m,n);
    printMat("R",R,m,n);
#endif
}

double *alloc1D(int m,int n){
    int nx=n, ny=m;
    double *E;
    // Ensures that allocated memory is aligned on a 16 byte boundary
    assert(E= (double*) memalign(16, sizeof(double)*nx*ny) );
    return(E);
}

void printMat(const char mesg[], double *E, int m, int n){
    int i;
#if 0
    if (m>8)
      return;
#else
    if (m>34)
      return;
#endif
    printf("%s\n",mesg);
    for (i=0; i < (m+2)*(n+2); i++){
       int rowIndex = i / (n+2);
       int colIndex = i % (n+2);
       if ((colIndex>0) && (colIndex<n+1))
          if ((rowIndex > 0) && (rowIndex < m+1))
            printf("%6.3f ", E[i]);
       if (colIndex == n+1)
	    printf("\n");
    }
}
