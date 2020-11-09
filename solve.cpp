/* 
 * Solves the Aliev-Panfilov model  using an explicit numerical scheme.
 * Based on code orginally provided by Xing Cai, Simula Research Laboratory
 * 
 * Modified and  restructured by Scott B. Baden, UCSD
 * 
 */


//just to check!
#include <assert.h>
#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <string>
#include <math.h>
#include "time.h"
#include "apf.h"
#include "Plotting.h"
#include "cblock.h"
#include <emmintrin.h>
#ifdef _MPI_
#include <mpi.h>
#endif
using namespace std;

void repNorms(double l2norm, double mx, double dt, int m,int n, int niter, int stats_freq);
void stats(double *E, int m, int n, double *_mx, double *sumSq);
void stats_block(double *E, int m, int n, double *_mx, double *sumSq, int row_skip);
extern control_block cb;

#ifdef SSE_VEC
// If you intend to vectorize using SSE instructions, you must
// disable the compiler's auto-vectorizer
__attribute__((optimize("no-tree-vectorize")))
#endif 

// The L2 norm of an array is computed by taking sum of the squares
// of each element, normalizing by dividing by the number of points
// and then taking the sequare root of the result
//
double L2Norm(double sumSq){
    double l2norm = sumSq /  (double) (cb.m*cb.n);
    l2norm = sqrt(l2norm);
    return l2norm;
}

void solve(double **_E, double **_E_prev, double *R, double alpha, double dt, Plotter *plotter, double &L2, double &Linf){
#ifdef _MPI_
 // Simulated time is different from the integer timestep number
    double t = 0.0;

    double *E = *_E, *E_prev = *_E_prev;
    double *R_tmp = R;
    double *E_tmp = *_E;
    double *E_prev_tmp = *_E_prev;
    double mx, sumSq;
    int niter;
    int m = cb.m, n=cb.n;


    int nprocs;
    MPI_Comm_size(MPI_COMM_WORLD,&nprocs);
    int myrank;
    MPI_Comm_rank(MPI_COMM_WORLD,&myrank);
    int x=cb.px;
    int y=cb.py;
    int idx=myrank/y;
    int idy=myrank-y*idx;


    //check if x,y divide m,n properly. and assign one extra row to each process with idx<m%x. assign one extra column to each process with idy<n%y;
    int mod_x=m%x;
    int mod_y=n%y;
    int added_x_prev=mod_x<idx?mod_x:idx;
    int added_y_prev=mod_y<idy?mod_y:idy;
    int to_add_x=added_x_prev<mod_x?1:0;
    int to_add_y=added_y_prev<mod_y?1:0;

    int my_row_start=idx*(m/x)+added_x_prev;
    int my_row_end=my_row_start+(m/x)+to_add_x-1;

    int my_col_start=idy*(n/y)+added_y_prev;
    int my_col_end=my_col_start+(n/y)+to_add_y-1;

    int innerBlockRowStartIndex = (n+2)*(my_row_start+1)+1;
    int innerBlockRowEndIndex = (n+2)*(my_row_end+1)+1;

    
    for (niter = 0; niter < cb.niters; niter++)
    {
        
        if  (cb.debug && (niter==0)){
            stats(E_prev,m,n,&mx,&sumSq);
            double l2norm = L2Norm(sumSq);
            repNorms(l2norm,mx,dt,m,n,-1, cb.stats_freq);
            if (cb.plot_freq)
                plotter->updatePlot(E,  -1, m+1, n+1);
        }
        
        /*
         * Copy data from boundary of the computational box to the
         * padding region, set up for differencing computational box's boundary
         *
         * These are physical boundary conditions, and are not to be confused
         * with ghost cells that we would use in an MPI implementation
         *
         * The reason why we copy boundary conditions is to avoid
         * computing single sided differences at the boundaries
         * which increase the running time of solve()
         *
         */
        
        // 4 FOR LOOPS set up the padding needed for the boundary conditions
        int i,j;
        
        
        /*
         * The following section handles communication with up to four adjacent processes
         * The counter keeps track of the number of requests that need to be completed
         * If the if statements are not satisfied, the program reverts to filling in the boundaries using values within the block
         */
        
        MPI_Request request[8];
        int counter = 0;
        
        int row_count=my_row_end-my_row_start+1;
        int stride=n+2;
        int blocklength=1;
        MPI_Datatype column_jumper;

        MPI_Type_vector(row_count,blocklength,stride,MPI_DOUBLE,&column_jumper);
        MPI_Type_commit(&column_jumper);

        if ((idy > 0)&&(!cb.noComm))
        {
            MPI_Isend((E_prev+innerBlockRowStartIndex+my_col_start), 1, column_jumper, myrank-1, 0, MPI_COMM_WORLD, &request[counter]);
            MPI_Irecv((E_prev+innerBlockRowStartIndex+my_col_start-1), 1, column_jumper, myrank-1, 0, MPI_COMM_WORLD, &request[counter+1]);
            counter += 2;
        }
        else
        {
            // Fills in the LEFT Ghost Cells
            for (i = my_row_start+1; i < (my_row_end+3)*(n+2); i+=(n+2)) {
                E_prev[i] = E_prev[i+2];
            }
        }
        
        
        if ((idy < y-1)&&(!cb.noComm))
        {
            MPI_Isend((E_prev+innerBlockRowStartIndex+my_col_end), 1, column_jumper, myrank+1, 0, MPI_COMM_WORLD, &request[counter]);
            MPI_Irecv((E_prev+innerBlockRowStartIndex+my_col_end+1), 1, column_jumper, myrank+1, 0, MPI_COMM_WORLD, &request[counter+1]);
            counter += 2;
            
        }
        else
        {
            // Fills in the RIGHT Ghost Cells
            for (i = (n+1)+my_row_start*(n+2); i < (my_row_end+3)*(n+2); i+=(n+2)) {
                E_prev[i] = E_prev[i-2];
            }
            
        }
        
        
        if ((idx > 0)&&(!cb.noComm))
        {
            j=innerBlockRowStartIndex+my_col_start;//j points to the begining of the top-most inner row
            MPI_Isend((E_prev+j), my_col_end-my_col_start+1, MPI_DOUBLE, myrank-y, 0, MPI_COMM_WORLD, &request[counter]);
            
            j=innerBlockRowStartIndex-(n+2)+my_col_start;//j points to the begining of the top-most ghost row
            MPI_Irecv((E_prev+j), my_col_end-my_col_start+1, MPI_DOUBLE, myrank-y, 0, MPI_COMM_WORLD, &request[counter+1]);
            counter += 2;
        }
        else
        {
            // Fills in the TOP Ghost Cells
            for (i = 0; i < (n+2); i++) {
                E_prev[i] = E_prev[i + (n+2)*2];
            }
        }
        
        
        if ((idx < x-1)&&(!cb.noComm))
        {
            j=innerBlockRowEndIndex+my_col_start;//j points to the begining of the bottom-most inner row
            MPI_Isend((E_prev+j), my_col_end-my_col_start+1, MPI_DOUBLE, myrank+y, 0, MPI_COMM_WORLD, &request[counter]);
            
            j=innerBlockRowEndIndex+(n+2)+my_col_start;//j points to the begining of the bottom-most ghost row
            MPI_Irecv((E_prev+j), my_col_end-my_col_start+1, MPI_DOUBLE, myrank+y, 0, MPI_COMM_WORLD, &request[counter+1]);
            counter += 2;
        }
        else
        {
            // Fills in the BOTTOM Ghost Cells
            for (i = ((m+2)*(n+2)-(n+2)); i < (m+2)*(n+2); i++) {
                E_prev[i] = E_prev[i - (n+2)*2];
            }
            
        }
        
        if(!cb.noComm)
            MPI_Waitall(counter, request, MPI_STATUSES_IGNORE);
        

        
        //////////////////////////////////////////////////////////////////////////////
        
        #define FUSED 1
        
        #ifdef FUSED
        // Solve for the excitation, a PDE
        //send the necessary gost cells to the neighbors:
        
        for(j = innerBlockRowStartIndex; j <= innerBlockRowEndIndex; j+=(n+2)) {
            E_tmp = E + j;
            E_prev_tmp = E_prev + j;
            R_tmp = R + j;
            for(i = my_col_start; i <= my_col_end; i++) {
                E_tmp[i] = E_prev_tmp[i]+alpha*(E_prev_tmp[i+1]+E_prev_tmp[i-1]-4*E_prev_tmp[i]+E_prev_tmp[i+(n+2)]+E_prev_tmp[i-(n+2)]);
                E_tmp[i] += -dt*(kk*E_prev_tmp[i]*(E_prev_tmp[i]-a)*(E_prev_tmp[i]-1)+E_prev_tmp[i]*R_tmp[i]);
                R_tmp[i] += dt*(epsilon+M1* R_tmp[i]/( E_prev_tmp[i]+M2))*(-R_tmp[i]-kk*E_prev_tmp[i]*(E_prev_tmp[i]-b-1));
            }
        }
        #else
        // Solve for the excitation, a PDE
        for(j = innerBlockRowStartIndex; j <= innerBlockRowEndIndex; j+=(n+2)) {
            E_tmp = E + j;
            E_prev_tmp = E_prev + j;
            for(i = 0; i < n; i++) {
                E_tmp[i] = E_prev_tmp[i]+alpha*(E_prev_tmp[i+1]+E_prev_tmp[i-1]-4*E_prev_tmp[i]+E_prev_tmp[i+(n+2)]+E_prev_tmp[i-(n+2)]);
            }
        }
        
        /*
         * Solve the ODE, advancing excitation and recovery variables
         *     to the next timtestep
         */
        
        for(j = innerBlockRowStartIndex; j <= innerBlockRowEndIndex; j+=(n+2)) {
            E_tmp = E + j;
            R_tmp = R + j;
            E_prev_tmp = E_prev + j;
            for(i = 0; i < n; i++) {
                E_tmp[i] += -dt*(kk*E_prev_tmp[i]*(E_prev_tmp[i]-a)*(E_prev_tmp[i]-1)+E_prev_tmp[i]*R_tmp[i]);
                R_tmp[i] += dt*(epsilon+M1* R_tmp[i]/( E_prev_tmp[i]+M2))*(-R_tmp[i]-kk*E_prev_tmp[i]*(E_prev_tmp[i]-b-1));
            }
        }
        #endif
        /////////////////////////////////////////////////////////////////////////////////
        
        if (cb.stats_freq){
            if ( !(niter % cb.stats_freq)){
                stats(E,m,n,&mx,&sumSq);
                double l2norm = L2Norm(sumSq);
                repNorms(l2norm,mx,dt,m,n,niter, cb.stats_freq);
            }
        }
        
        if (cb.plot_freq){
            if (!(niter % cb.plot_freq)){
                plotter->updatePlot(E,  niter, m, n);
            }
        }
        
        // Swap current and previous meshes
        double *tmp = E; E = E_prev; E_prev = tmp;
        
    } //end of 'niter' loop at the beginning

    
    //Handles the calculates the Linf and sum of squares for this partition 
    
    stats_block(E_prev+innerBlockRowStartIndex+my_col_start,  my_row_end-my_row_start+1,  my_col_end-my_col_start+1,  &Linf, &sumSq, n+2);
    
    /*
     * This handles the communication for the distributed statistics calculation
     * All worker processes send their result back to process 0
     * Proccess 0 computes the L2 norm and identifies the global Linf
     */
    L2 = sumSq;
    if (myrank==0)
    {
        MPI_Request requests[1000];
        double Linfs[1000];
        double sumSqs[1000];
        for (int i=1; i<nprocs; i++)
        {
            MPI_Irecv(&Linfs[i], 1, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, &requests[2*(i-1)]);
            MPI_Irecv(&sumSqs[i], 1, MPI_DOUBLE, i, 1, MPI_COMM_WORLD, &requests[2*(i-1)+1]);
        }
        MPI_Waitall(2*(nprocs-1), requests, MPI_STATUSES_IGNORE);
        for (int i=1; i<nprocs; i++)
        {
            if (Linfs[i] > Linf)
                Linf = Linfs[i];
            L2 += sumSqs[i];
        }
        L2 = L2Norm(L2);
    }
    else
    {
        MPI_Request requests[2];
        MPI_Isend(&Linf, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &requests[0]);
        MPI_Isend(&sumSq, 1, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD, &requests[1]);
        MPI_Waitall(2, requests, MPI_STATUSES_IGNORE);
    }
    
    // Swap pointers so we can re-use the arrays
    *_E = E;
    *_E_prev = E_prev;
#endif

#ifndef _MPI_
    double t = 0.0;
    double *E = *_E, *E_prev = *_E_prev;
    double *R_tmp = R;
    double *E_tmp = *_E;
    double *E_prev_tmp = *_E_prev;
    double mx, sumSq;
    int niter;
    int m = cb.m, n=cb.n;
    int innerBlockRowStartIndex = (n+2)+1;
    int innerBlockRowEndIndex = (((m+2)*(n+2) - 1) - (n)) - (n+2);


     // We continue to sweep over the mesh until the simulation has reached
     // the desired number of iterations
      for (niter = 0; niter < cb.niters; niter++){
      
          if  (cb.debug && (niter==0)){
          stats(E_prev,m,n,&mx,&sumSq);
              double l2norm = L2Norm(sumSq);
          repNorms(l2norm,mx,dt,m,n,-1, cb.stats_freq);
          if (cb.plot_freq)
              plotter->updatePlot(E,  -1, m+1, n+1);
          }

       /* 
        * Copy data from boundary of the computational box to the
        * padding region, set up for differencing computational box's boundary
        *
        * These are physical boundary conditions, and are not to be confused
        * with ghost cells that we would use in an MPI implementation
        *
        * The reason why we copy boundary conditions is to avoid
        * computing single sided differences at the boundaries
        * which increase the running time of solve()
        *
        */
        
        // 4 FOR LOOPS set up the padding needed for the boundary conditions
        int i,j;

        // Fills in the TOP Ghost Cells
        for (i = 0; i < (n+2); i++) {
            E_prev[i] = E_prev[i + (n+2)*2];
        }

        // Fills in the RIGHT Ghost Cells
        for (i = (n+1); i < (m+2)*(n+2); i+=(n+2)) {
            E_prev[i] = E_prev[i-2];
        }

        // Fills in the LEFT Ghost Cells
        for (i = 0; i < (m+2)*(n+2); i+=(n+2)) {
            E_prev[i] = E_prev[i+2];
        }   

        // Fills in the BOTTOM Ghost Cells
        for (i = ((m+2)*(n+2)-(n+2)); i < (m+2)*(n+2); i++) {
            E_prev[i] = E_prev[i - (n+2)*2];
        }

    //////////////////////////////////////////////////////////////////////////////

    #define FUSED 1

    #ifdef FUSED
        // Solve for the excitation, a PDE
        for(j = innerBlockRowStartIndex; j <= innerBlockRowEndIndex; j+=(n+2)) {
            E_tmp = E + j;
        E_prev_tmp = E_prev + j;
            R_tmp = R + j;
        for(i = 0; i < n; i++) {
            E_tmp[i] = E_prev_tmp[i]+alpha*(E_prev_tmp[i+1]+E_prev_tmp[i-1]-4*E_prev_tmp[i]+E_prev_tmp[i+(n+2)]+E_prev_tmp[i-(n+2)]);
                E_tmp[i] += -dt*(kk*E_prev_tmp[i]*(E_prev_tmp[i]-a)*(E_prev_tmp[i]-1)+E_prev_tmp[i]*R_tmp[i]);
                R_tmp[i] += dt*(epsilon+M1* R_tmp[i]/( E_prev_tmp[i]+M2))*(-R_tmp[i]-kk*E_prev_tmp[i]*(E_prev_tmp[i]-b-1));
            }
        }
    #else
        // Solve for the excitation, a PDE
        for(j = innerBlockRowStartIndex; j <= innerBlockRowEndIndex; j+=(n+2)) {
            E_tmp = E + j;
                E_prev_tmp = E_prev + j;
                for(i = 0; i < n; i++) {
                    E_tmp[i] = E_prev_tmp[i]+alpha*(E_prev_tmp[i+1]+E_prev_tmp[i-1]-4*E_prev_tmp[i]+E_prev_tmp[i+(n+2)]+E_prev_tmp[i-(n+2)]);
                }
        }

        /* 
         * Solve the ODE, advancing excitation and recovery variables
         *     to the next timtestep
         */

        for(j = innerBlockRowStartIndex; j <= innerBlockRowEndIndex; j+=(n+2)) {
            E_tmp = E + j;
            R_tmp = R + j;
        E_prev_tmp = E_prev + j;
            for(i = 0; i < n; i++) {
          E_tmp[i] += -dt*(kk*E_prev_tmp[i]*(E_prev_tmp[i]-a)*(E_prev_tmp[i]-1)+E_prev_tmp[i]*R_tmp[i]);
          R_tmp[i] += dt*(epsilon+M1* R_tmp[i]/( E_prev_tmp[i]+M2))*(-R_tmp[i]-kk*E_prev_tmp[i]*(E_prev_tmp[i]-b-1));
            }
        }
    #endif
         /////////////////////////////////////////////////////////////////////////////////

       if (cb.stats_freq){
         if ( !(niter % cb.stats_freq)){
            stats(E,m,n,&mx,&sumSq);
            double l2norm = L2Norm(sumSq);
            repNorms(l2norm,mx,dt,m,n,niter, cb.stats_freq);
        }
       }

       if (cb.plot_freq){
              if (!(niter % cb.plot_freq)){
            plotter->updatePlot(E,  niter, m, n);
            }
        }

       // Swap current and previous meshes
       double *tmp = E; E = E_prev; E_prev = tmp;

     } //end of 'niter' loop at the beginning

      // return the L2 and infinity norms via in-out parameters
      stats(E_prev,m,n,&Linf,&sumSq);
      L2 = L2Norm(sumSq);

      // Swap pointers so we can re-use the arrays
      *_E = E;
      *_E_prev = E_prev;
#endif
}