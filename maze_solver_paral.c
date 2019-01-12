/*
Student Name: Behlülcan Mert ÇOTUK
Student Number: 2011400294
Compile Status: Compiling
Program Status: Working
Notes: There might be differences in whitespaces between the expected output and the program's output.
If you use "diff" command when comparing outputs, it is recommended to use "diff -b" to ignore whitespaces.
*/

#include <stdio.h>
#include <stdbool.h>
#include "mpi.h"
#include <stdlib.h>

#define send_data_tag 2001
#define return_data_tag 2002


// counts the number of wall neighbours of a cell
int wall_neighbours(int rows, int cols, int A[][cols], int i, int j) {
    
    // number of walls around the cell
    int walls_around = 0;

    // check up
    if((i-1)>=0) {
        if(A[i-1][j]==0) walls_around++;
    }

    // check down
    if((i+1)<rows) {
        if(A[i+1][j]==0) walls_around++;
    }

    // check left
    if((j-1)>=0) {
        if(A[i][j-1]==0) walls_around++;
    }

    // check right
    if((j+1)<cols) {
        if(A[i][j+1]==0) walls_around++;
    }

    return walls_around;
}

// counts the number of wall neighbours of a cell also for boundary conditions
int wall_neighbours_boundary(int rows, int cols, int A[][cols], int i, int j
                    , int proc_id, int num_procs, int upper_neighbour[cols], int lower_neighbour[cols]) {
    
    // number of walls around the cell
    int walls_around = 0;

    // check upwards
    if((i-1)>=0) {
        if(A[i-1][j]==0) walls_around++;
    } else {
        if(proc_id > 1) {
            if(upper_neighbour[j]==0) walls_around++;
        }
    }

    // check downwards
    if((i+1)<rows) {
        if(A[i+1][j]==0) walls_around++;
    } else {
        if(proc_id < num_procs-1) {
            if(lower_neighbour[j]==0) walls_around++;
        }
    }

    // check left
    if((j-1)>=0) {
        if(A[i][j-1]==0) walls_around++;
    }

    // check right
    if((j+1)<cols) {
        if(A[i][j+1]==0) walls_around++;
    }

    return walls_around;
}

// checks whether there are remaining dead ends
bool dead_ends_present(int rows, int cols, int A[][cols]) {

    int i, j;
    for( i = 0; i < rows; i++) {
        for( j = 0; j < cols; j++ ) {
            if(A[i][j]==1 && wall_neighbours(rows,cols,A,i,j)==3) return true;
        }
    }
    return false;
}


int main(int argc, char **argv)
{   
    // initializations
    int rank, num_procs;
    int proc_id;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
    MPI_Status status;


    // master proc
    if(rank == 0) {

        // read size of the input from file
        FILE *input = fopen(argv[1], "r");
        size_t size;
        fscanf(input, "%zd", &size);

        // create the main array A
        size_t i, j;
        int A[size][size];

        // fill A from file
        for( i = 0; i < size; ++i ) {
            for( j = 0; j < size; ++j ) {
                fscanf( input, "%d", &A[i][j] );
            }
        }
        fclose(input);

        // partition array rows almost equally
        int distribution[num_procs-1];
        int start_rows[num_procs-1];
        int current_start_row = 0;
        int quotient = size / (num_procs-1);
        int remains = size % (num_procs-1);
        for(proc_id = 1; proc_id < num_procs; proc_id++) {
            start_rows[proc_id] = current_start_row;
            distribution[proc_id] = quotient;
            if(remains != 0) {
                distribution[proc_id]++;
                remains--;
            }
            current_start_row += distribution[proc_id];
        }


        // send next iteration to slave procs
        // next_itr is used as a boolean (i.e. 1 means true, 0 means false) to synchronize the loops of master and slave procs
        int next_itr;
        if(dead_ends_present(size,size,A)) next_itr = 1;
        else next_itr = 0;
        for(proc_id = 1; proc_id < num_procs; proc_id++) {
            MPI_Send( &next_itr, 1, MPI_INT, 
                            proc_id, send_data_tag, MPI_COMM_WORLD);
        }

        while(next_itr) {

            // send partitions to slave procs
            for(proc_id = 1; proc_id < num_procs; proc_id++) {
                MPI_Send( &distribution[proc_id], 1, MPI_INT, 
                            proc_id, send_data_tag, MPI_COMM_WORLD);
                MPI_Send( &size, 1, MPI_INT, 
                            proc_id, send_data_tag, MPI_COMM_WORLD);
                MPI_Send( &A[start_rows[proc_id]][0], distribution[proc_id]*size, MPI_INT, 
                            proc_id, send_data_tag, MPI_COMM_WORLD);
            }


            // recieve partitions and update A
            for(proc_id = 1; proc_id < num_procs; proc_id++) {
                MPI_Recv( &A[start_rows[proc_id]][0], distribution[proc_id]*size, MPI_INT,
                    proc_id, return_data_tag, MPI_COMM_WORLD, &status);
            }

            // send next iteration to slave procs
            // next_itr is used as a boolean (i.e. 1 means true, 0 means false) to synchronize the loops of master and slave procs
            if(dead_ends_present(size,size,A)) next_itr = 1;
            else next_itr = 0;
            for(proc_id = 1; proc_id < num_procs; proc_id++) {
                MPI_Send( &next_itr, 1, MPI_INT, 
                                proc_id, send_data_tag, MPI_COMM_WORLD);
            }
        }

        // write solved grid to output file
        FILE *output = fopen(argv[2], "w");
        for(i=0; i<size; i++) {
            for(j=0; j<size; j++) {
                fprintf(output, " %d", A[i][j]);
            }
            fprintf(output, "\n");
        }
        fclose(output);
    }




    // slave procs
    else {

        int next_itr;
        MPI_Recv( &next_itr, 1, MPI_INT,
                0, send_data_tag, MPI_COMM_WORLD, &status);


        while(next_itr) {

            // recieve partition from master proc
            int rows, cols;
            MPI_Recv( &rows, 1, MPI_INT,
                    0, send_data_tag, MPI_COMM_WORLD, &status);
            MPI_Recv( &cols, 1, MPI_INT,
                    0, send_data_tag, MPI_COMM_WORLD, &status);
            int part_array[rows][cols];
            MPI_Recv( &part_array[0][0], rows*cols, MPI_INT, 
                        0, send_data_tag, MPI_COMM_WORLD, &status);

            // create neighbour cells array for boundary conditions
            size_t i, j;
            int upper_neighbour[cols], lower_neighbour[cols];
            for(j = 0; j < cols; j++) {
                upper_neighbour[j] = 0;
                lower_neighbour[j] = 0;
            }

            // communicate slave procs to update neigbour cells array
            if(rank % 2 == 1) {
                if(rank + 1 < num_procs) {
                    MPI_Send( &part_array[rows-1][0], j, MPI_INT, 
                        rank + 1, send_data_tag, MPI_COMM_WORLD);
                    MPI_Recv( &lower_neighbour[0], j, MPI_INT, 
                        rank + 1, send_data_tag, MPI_COMM_WORLD, &status);
                }
                if(rank > 1) {
                    
                    MPI_Recv( &upper_neighbour[0], j, MPI_INT, 
                        rank - 1, send_data_tag, MPI_COMM_WORLD, &status);
                    MPI_Send( &part_array[0][0], j, MPI_INT, 
                        rank - 1, send_data_tag, MPI_COMM_WORLD);
                }
            } else {
                if(rank > 1) {
                    MPI_Recv( &upper_neighbour[0], j, MPI_INT, 
                        rank - 1, send_data_tag, MPI_COMM_WORLD, &status);
                    MPI_Send( &part_array[0][0], j, MPI_INT, 
                        rank - 1, send_data_tag, MPI_COMM_WORLD);
                }
                if(rank + 1 < num_procs) {
                    MPI_Send( &part_array[rows-1][0], j, MPI_INT, 
                        rank + 1, send_data_tag, MPI_COMM_WORLD);
                    MPI_Recv( &lower_neighbour[0], j, MPI_INT, 
                        rank + 1, send_data_tag, MPI_COMM_WORLD, &status);
                }
            }

            // mark as wall if it is a dead end
            for(i=0; i<rows; i++) {
                for(j=0; j<cols; j++) {
                    if(part_array[i][j]==1 && wall_neighbours_boundary(rows,cols,part_array,i,j
                    ,rank,num_procs,upper_neighbour,lower_neighbour)==3) {
                        part_array[i][j] = 0;
                    }
                }
            }

            // send partition to master proc
            MPI_Send( &part_array[0][0], rows*cols, MPI_INT, 
                        0, return_data_tag, MPI_COMM_WORLD);

            // recieve next iteration from master proc
            // next_itr is used as a boolean (i.e. 1 means true, 0 means false) to synchronize the loops of master and slave procs
            MPI_Recv( &next_itr, 1, MPI_INT,
                0, send_data_tag, MPI_COMM_WORLD, &status);
        }
    }

    MPI_Finalize();
    return 0;
}

