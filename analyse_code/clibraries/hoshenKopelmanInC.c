#include <stdlib.h>
#include <math.h>
#include <Python.h>
#include <assert.h>
#include <stdbool.h>


#define max(a,b) (a>b?a:b)


void print_double_array(const double *arr, int size) {
    if (arr == NULL || size <= 0) {
        printf("[]\n");
        return;
    }

    printf("[");
    for (int i = 0; i < size; i++) {
        printf("%g", arr[i]);  // %g prints doubles nicely (no unnecessary zeros)
        if (i < size - 1) {
            printf(", ");
        }
    }
    printf("]\n");
}


int cmp_int(const void *a, const void *b) {
    int ia = *(const int *)a;
    int ib = *(const int *)b;
    return (ia > ib) - (ia < ib);  // returns -1, 0, or 1
}

void count_frequencies(const int *array, int n) {
    if (n <= 0) return;

    // Make a copy so we don't modify the original array
    int *copy = malloc(n * sizeof *copy);
    if (!copy) {
        perror("malloc");
        return;
    }
    memcpy(copy, array, n * sizeof *copy);

    // Sort the copy
    qsort(copy, n, sizeof *copy, cmp_int);

    // Count occurrences of each distinct value
    int current = copy[0];
    int count   = 1;
    int total_count = 0;

    for (int i = 1; i < n; i++) {
        if (copy[i] == current) {
            count++;
        } else {
            printf("%d occurs %d times\n", current, count);
            current = copy[i];
            count   = 1;
            total_count += count;
        }
    }
    // Print the last value

    printf("%d occurs %d times\n", current, count);
    total_count += count;
    printf("total count is %d \n", total_count);
    free(copy);

}


double* find_neighbours(double *array, int rows, int cols, 
                                  int x, int y, int z) {
    for (int i = 0; i < rows; i++) {
        // For column-major: array[j * rows + i] gives row i, column j
        int col1_val = array[0 * rows + i];  // First column, row i
        int col2_val = array[1 * rows + i];  // Second column, row i  
        int col3_val = array[2 * rows + i];  // Third column, row i
        
        // Check all three conditions
        if (col1_val == x &&    // Column 1 == x
            col2_val == y &&    // Column 2 == y  
            col3_val == z) {    // Column 3 == z
            //printf("Match found at row %d\n", i);
            
            // Return pointer to the beginning of this row
            // In column-major, we need to create a row pointer
            double *row_ptr = (double *)malloc(cols * sizeof(double));
            for (int j = 0; j < cols; j++) {
                row_ptr[j] = array[j * rows + i];  // Column j, row i
            }
            return row_ptr;
        }
    }
    //printf("Nothing found for (%d, %d, %d)\n", x, y, z);
    return NULL;  // No match found
}

int check_merger_criteria(double *ev_array, double *row, float cryst_cutoff, double ndot_cutoff, bool skip_row) {
    if (skip_row) {
        return 0;
    }
    if (row[3] > cryst_cutoff) { // Both lattice points should be crystalline, crystalliiny central point checked in H-K main function
        double inproduct = fabs(row[4])*fabs(ev_array[0]) + fabs(row[5]) * fabs(ev_array[1]) + fabs(row[6]) * fabs(ev_array[2]);
        //double inproduct = fabs(row[4] * ev_array[0]) + fabs(row[5] * ev_array[1]) + fabs(row[6] * ev_array[2]); 
        //double inproduct = fabs(row[4]*ev_array[0] + row[5] * ev_array[1] + row[6] * ev_array[2]);
        if (inproduct >= ndot_cutoff) { 
            // Look up label
            return 1;
        }
    }
    return 0;
}

int hk_find(int *output_1d, int x) {
    // Finds where x is in output_1d
    int y = x;
    while (output_1d[y] != y) {
        y = output_1d[y];
    }
    while (output_1d[x] != x) {
        int z = output_1d[x];
        output_1d[x] = y;
        x = z;
    }

    return y;
}

int hk_union(int *output_1d, int x, int y) {
    //printf("cluster of position x = %i, cluster of position y = %i\n",x,y);
    //printf("find result for x is %i, for y is %i \n", hk_find(output_1d, x), hk_find(output_1d, y));
    return output_1d[hk_find(output_1d, x)] = hk_find(output_1d, y);
}

int hk_make_set(int *output_1d, int n_labels) {
    output_1d[0]++;
    assert(output_1d[0] <= n_labels);
    output_1d[output_1d[0]] = output_1d[0];
    return output_1d[0];
}

int hk_set_position_1(int check, int *label_matrix_position) {
    if (label_matrix_position == 0)  {
        return check;
    }
    return label_matrix_position;
}


// Implementaiton of Hoshen - Kopelman algorithm 
// Input: 2D array containing lattice id, 1D array containing crystallisation
// Output: 2D array containing lattice id + cluster id


void hoshen_kopelman_crystallisation(double *array, int rows, int cols, int nridges, float cryst_cutoff, float ndot_cutoff, int ***label_matrix, int ***new_label_matrix, 
    int *labels) {
    /* Note array structure is [xbox ybox zbox cryst xev yev zev]
    // Develop method to loop over array, access xbox ybox zbox +- 1 and return xev yev zev 
    Input: cryst_array, rows, cols, cryst_cutoff, ndot_cutoff
    Output: array containing cluster-id per lattice point, labels-id array
    Implemented from https://www.ocf.berkeley.edu/~fricke/projects/hoshenkopelman/hoshenkopelman.html*/
    int max_labels = nridges*nridges*nridges;
    int n_labels = max_labels;

    // Setup label array 
    int *new_labels = calloc(sizeof(int), max_labels);
    new_labels[0] = 0;

    

    for (int m = 0; m < 1; m ++ ) {
        //for (int l = 0; l < nridges*nridges*nridges; l++) {
        for (int l = 0; l < nridges*nridges *nridges; l++) {


        //for (int l = 0; l < nridges * 2; l ++)
        // Note that cryst_array is sorted on y-x-z order, so that y loops most often and z loops less often. 
                    //int l = i * rows + j * rows + k;
            if (array[l + rows * 3] > cryst_cutoff) {
                //printf("%i \n", l);
                //printf("%f \n", array[l]);
                int xbox = array[l];
                int ybox = array[l + rows];
                int zbox = array[l + rows * 2];
                //double cryst = array[l + rows * 3];
                double xev = array[l + rows * 4];
                double yev = array[l + rows * 5];
                double zev = array[l + rows * 6];
                //int box_array[3] = {xbox, ybox, zbox};
                double ev_array[3] = {xev, yev, zev};

                int x_neighbour = xbox-1;
                int y_neighbour = ybox-1;
                int z_neighbour = zbox - 1;

                bool skiprowx = false;
                bool skiprowy = false;
                bool skiprowz = false;
                if (m == 0){
                    if (xbox == 0){
                        skiprowx = true;
                        x_neighbour = 0;
                    }
                    if (ybox == 0) {
                        skiprowy = true;
                        y_neighbour = 0;
                    }
                    if (zbox == 0) {
                        skiprowz = true;
                        z_neighbour = 0;
                    }
                }
                if (m > 0) {
                    x_neighbour = (xbox-1 + nridges) % nridges;
                    y_neighbour = (ybox-1 + nridges) % nridges;
                    z_neighbour = (zbox - 1 + nridges) %nridges;
                }


                double *row_left = find_neighbours(array, rows, cols, (x_neighbour), ybox, zbox);
                double *row_before = find_neighbours(array, rows, cols, xbox, (y_neighbour), zbox);
                double *row_upper = find_neighbours(array, rows, cols, xbox, ybox, (z_neighbour));

                //double *row_right = find_neighbours(array, rows, cols, ((xbox+1 + nridges) % nridges), ybox, zbox);
                //double *row_after = find_neighbours(array, rows, cols, xbox, ((ybox+1 + nridges) % nridges), zbox);
                //double *row_under = find_neighbours(array, rows, cols, xbox, ybox, ((zbox+1 + nridges) % nridges));
                //printf("\n for item %i %i %i %f %f %f %f \n", xbox, ybox, zbox, cryst, xev, yev, zev);
                int i = xbox; int j = ybox; int k = zbox;
                //printf("%i %i %i", i,j,k);

                // Confirm if rows confirm to crystallinity and eigenvalue criteria

                int row_left_check = check_merger_criteria(ev_array, row_left, cryst_cutoff, ndot_cutoff, skiprowx);
                int row_before_check = check_merger_criteria(ev_array, row_before, cryst_cutoff, ndot_cutoff, skiprowy);
                int row_upper_check = check_merger_criteria(ev_array, row_upper, cryst_cutoff, ndot_cutoff, skiprowz);

                //int row_right_check = check_merger_criteria(ev_array, row_right, cryst_cutoff, ndot_cutoff);
                //int row_after_check = check_merger_criteria(ev_array, row_after, cryst_cutoff, ndot_cutoff);
                //int row_under_check = check_merger_criteria(ev_array, row_under, cryst_cutoff, ndot_cutoff);
                // printf("found for -x, -y, -z neighbour respectively\n");
                // print_double_array(row_left, 7);
                // print_double_array(row_before, 7);
                // print_double_array(row_upper, 7);

                // printf("found for +x, +y, +z neighbour respectively\n");
                // print_double_array(row_right, 7);
                // print_double_array(row_after, 7);
                // print_double_array(row_under, 7);
        

                // Set matrix positions to 1 if they do not have a cluster label 
                //label_matrix[(i-1 + nridges) % nridges][j][k] = hk_set_position_1(row_left_check,  label_matrix[(i-1 + nridges) % nridges][j][k]);
                //label_matrix[i][(j-1 + nridges) % nridges][k] = hk_set_position_1(row_before_check,  label_matrix[i][(j-1 + nridges) % nridges][k]);
                //label_matrix[i][j][(k-1 + nridges) %nridges] =  hk_set_position_1(row_upper_check,  label_matrix[i][j][(k-1 + nridges) %nridges]);
                

                // Also save actual matrix labels 
                int position_left =  label_matrix[x_neighbour][j][k];
                int position_before = label_matrix[i][y_neighbour][k];
                int position_upper = label_matrix[i][j][z_neighbour];


                //printf("check left %i actual position %i \n", row_left_check,position_left);
                // Dbug interface supplied below
                // if (xbox == 6) {
                //     if (ybox == 12) {
                //         if (zbox == 23) {
                //             printf("\n for item %i %i %i %f %f %f %f \n", xbox, ybox, zbox, cryst, xev, yev, zev);
                //             printf("found for -x, -y, -z neighbour respectively\n");
                //             print_double_array(row_left, 7);
                //             print_double_array(row_before, 7);
                //             print_double_array(row_upper, 7);

                //             printf("found for +x, +y, +z neighbour respectively\n");
                //             print_double_array(row_right, 7);
                //             print_double_array(row_after, 7);
                //             print_double_array(row_under, 7);
                //             printf("case %i selected \n", !! row_left_check + !!row_before_check + !!row_upper_check);
                //             printf("%i %i %i \n", row_left_check, row_before_check, row_upper_check);
                //             printf("%i %i %i \n", row_right_check, row_after_check, row_under_check);
                //             printf("%i %i %i \n", position_left, position_before, position_upper);
                //         }
                //     }
                // }


                switch (!! row_left_check + !!row_before_check + !!row_upper_check) {
                    case 0: //Make new cluster 
                        if (m == 0) {
                            //if (row_after_check + row_right_check + row_under_check > 0) {
                                //label_matrix[i][j][k] = hk_make_set(labels, n_labels);
                            //}
                            label_matrix[i][j][k] = hk_make_set(labels, n_labels);
                            
                            // printf("case 0, m = %i, at position %i %i %i, counter %i new cluster %i \n", m,i,j,k,l,label_matrix[i][j][k]);
                        }

                        break;
                    case 1: //Add to cluster left/above
                            // printf("row check %i %i %i \n", row_left_check, row_before_check, row_upper_check);
                            // printf("%i %i %i \n", position_left, position_before, position_upper);
                            //printf("case 1, at position %i %i %i, counter %i, added to cluster %i \n", i,j,k,l,max(position_before,max(position_upper,position_left)));
                            label_matrix[i][j][k] = max(position_before,max(position_upper,position_left));  
                        if (row_before_check) {
                            label_matrix[i][j][k] = position_before;
                        }
                        if (row_upper_check) {
                            label_matrix[i][j][k] = position_upper;
                        }
                        if (row_left_check) {
                            label_matrix[i][j][k] = position_left;
                        }
                        break;
                    case 2: //Combine two clusters 
                        //printf("case 2 \n");
                        // label_matrix[i][j][k] = (!! position_before == 0 ? hk_union(labels, position_left, position_upper) : (!! position_upper == 0 ?
                        //     hk_union(labels, position_before, position_left) : hk_union(labels, position_before, position_upper)));

                        if (!! position_before == 0) {
                            label_matrix[i][j][k] = hk_union(labels, position_left, position_upper);
                        }
                        else {
                            if (!! position_upper == 0) {
                                label_matrix[i][j][k] = hk_union(labels, position_left, position_before);
                            }
                            else {
                                label_matrix[i][j][k] = hk_union(labels, position_before, position_upper);
                            }
                        }
                        //printf("merged two clusters, counter %i, current position %i %i %i, new cluster %i\n",l,i,j,k,label_matrix[i][j][k]);
                        break;
                    case 3: //Combine three clusters
                    //printf("case 3 \n");
                        label_matrix[i][j][k] = hk_union(labels, position_left, position_upper);
                        label_matrix[i][j][k] = hk_union(labels, position_before, position_left);
                        label_matrix[i][j][k] = hk_union(labels, position_before, position_upper);
                        break;
                        //printf("merged two clusters, counter %i, current position %i %i %i, new cluster %i\n",l,i,j,k,label_matrix[i][j][k]);
                }
                //printf("label matrix %i\n", label_matrix[i][j][k]);

                // printf("hk assigned %i %i %i to label %i\n",i,j,k, label_matrix[i][j][k]);
                // printf("ev = {%f %f %f} \n", xev, yev, zev);
                // printf("neighbours:-x %i, -y %i, -z %i \n \n", position_left, position_before, position_upper);
            }
        }
    //printf("first loop done \n");
    }

    // for (int i = 0; i < nridges; i ++) {
    //     for (int j = 0; j < nridges; j ++) {
    //         for (int k = 0; k < nridges; k ++) {
    //             if (label_matrix[i][j][k] != 0) {
    //                 printf("label_matrix %i %i %i  = %d \n", i,j,k, label_matrix[i][j][k]);           
    //             }

                
    // }
    // }
    // }

    // Print labels 

    // for (int i = 0; i < 200; i ++){
    //     printf("labels[%i] %i \n", i, labels[i]);
    // }

    // printf("Frequencies as follows \n");
    // count_frequencies(labels, max_labels);


    // Do new label binding 
    // TODO fix recursive labelling
    //printf("%i \n", label_matrix[6][11][23]);

    for (int i = 0; i < nridges; i ++) {
        for (int j = 0; j < nridges; j ++) {
            for (int k = 0; k < nridges; k ++) {
                if (label_matrix[i][j][k] > 0) {
                    int x = hk_find(labels, label_matrix[i][j][k]);
                    //printf("x = %i for label_matrix[%i][%i][%i] = %i \n", x,i,j,k,label_matrix[i][j][k]);
                    if (new_labels[x] == 0) {
                        new_labels[0]++;
                        new_labels[x] = new_labels[0];
                    }
                    label_matrix[i][j][k] = new_labels[x];
                }
            }
        }
    } 





    // printf("New labeling as follows: \n\n\n");
    // for (int i = 0; i < nridges; i ++) {
    //     for (int j = 0; j < nridges; j ++) {
    //         for (int k = 0; k < nridges; k ++) {
    //             if (label_matrix[i][j][k] != 0) {
    //                 printf("new_label_matrix %i %i %i  = %d \n", i,j,k, label_matrix[i][j][k]);           
    //             }

                
    // }
    // }
    // }

    // for (int i = 0; i < 200; i ++){

    //     printf("new_labels[%i] %i \n", i, new_labels[i]);
    // }

    // printf("New labeling frequencies as follows \n");
    // count_frequencies(new_labels, max_labels);

}