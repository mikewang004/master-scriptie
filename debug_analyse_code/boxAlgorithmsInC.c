// Script to find the box id of an entire box 
// Import: 1d polymer position array 
// Returns: 1d polymer box id array 

#include <stdlib.h>
#include <math.h>
#include <Python.h>
#include <assert.h>


#define max(a,b) (a>b?a:b)






int* find_nearest_value(const double nearest_values[], size_t a_size, 
                                 const double data[], size_t n_size) {
    if (nearest_values == NULL || data == NULL || a_size == 0 || n_size == 0) {
        return NULL;
    }
    
    int* results = (int*)malloc(n_size * sizeof(int));
    if (results == NULL) return NULL;
    
    // If nearest_values is sorted, we could use binary search for O(n log a) instead of O(n*a)
    for (size_t i = 0; i < n_size; i++) {
        double current_data = data[i];
        int min_idx = 0;
        double min_diff = fabs(nearest_values[0] - current_data);
        
        // Find nearest value for current data point
        for (size_t j = 1; j < a_size; j++) {
            double diff = fabs(nearest_values[j] - current_data);
            if (diff < min_diff) {
                min_diff = diff;
                min_idx = j;
            }
        }
        results[i] = min_idx;
    }
    
    return results;
}

double inner_products_columnwise_array(double* array1, double* array2, int rows, int cols) {
    //Takes an array [rows x cols], returns inner product according to array1[1] * array2[1] + array[1][2] * array2[2] + ... + array1[n] * array2[n]
    double* results = (double*)malloc((rows) * sizeof(double));

    for (int i = 0; i < rows; i ++ ) {
        double dot_product = 0.0;
        for (int j = 0; j < cols; j ++ ) {
            int index = i + rows * j;
            dot_product += array1[index] * array2[index];
        }
        results[i] = dot_product;
    }
    return *results;
}

double* inner_products_per_polymer(double* array, int rows, int cols) {
    // Deprecated, rework later, corresponding python function is bond_bond_correlation
    if (array == NULL || rows < 3 || cols != 3) {
        printf("Error: Invalid input parameters\n");
        printf("  array: %s, rows: %d, cols: %d\n", 
               array == NULL ? "NULL" : "valid", rows, cols);
        return NULL;
    }
    
    // Allocate result array (rows-2 elements for consecutive bond pairs)
    double* results = (double*)malloc((rows - 1) * sizeof(double));
    if (results == NULL) {
        printf("Error: Memory allocation failed for %d results\n", rows - 2);
        return NULL;
    }
    
    // Compute dot products between consecutive normalized bond vectors
    for (int i = 0; i < rows - 1; i++) {
        // Get two consecutive positions
        double* pos_i = &array[i * cols];       // Position i
        double* pos_i1 = &array[(i + 1) * cols]; // Position i+1
        results[i] = (pos_i[0] * pos_i1[0]) + 
                     (pos_i[1] * pos_i1[1]) + 
                     (pos_i[2] * pos_i1[2]);
    }
    
    return results;
}


double* find_neighbours(double *array, int rows, int cols, 
                                  int x, int y, int z) {
    for (int i = 0; i < rows; i++) {
        // For column-major: array[j * rows + i] gives row i, column j
        double col1_val = array[0 * rows + i];  // First column, row i
        double col2_val = array[1 * rows + i];  // Second column, row i  
        double col3_val = array[2 * rows + i];  // Third column, row i
        
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

double* check_merger_criteria(double *ev_array, double *row, float cryst_cutoff, float ndot_cutoff) {
    if (row[3] > cryst_cutoff) { // Both lattice points should be crystalline, crystalliiny central point checked in H-K main function
        double inproduct = fabs(row[4] * ev_array[0]) + fabs(row[5] * ev_array[1]) + fabs(row[6] * ev_array[2]); 
        if (inproduct >= ndot_cutoff) { 
            // Look up label
            return 1;
        }
    }
    return 0;
}


int find_neighbours_index(double *array, int rows, int cols, 
                                  int x, int y, int z) {
    for (int i = 0; i < rows; i++) {
        // For column-major: array[j * rows + i] gives row i, column j
        double col1_val = array[0 * rows + i];  // First column, row i
        double col2_val = array[1 * rows + i];  // Second column, row i  
        double col3_val = array[2 * rows + i];  // Third column, row i
        
        // Check all three conditions
        if (col1_val == x &&    // Column 1 == x
            col2_val == y &&    // Column 2 == y  
            col3_val == z) {    // Column 3 == z
            return i;
        }
    }
    //printf("Nothing found for (%d, %d, %d)\n", x, y, z);
    return NULL;  // No match found
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

int hk_make_set(int *output_1d, int *n_labels) {
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
    int *new_labels) {
    /* Note array structure is [xbox ybox zbox cryst xev yev zev]
    // Develop method to loop over array, access xbox ybox zbox +- 1 and return xev yev zev 
    Input: cryst_array, rows, cols, cryst_cutoff, ndot_cutoff
    Output: array containing cluster-id per lattice point, labels-id array
    Implemented from https://www.ocf.berkeley.edu/~fricke/projects/hoshenkopelman/hoshenkopelman.html*/
    int cluster_id = 1;
    int steps_since_last_cluster_change = 0;
    int max_labels = nridges*nridges;
    int n_labels = max_labels;

    // Setup label array 
    int *labels = calloc(sizeof(int), max_labels);
    labels[0] = 0;

    

    int l = 0;
    for (int m = 0; m < 2; m ++ ) {
        for (int l = 0; l < nridges*nridges*nridges; l++) {
        // Note that cryst_array is sorted on y-x-z order, so that y loops most often and z loops less often. 
                    //int l = i * rows + j * rows + k;
            if (array[l + rows * 4] > cryst_cutoff) {
                //printf("%f \n", array[l]);
                int xbox = array[l];
                int ybox = array[l + rows];
                int zbox = array[l + rows * 2];
                double xev = array[l + rows * 4];
                double yev = array[l + rows * 5];
                double zev = array[l + rows * 6];
                int box_array[3] = {xbox, ybox, zbox};
                double ev_array[3] = {xev, yev, zev};

                double *row_left = find_neighbours(array, rows, cols, ((xbox-1 + nridges) % nridges), ybox, zbox);
                double *row_before = find_neighbours(array, rows, cols, xbox, ((ybox-1 + nridges) % nridges), zbox);
                double *row_upper = find_neighbours(array, rows, cols, xbox, ybox, ((zbox-1 + nridges) % nridges));

                int i = xbox; int j = ybox; int k = zbox;
                //printf("%i %i %i", i,j,k);

                // Confirm if rows confirm to crystallinity and eigenvalue criteria

                int row_left_check = check_merger_criteria(ev_array, row_left, cryst_cutoff, ndot_cutoff);
                int row_before_check = check_merger_criteria(ev_array, row_before, cryst_cutoff, ndot_cutoff);
                int row_upper_check = check_merger_criteria(ev_array, row_upper, cryst_cutoff, ndot_cutoff);
        

                // Set matrix positions to 1 if they do not have a cluster label 
                label_matrix[(i-1 + nridges) % nridges][j][k] = hk_set_position_1(row_left_check,  label_matrix[(i-1 + nridges) % nridges][j][k]);
                label_matrix[i][(j-1 + nridges) % nridges][k] = hk_set_position_1(row_before_check,  label_matrix[i][(j-1 + nridges) % nridges][k]);
                label_matrix[i][j][(k-1 + nridges) %nridges] =  hk_set_position_1(row_upper_check,  label_matrix[i][j][(k-1 + nridges) %nridges]);
                

                // Also save actual matrix labels 
                int position_left =  label_matrix[(i-1 + nridges) % nridges][j][k];
                int position_before = label_matrix[i][(j-1 + nridges) % nridges][k];
                int position_upper = label_matrix[i][j][(k-1 + nridges) %nridges];

                //printf("check left %i actual position %i \n", row_left_check,position_left);
                



                switch (!! row_left_check + !!row_before_check + !!row_upper_check) {
                    case 0: //Make new cluster 
                        if (m == 0) {
                            label_matrix[i][j][k] = hk_make_set(labels, n_labels);
                            //printf("at position %i %i %i, counter %i new cluster %i \n", i,j,k,l,label_matrix[i][j][k]);
                        }

                        break;
                    case 1: //Add to cluster left/above
                        if (m == 0){
                            label_matrix[i][j][k] = max(position_before,max(position_upper,position_left));  
                            //printf("at position %i %i %i, added to cluster %i \n", i,j,k,label_matrix[i][j][k]);
                        }
                        else {
                            if (position_before) {
                                label_matrix[i][j][k] = position_before;
                            }
                            if (position_upper) {
                                label_matrix[i][j][k] = position_upper;
                            }
                            if (position_left) {
                                label_matrix[i][j][k] = position_left;
                            }
                        }
                        break;
                    case 2: //Combine two clusters 
                        //printf("case 2 \n");
                        // label_matrix[i][j][k] = (!! position_before == 0 ? hk_union(labels, position_left, position_upper) : (!! position_upper == 0 ?
                        //     hk_union(labels, position_before, position_left) : hk_union(labels, position_before, position_upper)));
                        //printf("merged two clusters, current position %i %i %i, new cluster %i\n",i,j,k,label_matrix[i][j][k] );

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
                        break;
                    case 3: //Combine three clusters
                        label_matrix[i][j][k] = hk_union(labels, position_left, position_upper);
                        label_matrix[i][j][k] = hk_union(labels, position_before, position_left);
                        label_matrix[i][j][k] = hk_union(labels, position_before, position_upper);
                        break;
                }
                //printf("label matrix %i\n", label_matrix[i][j][k]);
            }
        }
    //printf("first loop done \n");
    }

    // for (int i = 0; i < 20; i ++) {
    //     printf("labels[%i] = %i \n", i, labels[i]);
    //     //new_labels[i] = labels[i];
    // }
 
    // for (int i = 0; i < nridges; i ++) {
    //     for (int j = 0; j < nridges; j ++) {
    //         for (int k = 0; k < nridges; k ++) {
    //             //printf("label matrix value %d \n", label_matrix[i][j][k]);
                
    // }
    // }
    // }


    // Do new label binding 
    // TODO fix recursive labelling

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





}







