// Script to find the box id of an entire box 
// Import: 1d polymer position array 
// Returns: 1d polymer box id array 

#include <stdlib.h>
#include <math.h>
#include <Python.h>




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

// double* inner_products_per_polymer(double* array, int rows, int cols) {
//     if (array == NULL || rows <= 1 || cols != 3) {
//         printf("Error: Invalid input parameters\n");
//         printf("  array: %s, rows: %d, cols: %d\n", 
//                array == NULL ? "NULL" : "valid", rows, cols);
//         return NULL;
//     }
    
//     // Allocate result array (rows-1 elements for i-th and (i+1)-th pairs)
//     double* results = (double*)malloc((rows - 1) * sizeof(double));
//     if (results == NULL) {
//         printf("Error: Memory allocation failed for %d results\n", rows - 1);
//         return NULL;
//     }
    
//     // Compute inner products between consecutive rows
//     for (int i = 0; i < rows - 1; i++) {
//         double* row_i = &array[i * cols];      // i-th row
//         double* row_i1 = &array[(i + 1) * cols]; // (i+1)-th row
        
//         // Inner product: x1*x2 + y1*y2 + z1*z2
//         results[i] = (row_i[0] * row_i1[0]) + 
//                      (row_i[1] * row_i1[1]) + 
//                      (row_i[2] * row_i1[2]);
//     }
    
//     return results;
// }

double* inner_products_per_polymer(double* array, int rows, int cols) {
    if (array == NULL || rows < 3 || cols != 3) {
        printf("Error: Invalid input parameters\n");
        printf("  array: %s, rows: %d, cols: %d\n", 
               array == NULL ? "NULL" : "valid", rows, cols);
        return NULL;
    }
    
    // Allocate result array (rows-2 elements for consecutive bond pairs)
    double* results = (double*)malloc((rows - 2) * sizeof(double));
    if (results == NULL) {
        printf("Error: Memory allocation failed for %d results\n", rows - 2);
        return NULL;
    }
    
    // Compute dot products between consecutive normalized bond vectors
    for (int i = 0; i < rows - 2; i++) {
        // Get three consecutive positions
        double* pos_i = &array[i * cols];       // Position i
        double* pos_i1 = &array[(i + 1) * cols]; // Position i+1
        double* pos_i2 = &array[(i + 2) * cols]; // Position i+2
        
        // Calculate bond vector from i to i+1
        double bond1_x = pos_i1[0] - pos_i[0];
        double bond1_y = pos_i1[1] - pos_i[1];
        double bond1_z = pos_i1[2] - pos_i[2];
        
        // Calculate bond vector from i+1 to i+2
        double bond2_x = pos_i2[0] - pos_i1[0];
        double bond2_y = pos_i2[1] - pos_i1[1];
        double bond2_z = pos_i2[2] - pos_i1[2];
        
        // Calculate magnitudes
        double mag1 = sqrt(bond1_x*bond1_x + bond1_y*bond1_y + bond1_z*bond1_z);
        double mag2 = sqrt(bond2_x*bond2_x + bond2_y*bond2_y + bond2_z*bond2_z);
        
        // Avoid division by zero for zero-length bonds
        if (mag1 < 1e-12 || mag2 < 1e-12) {
            results[i] = 0.0;
            continue;
        }
        
        // Normalize bond vectors
        double bond1_hat_x = bond1_x / mag1;
        double bond1_hat_y = bond1_y / mag1;
        double bond1_hat_z = bond1_z / mag1;
        
        double bond2_hat_x = bond2_x / mag2;
        double bond2_hat_y = bond2_y / mag2;
        double bond2_hat_z = bond2_z / mag2;
        
        // Dot product of normalized bond vectors: b̂_i • b̂_{i+1}
        results[i] = (bond1_hat_x * bond2_hat_x) + 
                     (bond1_hat_y * bond2_hat_y) + 
                     (bond1_hat_z * bond2_hat_z);
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


int check_merger(double *row, double *output_2d, int *box_array, int *cluster_id, float *ev_array, int i, int rows, int cols, float cryst_cutoff, float ndot_cutoff) {
    //int i indicates current position of output_2d, row is offset, box_array is location of current center lattice point
    if (row[3] > cryst_cutoff) { // Both lattice points should be crystalline, crystalliiny central point checked in H-K main function
        double inproduct = fabs(row[4] * ev_array[0]) + fabs(row[5] * ev_array[1]) + fabs(row[6] * ev_array[2]); 
        printf("%f * %f  + %f * %f  +%f * %f  = %f \n", row[4], ev_array[0], row[5], ev_array[1], row[6], ev_array[2], inproduct);
        //printf("%f \n", inproduct);
        if (inproduct >= ndot_cutoff) {
            //printf("Index considered \n");
            // First check if current box values are not already in output array 
            int row_index_already_in_output = find_neighbours_index(output_2d, rows, cols, row[0], row[1], row[2]);
            if (row_index_already_in_output == NULL) {
                output_2d[i] = row[0]; 
                output_2d[i + rows] = row[1];
                output_2d[i + rows * 2] = row[2];
                output_2d[i + rows * 3] = *cluster_id;
                i = i + 1;
            }
            else {
            // Also check whether item is already in cluster. If so, keep current cluster and advance counter by one 
            // TODO implement proper union algorithm
                if (output_2d[row_index_already_in_output + rows * 4] == 0) {
                    output_2d[row_index_already_in_output + rows * 3] = *cluster_id;
                }
                else {
                    *cluster_id = *cluster_id + 1;
                    printf("current cluster id %i \n", *cluster_id);
                }
            }
        }
    }
    return i;
}





// Implementaiton of Hoshen - Kopelman algorithm 
// Input: 2D array containing lattice id, 1D array containing crystallisation
// Output: 2D array containing lattice id + cluster id, 1d array containing cluster id sizes


void hoshen_kopelman_crystallisation(double *array, int rows, int cols, double *output_2d, double *output_1d, int *actual_1d_size, int no_boxes, float cryst_cutoff, float ndot_cutoff) {
    // Note array structure is [xbox ybox zbox cryst xev yev zev]
    // Develop method to loop over array, access xbox ybox zbox +- 1 and return xev yev zev 
    int last_filled_output_row = 0;
    int cluster_id = 1;
    for (int i = 0; i < rows; i ++) {
    //for (int i = 0; i < 50; i ++) {
        // Note column based indexing
        // for (int j = 0; j < cols; j ++) {
        //     printf("%f \n", array[i + j]);
        // }
        if (array[i + rows * 4] > cryst_cutoff) {


            int xbox = array[i];
            int ybox = array[i + rows];
            int zbox = array[i + rows * 2];
            float xev = array[i + rows * 4];
            float yev = array[i + rows * 5];
            float zev = array[i + rows * 6];
            int box_array[3] = {xbox, ybox, zbox};
            float ev_array[3] = {xev, yev, zev};

            // // Find six neighbors 
            // printf("%i %i %i \n", xbox,ybox,zbox);
            double *row_left = find_neighbours(array, rows, cols, ((xbox-1 + no_boxes) % no_boxes), ybox, zbox);
            double *row_right = find_neighbours(array, rows, cols, ((xbox+1) % no_boxes), ybox, zbox);
            double *row_before = find_neighbours(array, rows, cols, xbox, ((ybox-1 + no_boxes) % no_boxes), zbox);
            double *row_behind = find_neighbours(array, rows, cols, xbox, ((ybox+1) % no_boxes), zbox);
            double *row_lower = find_neighbours(array, rows, cols, xbox, ybox, ((zbox-1 + no_boxes) % no_boxes));
            double *row_upper = find_neighbours(array, rows, cols, xbox, ybox, ((zbox+1) % no_boxes));

            // // Check for crystallinity 
            int og_last_filled_output_row = last_filled_output_row;
            last_filled_output_row =  check_merger(row_left, output_2d, box_array, &cluster_id,  ev_array, last_filled_output_row, rows, cols,cryst_cutoff, ndot_cutoff);
            last_filled_output_row = check_merger(row_right, output_2d, box_array, &cluster_id, ev_array, last_filled_output_row, rows, cols, cryst_cutoff, ndot_cutoff); 
            last_filled_output_row =  check_merger(row_before, output_2d, box_array, &cluster_id,  ev_array, last_filled_output_row, rows, cols,cryst_cutoff, ndot_cutoff);
            last_filled_output_row = check_merger(row_behind, output_2d, box_array, &cluster_id, ev_array, last_filled_output_row, rows, cols,cryst_cutoff, ndot_cutoff); 
            last_filled_output_row =  check_merger(row_lower, output_2d, box_array, &cluster_id,  ev_array, last_filled_output_row, rows, cols,cryst_cutoff, ndot_cutoff);
            last_filled_output_row = check_merger(row_upper, output_2d, box_array, &cluster_id, ev_array, last_filled_output_row, rows, cols,cryst_cutoff, ndot_cutoff); 
            //printf("%i \n" ,last_filled_output_row);
            // If there is a match, also apply cluster id to center spot 
            if (og_last_filled_output_row != last_filled_output_row) {
                printf("Apply center cluster id \n");
                // Check if row is already present in output_2d
                int row_index_already_in_output = find_neighbours_index(output_2d, rows, cols, xbox, ybox, zbox);
                if (row_index_already_in_output == NULL) {
                    output_2d[last_filled_output_row] = xbox; 
                    output_2d[last_filled_output_row + rows] = ybox;
                    output_2d[last_filled_output_row + rows * 2] = zbox;
                    output_2d[last_filled_output_row + rows * 3] = cluster_id;
                }
                else {
                    output_2d[row_index_already_in_output + rows * 3] = cluster_id;
                }

                
            }
            
        }
    }
}







