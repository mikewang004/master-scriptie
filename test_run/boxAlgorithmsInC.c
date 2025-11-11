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

double* find_row_by_three_columns(double *array, int rows, int cols, 
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
            printf("Match found at row %d\n", i);
            
            // Return pointer to the beginning of this row
            // In column-major, we need to create a row pointer
            double *row_ptr = (double *)malloc(cols * sizeof(double));
            for (int j = 0; j < cols; j++) {
                row_ptr[j] = array[j * rows + i];  // Column j, row i
            }
            return row_ptr;
        }
    }
    printf("Nothing found for (%d, %d, %d)\n", x, y, z);
    return NULL;  // No match found
}





// Implementaiton of Hoshen - Kopelman algorithm 
// Input: 2D array containing lattice id, 1D array containing crystallisation
// Output: 2D array containing lattice id + cluster id, 1d array containing cluster id sizes


void hoshen_kopelman_crystallisation(double *array, int rows, int cols, double *output_2d, double *output_1d, int *actual_1d_size, int no_boxes) {
    // Note array structure is [xbox ybox zbox cryst xev yev zev]
    // Develop method to loop over array, access xbox ybox zbox +- 1 and return xev yev zev 
    for (int i = 0; i < 3; i ++) {
        // Note column based indexing
        // for (int j = 0; j < cols; j ++) {
        //     printf("%f \n", array[i + j]);
        // }
        int xbox = array[i];
        int ybox = array[i + rows];
        int zbox = array[i + rows * 2];
        //printf("array[%d] = %f\n", i+150000, array[i+150000]);

        // // Find six neighbors 
        // //printf("%i \n",((xbox-2 + no_boxes) % no_boxes));
        printf("%i %i %i \n", (xbox-1 + no_boxes) % no_boxes,ybox,zbox);
        // //printf("%i \n",((xbox+1) % no_boxes + no_boxes) % no_boxes);
        double *row_left = find_row_by_three_columns(array, rows, cols, ((xbox-1 + no_boxes) % no_boxes), ybox, zbox);
        double *row_right = find_row_by_three_columns(array, rows, cols, ((xbox+1) % no_boxes), ybox, zbox);
        double *row_before = find_row_by_three_columns(array, rows, cols, xbox, ((ybox-1 + no_boxes) % no_boxes), zbox);
        double *row_behind = find_row_by_three_columns(array, rows, cols, xbox, ((ybox+1) % no_boxes), zbox);
        double *row_lower = find_row_by_three_columns(array, rows, cols, xbox, ybox, ((zbox-1 + no_boxes) % no_boxes));
        double *row_upper = find_row_by_three_columns(array, rows, cols, xbox, ybox, ((zbox+1) % no_boxes));
        // //double *same_row = find_row_by_three_columns(array, rows, cols, xbox, ybox, zbox);
        for (int k = 0; k < cols; k ++) {
            printf("%f \n", row_left[k]);
            //printf("%f row array elements \n", row_left[k]);
            
        }
        for (int j = 0; j < cols; j ++) {
            printf("%f \n", row_right[j]);
        }

        // // Check for crystallinity 
    }
}







