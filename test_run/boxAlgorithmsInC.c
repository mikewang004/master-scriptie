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
        double *current_row = &array[i * cols];
        
        // Check all three conditions
        if (current_row[0] == x &&    // Column 1 == x
            current_row[1] == y &&    // Column 2 == y  
            current_row[2] == z) {    // Column 3 == z
            printf("Match found \n");
            return current_row;  // Found matching row
        }
    }
    return NULL;  // No match found
}



// Implementaiton of Hoshen - Kopelman algorithm 
// Input: 2D array containing lattice id, 1D array containing crystallisation
// Output: 2D array containing lattice id + cluster id, 1d array containing cluster id sizes


void hoshen_kopelman_crystallisation(double *array, int rows, int cols, double *output_2d, double *output_1d, int *actual_1d_size, int no_boxes) {
    // Note array structure is [xbox ybox zbox cryst xev yev zev]
    // Develop method to loop over array, access xbox ybox zbox +- 1 and return xev yev zev 
    for (int i = 0; i < 3; i ++) {
        for (int j = 0; j < cols; j ++) {
            printf("%f \n", array[i + j]);
        }
        int xbox = array[i * cols];
        int ybox = array[i * cols + 1];
        int zbox = array[i * cols + 2];
        // Find six neighbors 
        printf("%i \n",((xbox-1) % no_boxes + no_boxes) % no_boxes);
        printf("%i \n",((xbox+1) % no_boxes + no_boxes) % no_boxes);
        //double *row_left = find_row_by_three_columns(array, rows, cols, ((xbox-1) % no_boxes + no_boxes) % no_boxes, ybox, zbox);
        //double *row_right = find_row_by_three_columns(array, rows, cols, ((xbox+1) % no_boxes + no_boxes) % no_boxes, ybox, zbox);
        double *same_row = find_row_by_three_columns(array, rows, cols, xbox, ybox, zbox);
        for (int j = 0; j < cols; j ++) {
            printf("%f \n", same_row[j]);
        }
        // for (int j = 0; j < cols; j ++) {
        //     printf("%f \n", row_right[j]);
        // }

        // Check for crystallinity 
    }
}







