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

void bond_bond_correlation(double* bond_array, double *bond_bond_triag_arr, int no_bonds, int dims) {
    // NB no_bonds is amount of bond vectors in a given polymer. This is the same as the amount of monomers in a polymer - 1
    // no_bonds is also the amount of rows
    for (int i = 0; i < no_bonds; i ++) {
        double pos_i0 = bond_array[i];
        double pos_i1 = bond_array[i + no_bonds];
        double pos_i2 = bond_array[i + (2 * no_bonds)];
        // printf("%f %f %f\n", pos_i0, pos_i1, pos_i2);
        for (int j = i; j < no_bonds; j ++) {
            double pos_j0 = bond_array[j];
            double pos_j1 = bond_array[j + no_bonds];
            double pos_j2 = bond_array[j + (2 * no_bonds)];
            bond_bond_triag_arr[i * no_bonds+ j ] = (pos_i0 *pos_j0) + 
                    (pos_i1 * pos_j1) + 
                    (pos_i2 * pos_j2);
        }
    }
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














