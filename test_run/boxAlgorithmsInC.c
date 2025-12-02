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
        // Get three consecutive positions
        double* pos_i = &array[i * cols];       // Position i
        double* pos_i1 = &array[(i + 1) * cols]; // Position i+1
        // double* pos_i2 = &array[(i + 2) * cols]; // Position i+2
        
        // Calculate bond vector from i to i+1
        // double bond1_x = pos_i1[0] - pos_i[0];
        // double bond1_y = pos_i1[1] - pos_i[1];
        // double bond1_z = pos_i1[2] - pos_i[2];
        
        // // Calculate bond vector from i+1 to i+2
        // double bond2_x = pos_i2[0] - pos_i1[0];
        // double bond2_y = pos_i2[1] - pos_i1[1];
        // double bond2_z = pos_i2[2] - pos_i1[2];
        
        // // Calculate magnitudes
        // double mag1 = sqrt(bond1_x*bond1_x + bond1_y*bond1_y + bond1_z*bond1_z);
        // double mag2 = sqrt(bond2_x*bond2_x + bond2_y*bond2_y + bond2_z*bond2_z);
        
        // // Avoid division by zero for zero-length bonds
        // if (mag1 < 1e-12 || mag2 < 1e-12) {
        //     results[i] = 0.0;
        //     continue;
        // }
        
        // // Normalize bond vectors
        // double bond1_hat_x = bond1_x / mag1;
        // double bond1_hat_y = bond1_y / mag1;
        // double bond1_hat_z = bond1_z / mag1;
        
        // double bond2_hat_x = bond2_x / mag2;
        // double bond2_hat_y = bond2_y / mag2;
        // double bond2_hat_z = bond2_z / mag2;
        
        // Dot product of normalized bond vectors: b̂_i • b̂_{i+1}
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




int hk_find(int *output_1d, int cluster_id) {
    int y = cluster_id;
    while (output_1d[y] != y) {
        y = output_1d[y];
    }
    while (output_1d[cluster_id] != cluster_id) {
        int z = output_1d[cluster_id];
        output_1d[cluster_id] = y;
        cluster_id = z;
    }

    return y;
}

void hk_union(int *output_1d, int x, int y) {
    output_1d[hk_find(output_1d, x)] = hk_find(output_1d, y);
}

int hk_make_set(int *output_1d, int *n_labels) {
    output_1d[0]++;
    assert(output_1d[0] <= n_labels);
    output_1d[output_1d[0]] = output_1d[0];
    return output_1d[0];
}



// Implementaiton of Hoshen - Kopelman algorithm 
// Input: 2D array containing lattice id, 1D array containing crystallisation
// Output: 2D array containing lattice id + cluster id


void hoshen_kopelman_crystallisation(double *array, int rows, int cols, int nridges, float cryst_cutoff, float ndot_cutoff) {
    /* Note array structure is [xbox ybox zbox cryst xev yev zev]
    // Develop method to loop over array, access xbox ybox zbox +- 1 and return xev yev zev 
    Input: cryst_array, rows, cols, cryst_cutoff, ndot_cutoff
    Output: array containing cluster-id per lattice point, labels-id array
    Implemented from https://www.ocf.berkeley.edu/~fricke/projects/hoshenkopelman/hoshenkopelman.html*/
    int cluster_id = 1;
    int steps_since_last_cluster_change = 0;
    int max_labels = nridges*nridges*nridges;
    int n_labels = max_labels;
    int init_matrix_counter = 0;

    // Setup label array 
    int* labels = calloc(sizeof(int), max_labels);
    labels[0] = 0;

    //Initialise matrix with cluster labels 
    int ***label_matrix; 
    	if (nridges)  
	{
		//memory management - allocate memory to matrix data structure
		label_matrix = (int ***)calloc(nridges, sizeof(int**));
		for (int i=0; i<nridges; i++)
		{
			label_matrix[i] = (int **)calloc(nridges, sizeof(int*));
			for (int j=0; j<nridges; j++)
			{
				label_matrix[i][j] = (int *)calloc(nridges,sizeof(int));
				for (int k=0; k<nridges; k++)
				{
                    label_matrix[i][j][k] = init_matrix_counter;
                    init_matrix_counter = init_matrix_counter + 1;
                    printf("%i \n", label_matrix[i][j][k]);
				}
			}
		}
    }


    for (int i = 0; i < nridges; i ++) {
        for (int j = 0; j < nridges; j ++) {
            for (int k = 0; k < nridges; k ++) {
    // Note k corresponds to y-dimension, j to x-dimension, i to z-dimension
                int l = i * rows + j * rows + k;
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
                    double *row_upper = find_neighbours(array, rows, cols, xbox, ybox, ((zbox+1) % nridges));

                    // Confirm if rows confirm to crystallinity and eigenvalue criteria

                    int row_left_check = check_merger_criteria(ev_array, row_left, cryst_cutoff, ndot_cutoff);
                    int row_before_check = check_merger_criteria(ev_array, row_before, cryst_cutoff, ndot_cutoff);
                    int row_upper_check = check_merger_criteria(ev_array, row_upper, cryst_cutoff, ndot_cutoff);
                }
                // TODO implement different cases what to do if rows confirm
            }
        }
    }

}







