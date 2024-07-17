#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "dumbmatrix.h"
#include "matrix.h"

int main(int argc, char* argv[]) {

	int seed = 10;
	int low = -100;
	int high = 100;
	int row = 5000;
	int col = 5000;


	printf("Let's generate a randomized array for dumb.\n");
    matrix_dumb *dumb1, *dumb2, *result_dumb;
	allocate_matrix_dumb(&(dumb1), row, col);
	allocate_matrix_dumb(&(dumb2), row, col);
	allocate_matrix_dumb(&(result_dumb), row, col);
	rand_matrix_dumb(dumb1, seed, low, high);
	rand_matrix_dumb(dumb2, seed + 1, low, high);

	
	printf("Starting randomized add_dumb.\n");
	clock_t start_dumb = clock();
	add_matrix_dumb(result_dumb, dumb1, dumb2);
	clock_t end_dumb = clock();
	printf("Time taken: %Lf s\n", (long double)(end_dumb - start_dumb) / CLOCKS_PER_SEC);
	deallocate_matrix_dumb(dumb1);
	deallocate_matrix_dumb(dumb2);


	printf("Starting randomized array for matrix.\n");
	matrix *mat1, *mat2, *result;
	allocate_matrix(&(mat1), row, col);
	allocate_matrix(&(mat2), row, col);
	allocate_matrix(&(result), row, col);

	for (int i = 0; i < row; i++) {
		for (int j = 0; j < col; j++) {
			set(mat1, i, j, get_matrix_dumb(dumb1, i, j));
			set(mat2, i, j, get_matrix_dumb(dumb2, i, j));
		}
	}

	printf("Starting randomized add_matrix.\n");
	clock_t start_mat = clock();
	add_matrix(result, mat1, mat2);
	clock_t end_mat = clock();
	printf("Time taken: %Lf s\n", (long double)(end_mat - start_mat) / CLOCKS_PER_SEC);

	printf("Test equal\n");
	for (int i = 0; i < row; i++) {
		for (int j = 0; j < col; j++) {
			if (get(result, i, j) != get_matrix_dumb(result_dumb, i, j)) {
				printf("The values of row and col do not equal: row = %d, col = %d\n, value_dumb = %f\n, value = %f\n", i, j, get_matrix_dumb(result_dumb, i, j), get(result, i, j));
			}
			
		}
		
	}
	
	
}
