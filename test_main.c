#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "dumbmatrix.h"
#include "matrix.h"

double round(double x) {
    if (x >= 0) {
        return (int)(x + 0.5);
    } else {
        return (int)(x - 0.5);
    }
}

double my_pow(double base, int exponent) {
    double result = 1.0;
    int exp = exponent > 0 ? exponent : -exponent; // 处理负指数

    for (int i = 0; i < exp; i++) {
        result *= base;
    }

    // 如果指数是负数，则取倒数
    if (exponent < 0) {
        return 1.0 / result;
    }
    return result;
}

int main(int argc, char* argv[]) {

	int seed = 10;
	int low = -100;
	int high = 100;
	int row = 1000;
	int col = 1000;
	int precision = 4;
	float factor = my_pow(10.0, precision);

	clock_t start_dumb, end_dumb, delta_dumb, start_mat, end_mat, delta_mat;
	printf("Let's generate a randomized array for dumb.\n");
    matrix_dumb *dumb1, *dumb2, *result_dumb;
	allocate_matrix_dumb(&(dumb1), row, col);
	allocate_matrix_dumb(&(dumb2), row, col);
	allocate_matrix_dumb(&(result_dumb), row, col);
	rand_matrix_dumb(dumb1, seed, low, high);
	rand_matrix_dumb(dumb2, seed + 1, low, high);

	
	printf("Starting randomized mul_dumb.\n");
	start_dumb = clock();
	//mul_matrix_dumb(result_dumb, dumb1, dumb2);
	end_dumb = clock();
	delta_dumb = end_dumb - start_dumb;
	printf("Time taken: %Lf s\n", (long double)(delta_dumb) / CLOCKS_PER_SEC);


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

	printf("Starting randomized mul_matrix.\n");
	start_mat = clock();
	//mul_matrix(result, mat1, mat2);
	end_mat = clock();
	delta_mat = end_mat - start_mat;
	printf("Time taken: %Lf s\n", (long double)(delta_mat) / CLOCKS_PER_SEC);

	printf("Test mul equal\n");

	for (int i = 0; i < row; i++) {
		for (int j = 0; j < col; j++) {
			double dumb = get_matrix_dumb(result_dumb, i, j);
			double num = get(result, i, j);
			dumb = round(dumb * factor) / factor;
			num = round(num * factor) / factor;
			if (dumb != num) {
				printf("The values of row and col do not equal: row = %d, col = %d, value_dumb = %f, value = %f\n", i, j, dumb, num);
			}
			
		}
	}
	printf("Mul Speedup is %f\n",  (double)(delta_dumb) / (double)(delta_mat));


	printf("Starting randomized abs_dumb.\n");
	start_dumb = clock();
	abs_matrix_dumb(result_dumb, dumb1);
	end_dumb = clock();
	delta_dumb = end_dumb - start_dumb;
	printf("Time taken: %Lf s\n", (long double)(delta_dumb) / CLOCKS_PER_SEC);

	printf("Starting randomized abs_matrix.\n");
	start_mat = clock();
	abs_matrix(result, mat1);
	end_mat = clock();
	delta_mat = end_mat - start_mat;
	printf("Time taken: %Lf s\n", (long double)(delta_mat) / CLOCKS_PER_SEC);

	printf("Test abs equal\n");

	for (int i = 0; i < row; i++) {
		for (int j = 0; j < col; j++) {
			double dumb = get_matrix_dumb(result_dumb, i, j);
			double num = get(result, i, j);
			dumb = round(dumb * factor) / factor;
			num = round(num * factor) / factor;
			if (dumb != num) {
				printf("The values of row and col do not equal: row = %d, col = %d, value_dumb = %f, value = %f\n", i, j, dumb, num);
			}
			
		}
	}
	printf("Abs Speedup is %f\n",  (double)(delta_dumb) / (double)(delta_mat));


	printf("Starting randomized add_dumb.\n");
	start_dumb = clock();
	add_matrix_dumb(result_dumb, dumb1, dumb2);
	end_dumb = clock();
	delta_dumb = end_dumb - start_dumb;
	printf("Time taken: %Lf s\n", (long double)(delta_dumb) / CLOCKS_PER_SEC);

	printf("Starting randomized add_matrix.\n");
	start_mat = clock();
	add_matrix(result, mat1, mat2);
	end_mat = clock();
	delta_mat = end_mat - start_mat;
	printf("Time taken: %Lf s\n", (long double)(delta_mat) / CLOCKS_PER_SEC);

	printf("Test add equal\n");

	for (int i = 0; i < row; i++) {
		for (int j = 0; j < col; j++) {
			double dumb = get_matrix_dumb(result_dumb, i, j);
			double num = get(result, i, j);
			dumb = round(dumb * factor) / factor;
			num = round(num * factor) / factor;
			if (dumb != num) {
				printf("The values of row and col do not equal: row = %d, col = %d, value_dumb = %f, value = %f\n", i, j, dumb, num);
			}
			
		}
	}
	printf("Add Speedup is %f\n",  (double)(delta_dumb) / (double)(delta_mat));
}


