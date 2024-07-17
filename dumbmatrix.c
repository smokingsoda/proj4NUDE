#include "dumbmatrix.h"
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <omp.h>

// Include SSE intrinsics
#if defined(_MSC_VER)
#include <intrin.h>
#elif defined(__GNUC__) && (defined(__x86_64__) || defined(__i386__))
#include <immintrin.h>
#include <x86intrin.h>
#endif

/* Below are some intel intrinsics that might be useful
 * void _mm256_storeu_pd (double * mem_addr, __m256d a)
 * __m256d _mm256_set1_pd (double a)
 * __m256d _mm256_set_pd (double e3, double e2, double e1, double e0)
 * __m256d _mm256_loadu_pd (double const * mem_addr)
 * __m256d _mm256_add_pd (__m256d a, __m256d b)
 * __m256d _mm256_sub_pd (__m256d a, __m256d b)
 * __m256d _mm256_fmadd_pd (__m256d a, __m256d b, __m256d c)
 * __m256d _mm256_mul_pd (__m256d a, __m256d b)
 * __m256d _mm256_cmp_pd (__m256d a, __m256d b, const int imm8)
 * __m256d _mm256_and_pd (__m256d a, __m256d b)
 * __m256d _mm256_max_pd (__m256d a, __m256d b)
 */

/* Generates a random double between low and high */
double rand_double_dumb(double low, double high)
{
    double range = (high - low);
    double div = RAND_MAX / range;
    return low + (rand() / div);
}

/* Generates a random matrix_dumb */
void rand_matrix_dumb(matrix_dumb *result, unsigned int seed, double low, double high)
{
    srand(seed);
    for (int i = 0; i < result->rows; i++)
    {
        for (int j = 0; j < result->cols; j++)
        {
            set_matrix_dumb(result, i, j, rand_double_dumb(low, high));
        }
    }
}

/*
 * Returns the double value of the matrix_dumb_dumb at the given row and column.
 * You may assume `row` and `col` are valid. Note that the matrix_dumb_dumb is in row-major order.
 */
double get_matrix_dumb(matrix_dumb *mat, int row, int col)
{
    // Task 1.1 TODO
    int index = mat->cols * row + col;
    return mat->data[index];
}

/*
 * Sets the value at the given row and column to val. You may assume `row` and
 * `col` are valid. Note that the matrix_dumb is in row-major order.
 */
void set_matrix_dumb(matrix_dumb *mat, int row, int col, double val)
{
    int index = mat->cols * row + col;
    mat->data[index] = val;
}

/*
 * Allocates space for a matrix_dumb struct pointed to by the double pointer mat with
 * `rows` rows and `cols` columns. You should also allocate memory for the data array
 * and initialize all entries to be zeros. `parent` should be set to NULL to indicate that
 * this matrix_dumb is not a slice. You should also set `ref_cnt` to 1.
 * You should return -1 if either `rows` or `cols` or both have invalid values. Return -2 if any
 * call to allocate memory in this function fails.
 * Return 0 upon success.
 */
int allocate_matrix_dumb(matrix_dumb **mat, int rows, int cols)
{
    // Task 1.2 TODO
    // HINTS: Follow these steps.
    // 1. Check if the dimensions are valid. Return -1 if either dimension is not positive.
    // 2. Allocate space for the new matrix_dumb struct. Return -2 if allocating memory failed.
    // 3. Allocate space for the matrix_dumb data, initializing all entries to be 0. Return -2 if allocating memory failed.
    // 4. Set the number of rows and columns in the matrix_dumb struct according to the arguments provided.
    // 5. Set the `parent` field to NULL, since this matrix_dumb was not created from a slice.
    // 6. Set the `ref_cnt` field to 1.
    // 7. Store the address of the allocated matrix_dumb struct at the location `mat` is pointing at.
    // 8. Return 0 upon success.
    if (rows < 1 || cols < 1)
    {
        return -1;
    }
    matrix_dumb *new_matrix = (matrix_dumb *)malloc(sizeof(matrix_dumb));
    if (new_matrix == NULL)
    {
        return -2;
    }

    // constructor
    new_matrix->rows = rows;
    new_matrix->cols = cols;
    new_matrix->data = (double *)calloc(rows * cols, sizeof(double));
    if (new_matrix->data == NULL)
    {
        return -2;
    }
    new_matrix->ref_cnt = 1;
    new_matrix->parent = NULL;
    *mat = new_matrix;
    return 0;
}

/*
 * You need to make sure that you only free `mat->data` if `mat` is not a slice and has no existing slices,
 * or that you free `mat->parent->data` if `mat` is the last existing slice of its parent matrix_dumb and its parent
 * matrix_dumb has no other references (including itself).
 */
void deallocate_matrix_dumb(matrix_dumb *mat)
{
    // Task 1.3 TODO
    // HINTS: Follow these steps.
    // 1. If the matrix_dumb pointer `mat` is NULL, return.
    // 2. If `mat` has no parent: decrement its `ref_cnt` field by 1. If the `ref_cnt` field becomes 0, then free `mat` and its `data` field.
    // 3. Otherwise, recursively call `deallocate_matrix_dumb` on `mat`'s parent, then free `mat`.
    if (mat == NULL)
    {
        return;
    }
    if (mat->parent == NULL)
    {
        mat->ref_cnt--;
        if (mat->ref_cnt == 0)
        {
            free(mat->data);
            free(mat);
        }
    }
    else
    {
        deallocate_matrix_dumb(mat->parent);
        free(mat);
    }
}

/*
 * Allocates space for a matrix_dumb struct pointed to by `mat` with `rows` rows and `cols` columns.
 * Its data should point to the `offset`th entry of `from`'s data (you do not need to allocate memory)
 * for the data field. `parent` should be set to `from` to indicate this matrix_dumb is a slice of `from`
 * and the reference counter for `from` should be incremented. Lastly, do not forget to set the
 * matrix_dumb's row and column values as well.
 * You should return -1 if either `rows` or `cols` or both have invalid values. Return -2 if any
 * call to allocate memory in this function fails.
 * Return 0 upon success.
 * NOTE: Here we're allocating a matrix_dumb struct that refers to already allocated data, so
 * there is no need to allocate space for matrix_dumb data.
 */
int allocate_matrix_dumb_ref(matrix_dumb **mat, matrix_dumb *from, int offset, int rows, int cols)
{
    // Task 1.4 TODO
    // HINTS: Follow these steps.
    // 1. Check if the dimensions are valid. Return -1 if either dimension is not positive.
    // 2. Allocate space for the new matrix_dumb struct. Return -2 if allocating memory failed.
    // 3. Set the `data` field of the new struct to be the `data` field of the `from` struct plus `offset`.
    // 4. Set the number of rows and columns in the new struct according to the arguments provided.
    // 5. Set the `parent` field of the new struct to the `from` struct pointer.
    // 6. Increment the `ref_cnt` field of the `from` struct by 1.
    // 7. Store the address of the allocated matrix_dumb struct at the location `mat` is pointing at.
    // 8. Return 0 upon success.

    if (rows <= 0 || cols <= 0)
    {
        return -1;
    }

    (*mat) = (matrix_dumb *)malloc(sizeof(matrix_dumb));
    if (!(*mat))
    {
        return -2;
    }

    (*mat)->data = from->data + offset;
    (*mat)->rows = rows;
    (*mat)->cols = cols;
    (*mat)->parent = from;
    from->ref_cnt++;

    return 0;
}

/*
 * Sets all entries in mat to val. Note that the matrix_dumb is in row-major order.
 */
void fill_matrix_dumb(matrix_dumb *mat, double val)
{
    // Task 1.5 TODO
    for (int i = 0; i < mat->rows * mat->cols; i++)
    {
        mat->data[i] = val;
    }
}

/*
 * Store the result of taking the absolute value element-wise to `result`.
 * Return 0 upon success.
 * Note that the matrix_dumb is in row-major order.
 */
int abs_matrix_dumb(matrix_dumb *result, matrix_dumb *mat)
{
    // Task 1.5 TODO
    if (result->cols != mat->cols || result->rows != mat->rows || result->data == NULL)
    {
        return -1;
    }
    for (int i = 0; i < mat->rows; i++) {
        for (int j = 0; j < mat->rows; j++) {
            int index = i * mat->cols;
            if (mat->data[j + index] < 0) {
                result->data[j + index] = -mat->data[j + index];
            } else {
                result->data[j + index] = mat->data[j + index];
            }
        }
        
    }
    return 0;
}

/*
 * (OPTIONAL)
 * Store the result of element-wise negating mat's entries to `result`.
 * Return 0 upon success.
 * Note that the matrix_dumb is in row-major order.
 */
int neg_matrix_dumb(matrix_dumb *result, matrix_dumb *mat)
{
    if (result->cols != mat->cols || result->rows != mat->rows || result->data == NULL)
    {
        return -1;
    }
    for (int j = 0; j < mat->cols; j++) {
        for (int i = 0; i < mat->rows; i++) {
            int index = i * mat->cols;
            result->data[j + index] = -mat->data[j + index];
        }
    }
    return 0;
}

/*
 * Store the result of adding mat1 and mat2 to `result`.
 * Return 0 upon success.
 * You may assume `mat1` and `mat2` have the same dimensions.
 * Note that the matrix_dumb is in row-major order.
 */
int add_matrix_dumb(matrix_dumb *result, matrix_dumb *mat1, matrix_dumb *mat2)
{
    // Task 1.5 TODO
    if (result->cols != mat1->cols || result->rows != mat1->rows || result->data == NULL)
    {
        return -1;
    }
    if (mat1->cols != mat2->cols || mat1->rows != mat2->rows)
    {
        return -1;
    }
    for (int j = 0; j < mat1->rows; j++) {
        for (int i = 0; i < mat1->cols; i++) {
            result->data[i * mat1->cols + j] = mat1->data[i * mat1->cols + j] + mat2->data[i * mat1->cols + j];
        }
        
    }
    
    return 0;
}

/*
 * (OPTIONAL)
 * Store the result of subtracting mat2 from mat1 to `result`.
 * Return 0 upon success.
 * You may assume `mat1` and `mat2` have the same dimensions.
 * Note that the matrix_dumb is in row-major order.
 */
int sub_matrix_dumb(matrix_dumb *result, matrix_dumb *mat1, matrix_dumb *mat2)
{
    // Task 1.5 TODO
    if (result->cols != mat1->cols || result->rows != mat1->rows || result->data == NULL)
    {
        return -1;
    }
    if (mat1->cols != mat2->cols || mat1->rows != mat2->rows)
    {
        return -1;
    }
    for (int j = 0; j < mat1->cols; j++) {
        for (int i = 0; i < mat1->rows; i++) {
            int index = i * mat1->cols;
            result->data[j + index] = mat1->data[j + index] - mat2->data[j + index];
        }
        
    }
    return 0;
}

/*
 * Store the result of multiplying mat1 and mat2 to `result`.
 * Return 0 upon success.
 * Remember that matrix_dumb multiplication is not the same as multiplying individual elements.
 * You may assume `mat1`'s number of columns is equal to `mat2`'s number of rows.
 * Note that the matrix_dumb is in row-major order.
 */
int mul_matrix_dumb(matrix_dumb *result, matrix_dumb *mat1, matrix_dumb *mat2)
{
    // Task 1.6 TODO
    if (mat1->rows != result->rows || mat2->cols != result->cols || mat1->cols != mat2->rows || result->data == NULL)
    {
        return -1;
    }
    for (int i = 0; i < mat1->rows; i++)
    {
        for (int j = 0; j < mat2->cols; j++)
        {
            for (int k = 0; k < mat1->cols; k++)
            {
                result->data[i * mat2->cols + j] += mat1->data[i * mat1->cols + k] * mat2->data[k * mat2->cols + j];
            }
        }
    }
    return 0;
}

/*
 * Store the result of raising mat to the (pow)th power to `result`.
 * Return 0 upon success.
 * Remember that pow is defined with matrix_dumb multiplication, not element-wise multiplication.
 * You may assume `mat` is a square matrix_dumb and `pow` is a non-negative integer.
 * Note that the matrix_dumb is in row-major order.
 */
int pow_matrix_dumb(matrix_dumb *result, matrix_dumb *mat, int pow)
{
    // Task 1.6 TODO
    if (mat->cols != mat->rows || result->rows != result->cols || result->data == NULL || pow < 0)
    {
        return -1;
    }
    if (pow == 0)
    {
        for (int i = 0; i < mat->rows; i++)
        {
            for (int j = 0; j < mat->cols; j++)
            {
                if (i == j)
                {
                    result->data[i * mat->cols + j] = 1.0;
                }
                else
                {
                    result->data[i * mat->cols + j] = 0.0;
                }
            }
        }
    }
    else if (pow == 1)
    {
        add_matrix_dumb(result, result, mat);
    }
    else
    {
        matrix_dumb *temp = NULL;
        matrix_dumb *temp_zero = NULL;
        allocate_matrix_dumb(&temp, mat->rows, mat->cols);
        allocate_matrix_dumb(&temp_zero, mat->rows, mat->cols);
        add_matrix_dumb(result, temp_zero, mat);
        for (int i = 1; i < pow; i++)
        {
            add_matrix_dumb(temp, temp_zero, result);
            fill_matrix_dumb(result, 0);
            mul_matrix_dumb(result, temp, mat);
        }
        deallocate_matrix_dumb(temp);
        deallocate_matrix_dumb(temp_zero);
    }


    return 0;
}