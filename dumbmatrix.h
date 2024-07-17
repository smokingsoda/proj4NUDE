
typedef struct matrix_dumb
{
    int rows;              // number of rows
    int cols;              // number of columns
    double *data;          // pointer to rows * columns doubles
    int ref_cnt;           // How many slices/matrices are referring to this matrix's data
    struct matrix_dumb *parent; // NULL if matrix is not a slice, else the parent matrix of the slice
} matrix_dumb;

double rand_double_dumb(double low, double high);
void rand_matrix_dumb(matrix_dumb *result, unsigned int seed, double low, double high);
int allocate_matrix_dumb(matrix_dumb **mat, int rows, int cols);
int allocate_matrix_dumb_ref(matrix_dumb **mat, matrix_dumb *from, int offset, int rows, int cols);
void deallocate_matrix_dumb(matrix_dumb *mat);
double get_matrix_dumb(matrix_dumb *mat, int row, int col);
void set_matrix_dumb(matrix_dumb *mat, int row, int col, double val);
void fill_matrix_dumb(matrix_dumb *mat, double val);
int add_matrix_dumb(matrix_dumb *result, matrix_dumb *mat1, matrix_dumb *mat2);
int sub_matrix_dumb(matrix_dumb *result, matrix_dumb *mat1, matrix_dumb *mat2);
int mul_matrix_dumb(matrix_dumb *result, matrix_dumb *mat1, matrix_dumb *mat2);
int pow_matrix_dumb(matrix_dumb *result, matrix_dumb *mat, int pow);
int neg_matrix_dumb(matrix_dumb *result, matrix_dumb *mat);
int abs_matrix_dumb(matrix_dumb *result, matrix_dumb *mat);