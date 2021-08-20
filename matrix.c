#include "matrix.h"
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

/*
 * Generates a random double between `low` and `high`.
 */
double rand_double(double low, double high) {
    double range = (high - low);
    double div = RAND_MAX / range;
    return low + (rand() / div);
}

/*
 * Generates a random matrix with `seed`.
 */
void rand_matrix(matrix *result, unsigned int seed, double low, double high) {
    srand(seed);
    for (int i = 0; i < result->rows; i++) {
        for (int j = 0; j < result->cols; j++) {
            set(result, i, j, rand_double(low, high));
        }
    }
}

/*
 * Allocate space for a matrix struct pointed to by the double pointer mat with
 * `rows` rows and `cols` columns. You should also allocate memory for the data array
 * and initialize all entries to be zeros. Remember to set all fieds of the matrix struct.
 * `parent` should be set to NULL to indicate that this matrix is not a slice.
 * You should return -1 if either `rows` or `cols` or both have invalid values, or if any
 * call to allocate memory in this function fails. If you don't set python error messages here upon
 * failure, then remember to set it in numc.c.
 * Return 0 upon success and non-zero upon failure.
 */

/* Allocates memory for a single matrix pointer */
int allocate_single(matrix *mat, int rows, int cols) {
    if (rows <= 0 || cols <= 0)
    {
        return -1;
    }
    double **temp = (double **) malloc ( rows * sizeof(double *));
    if(temp == NULL) {
        return -1;
    }
    /* malloc */
    for (int i = 0; i < rows; ++i)
    {
        temp[i] = (double *) calloc (cols, sizeof(double));
        if (temp[i] == NULL) {
            return -1;
        }
    } 
    if (rows == 1 || cols == 1)
    {
        mat->is_1d = 1;
    }
    else {
        mat->is_1d = 0;
    }
    mat->rows = rows;
    mat->cols = cols;
    mat->data = temp;
    mat->parent = NULL;
    mat->ref_cnt = 1;
    return 0;
}

int allocate_matrix(matrix **mat, int rows, int cols) {
    /* TODO: YOUR CODE HERE */
    matrix *temp = (matrix *) malloc(sizeof(matrix));
    if (temp == NULL) {
        return -1;
    }
    int status = allocate_single(temp, rows, cols);
    *mat = temp;
    return status;
}

/*
 * Allocate space for a matrix struct pointed to by `mat` with `rows` rows and `cols` columns.
 * This is equivalent to setting the new matrix to be
 * from[row_offset:row_offset + rows, col_offset:col_offset + cols]
 * If you don't set python error messages here upon failure, then remember to set it in numc.c.
 * Return 0 upon success and non-zero upon failure.
 */
int allocate_matrix_ref(matrix **mat, matrix *from, int row_offset, int col_offset,
                        int rows, int cols) {
    /* TODO: YOUR CODE HERE */
    int end_r = rows + row_offset;
    int end_c = cols + col_offset;
    if (row_offset < 0 || col_offset < 0 || rows < 1 || cols < 1) {
        // invalid inputs.
        return -1;
    }
    if (end_r > from->rows || end_c > from->cols) {
        // out of bounds.
        return -1;
    }
    matrix *temp = (matrix *) malloc(sizeof(matrix));
    if (temp == NULL) {printf("malloc temp failed\n"); return -1;}
    temp->is_1d = (end_r == 1 || end_c == 1);
    temp->rows = rows; temp->cols = cols;
    temp->parent = from; temp->ref_cnt = 1;
    double **temp_rows = (double **) malloc (sizeof(double *) * rows);
    if (temp_rows == NULL) {printf("malloc temp failed.\n"); return -1;}

    for (int i = 0; i < rows; i++)
    {
        temp_rows[i] = from->data[i+ row_offset];
    }
    temp->data = temp_rows;
    from->ref_cnt += 1;
    *mat = temp;

    return 0;

    /* Kenny's code
    int r = rows + row_offset;
    int c = cols + col_offset;

    if (r > from->rows || c > from->cols)
    {
        return -1;
    }
    matrix *temp = (matrix *) malloc(sizeof(matrix) * rows + 1);
    if (temp == NULL) {
        return -1;
    }
    if (r == 1 || c == 1)
    {
        temp->is_1d = 1;
    }
    else 
    {
        temp->is_1d = 0;
    }
    temp->rows = rows; temp->cols = cols;
    temp->parent = from; temp->ref_cnt = 1;

    temp->data = (double **) malloc(sizeof(double *) * rows);
    if(temp->data == NULL) {
        return -1;
    }
    for (int i = 0; i < rows; ++i)
    {
        temp->data[i] = (double *) malloc(sizeof(double) * cols);
        if (temp->data[i] == NULL) {
            free(temp->data);
            return -1;
        }
        temp->data[i] = from->data[i + row_offset];
    }

    from->ref_cnt += 1;
    *mat = temp;
    return 0;
     */ // KENNYS ORIGINAL CODE.
}

/*
 * This function will be called automatically by Python when a numc matrix loses all of its
 * reference pointers.
 * You need to make sure that you only free `mat->data` if no other existing matrices are also
 * referring this data array.
 * See the spec for more information.
 */
void deallocate_matrix(matrix *mat) {
    /* TODO: YOUR CODE HERE */
    if (mat == NULL) {
        return;
    }
    if (mat->parent == NULL && mat->ref_cnt <= 1)
    {
        for (int i = 0; i < mat->rows; ++i)
        {
            free(mat->data[i]);
        }
    }
    free(mat->data);
    free(mat);
}
/** Copy value of mat to result */
void mat_copy(matrix *result, matrix *mat) {
    for (int i = 0; i < result->rows; ++i) {
        for (int j = 0; j < result->cols; ++j) {
            result->data[i][j] = mat->data[i][j];
        }
    }
}

/*
 * Return the double value of the matrix at the given row and column.
 * You may assume `row` and `col` are valid.
 */
double get(matrix *mat, int row, int col) {
    /* TODO: YOUR CODE HERE */
    return mat->data[row][col];
}

/*
 * Set the value at the given row and column to val. You may assume `row` and
 * `col` are valid
 */
void set(matrix *mat, int row, int col, double val) {
    /* TODO: YOUR CODE HERE */
    mat->data[row][col] = val;
}

/*
 * Set all entries in mat to val
 */
void fill_matrix(matrix *mat, double val) {
    /* TODO: YOUR CODE HERE
    if (mat->rows + mat->cols < 200)
    {
        for (int i = 0; i < mat->rows; ++i) {
            for (int j = 0; j < mat->cols; ++j) {
            mat->data[i][j] = val;
            }
        }   
    }
    int r = mat->rows / 4 * 4;
    int c = mat->cols / 4 * 4;
    int full = r * c;
    int i; int j;
    __m256d vals = _mm256_set1_pd (val);
    #pragma omp parallel for
    for (int z = 0; z < full; z+=4) {
        i = z / r;
        j = z % r;
        _mm256_storeu_pd((mat->data[i] + j), vals);
        }
    for (int x = r; x < mat->rows; ++x)
    {
        for (int y = c; y < mat->cols; ++y)
        {
            mat->data[x][y] = val;
        }
    } */
    int mat_r = mat->rows;
    int mat_c = mat->cols;
    if (mat_r + mat_c < 200) {
        for (int i = 0; i < mat_r; ++i) {
            for (int j = 0; j < mat_c; ++j) {
                mat->data[i][j] = val;
                }
            }
    }
    int i; int j;
    #pragma omp parallel for
    for (int x = 0; x < mat_r * mat_c; ++x)
    {
        i = x / mat_c;
        j = x % mat_c;
        mat->data[i][j] = val;
    }
}

/*
 * Store the result of adding mat1 and mat2 to `result`.
 * Return 0 upon success and a nonzero value upon failure.
 */
int add_matrix(matrix *result, matrix *mat1, matrix *mat2) {
    /* TODO: YOUR CODE HERE */
    /**
    if (mat1->rows % 4 != 0 || mat1->cols % 4 != 0 || mat1->rows + mat1->cols < 1000)
    {
        for (int i = 0; i < mat1->rows; ++i) {
            for (int j = 0; j < mat1->cols; ++j) {
                result->data[i][j] = mat1->data[i][j] + mat2->data[i][j];
            }
        }   
        return 0;
    }
    int orig_r = mat1->rows;
    int orig_c = mat1->cols;
    int r = orig_r / 4 * 4;
    int c = orig_c / 4 * 4;
    int full = r * c;
    int i; int j;

    #pragma omp parallel for
    for (int z = 0; z < full; z+=4) {
        i = z / orig_c;
        j = z % orig_c;
        _mm256_storeu_pd((result->data[i] + j), _mm256_add_pd(_mm256_loadu_pd(mat1->data[i] + j), _mm256_loadu_pd(mat2->data[i] + j)));
        }

    for (int f = 0; f < orig_r; ++f) {
            for (int g = c; g < orig_c; ++g) {
                printf("%d %d\n", f, g);
                result->data[f][g] = mat1->data[f][g] + mat2->data[f][g];
            }
        }
   printf("%f", mat1->data[996][998]); 
    for (int x = r; x < orig_r; ++x)
    {

        for (int y = 0; y < c; ++y)
        {
        printf("%d %d\n", x, y);
            result->data[x][y] = mat1->data[x][y] + mat2->data[x][y];
        }
    }
    */
    int mat_r = mat1->rows;
    int mat_c = mat1->cols;
    if (mat_r + mat_c < 200)
    {
            for (int i = 0; i < mat_r; ++i) {
                for (int j = 0; j < mat_c; ++j) {
                    result->data[i][j] = mat1->data[i][j] + mat2->data[i][j];
                }
            }
        return 0;
    }
    int i; int j;
    #pragma omp parallel for
    for (int x = 0; x < mat_r * mat_c; ++x)
    {
        i = x / mat_c;
        j = x % mat_c;
        result->data[i][j] = mat1->data[i][j] + mat2->data[i][j];
    }
    return 0;
}

/*
 * Store the result of subtracting mat2 from mat1 to `result`.
 * Return 0 upon success and a nonzero value upon failure.
 */
int sub_matrix(matrix *result, matrix *mat1, matrix *mat2) {
    int mat_r = mat1->rows;
    int mat_c = mat1->cols;
    if (mat_r + mat_c < 200)
    {
            for (int i = 0; i < mat_r; ++i) {
                for (int j = 0; j < mat_c; ++j) {
                    result->data[i][j] = mat1->data[i][j] - mat2->data[i][j];
                }
            }
        return 0;
    }
    int i; int j;
    #pragma omp parallel for
    for (int x = 0; x < mat_r * mat_c; ++x)
    {
        i = x / mat_c;
        j = x % mat_c;
        result->data[i][j] = mat1->data[i][j] - mat2->data[i][j];
    }
    return 0;
}

/** Compute a single multiplication of one row and one col with length len */
double single_product(double *n1, double*n2, int len) {
    double result[4];
    __m256d prod_v = _mm256_set1_pd (0);
    for(int i = 0; i < len/4 * 4; i += 4) {
        prod_v = _mm256_fmadd_pd(_mm256_loadu_pd(n2 + i), _mm256_loadu_pd(n1 + i), prod_v);
    }
    _mm256_storeu_pd(result, prod_v);
    for(int i = len/4 * 4; i < len; i++) {
        result[0] += n1[i] * n2[i];
    }
    return result[0] + result[1] + result[2] + result[3];
}

double **transpose(double **data, int data_rows, int data_cols) {
    double **temp = (double **) malloc (data_cols * sizeof(double *));
    for (int i = 0; i < data_cols; ++i)
    {
        temp[i] = (double *) malloc (data_rows * sizeof(double));
        
        for (int j = 0; j < data_rows; ++j)
        {
            temp[i][j] = data[j][i];
        }
    }
    return temp;
}

/*
 * Store the result of multiplying mat1 and mat2 to `result`.
 * Return 0 upon success and a nonzero value upon failure.
 * Remember that matrix multiplication is not the same as multiplying individual elements.
 */
int mul_matrix(matrix *result, matrix *mat1, matrix *mat2) {
    /* TODO: YOUR CODE HERE */
    int r1 = mat1->rows;
    int c1 = mat1->cols;
    int r2 = mat2->rows;
    int c2 = mat2->cols;
    if (r1 + c1 < 200)
    {
        double temp;
        for (int i = 0; i < r1; i++) {
            for (int j = 0; j < c2; j++) {
                temp = 0;
                for (int k = 0; k < c1; k++) {
                    temp += mat1->data[i][k] * mat2->data[k][j];
                }
                result->data[i][j] = temp;
            }
        }
        return 0;
    }
    double **mat2_t = transpose(mat2->data, r2, c2);
    int i; int j;
    #pragma omp parallel for
    for (int x = 0; x < r1 * c2; ++x)
    {
        i = x / c2;
        j = x % c2;
        result->data[i][j] = single_product(mat1->data[i], mat2_t[j], c1);
    }
    return 0;
}

/*
 * Store the result of raising mat to the (pow)th power to `result`.
 * Return 0 upon success and a nonzero value upon failure.
 * Remember that pow is defined with matrix multiplication, not element-wise multiplication.
 */
int pow_matrix(matrix *result, matrix *mat, int pow) {
    /* TODO: YOUR CODE HERE */
    if (pow == 1)
    {
        /* return copy of original matrix */
        mat_copy(result, mat);
        return 0;
    }
    int mat_r = mat->rows;
    int mat_c = mat->cols;
    int res_r = result->rows;
    //int res_c = result->cols;
    for (int i = 0; i < res_r; ++i)
        {
            result->data[i][i] = 1;
        }
    if (pow == 0)
    {
        /* identity matrix */
        return 0;
    }
    matrix *temp1 = NULL;
    allocate_matrix(&temp1, mat_r, mat_c);
    matrix *temp2 = NULL;
    allocate_matrix(&temp2, mat_r, mat_c);
    double **change;
    mat_copy(temp1, mat); //initialize matrix
    while (pow != 0) {
        if (pow & 1)
        {
            matrix *temp = NULL;
            allocate_matrix(&temp, mat_r, mat_c);
            mat_copy(temp, result);
            mul_matrix(result, temp, temp1);
        }
        mul_matrix(temp2, temp1, temp1);
        change = temp1->data;
        temp1->data = temp2->data;
        temp2->data = change;
        pow /= 2;
    }
    return 0;
}

/*
 * Store the result of element-wise negating mat's entries to `result`.
 * Return 0 upon success and a nonzero value upon failure.
 */
int neg_matrix(matrix *result, matrix *mat) {
    /**
    if (mat->rows % 4 != 0 || mat->cols != mat->rows || mat->rows + mat->cols < 1000)
    {
        for (int i = 0; i < mat->rows; ++i) {
            for (int j = 0; j < mat->cols; ++j) {
            result->data[i][j] = -mat->data[i][j];
            }
        }   
        return 0;
    }
    int r = mat->rows / 4 * 4;
    int c = mat->cols / 4 * 4;
    int full = r * c;
    int i; int j;
    __m256d zeros = _mm256_set1_pd (0);

    #pragma omp parallel for
    for (int z = 0; z < full; z+=4) {
        i = z / r;
        j = z % r;
        _mm256_storeu_pd((result->data[i] + j), _mm256_sub_pd(zeros, _mm256_loadu_pd(mat->data[i] + j)));
        }
    return 0; */
    int mat_r = mat->rows;
    int mat_c = mat->cols;
    if (mat_r + mat_c < 200)
    {
        for (int i = 0; i < mat_r; ++i)
        {
            for (int j = 0; j < mat_c; ++j)
            {
            result->data[i][j] = -mat->data[i][j];
            }
        }
    }
    int i; int j;
    #pragma omp parallel for
    for (int x = 0; x < mat_r * mat_c; ++x)
    {
        i = x / mat_c;
        j = x % mat_c;
        result->data[i][j] = -mat->data[i][j];
    }
    return 0;
}

/*
 * Store the result of taking the absolute value element-wise to `result`.
 * Return 0 upon success and a nonzero value upon failure.
 */
int abs_matrix(matrix *result, matrix *mat) {
    /* TODO: YOUR CODE HERE */
    for (int i = 0; i < mat->rows; ++i)
        {
            for (int j = 0; j < mat->cols; ++j)
            {
                double d = mat->data[i][j];
                if(d < 0) {
                    result->data[i][j] = -d;
                } else {
                    result->data[i][j] = d;
                }
            }
        }
        return 0;
}
