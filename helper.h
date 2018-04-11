#ifndef MATRIX_HELPER_H
#define MATRIX_HELPER_H
#include "structs.h"

//Utility functions
matrix_group* form_group(const float* a, const float* b, float* result,float value);
void Merge(float* left, float* right, float* dst, int leftLength, int rightLength);
int compare(const void* a, const void* b);
float* recursive_power(float* matrix, int exponent);

//LU decomposition
void LU_decomposition (const float* A, double* L, double* U);
float dete_triangle(float* triangle);
void load_matrix(float* small, const float* big, int deletedcolumn, int length);

//Matrix scalar operations
void* para_scalar_add(void* args);
void* para_scalar_mul(void* args);

//Matrix element-wise operations
void* para_matrix_add(void* args);

//Matrix attributes
void* para_sum(void* matrix);
void* para_min(void* matrix);
void* para_max(void* matrix);
void* para_trace(void* args);
void* para_freq(void* args);

//Matrix generation
void* para_identity(void* args);
void* para_reversed(void* args);
void* para_uniform(void* args);
void* para_clone(void* args);
void* para_transposed(void* args);

//Matrix computations
void* para_determinant_triangle(void* args);
void* para_cannon_mul(void* args);
void* para_qsort(void* args);
void* para_matrix_mul(void* args);

//Useful macro
#define IDX(x, y) ((x) * g_width + y)

#endif
