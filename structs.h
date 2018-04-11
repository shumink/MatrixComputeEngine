#ifndef STRUCTS_H
#define STRUCTS_H
struct matrix_group {
	float value;
	const float* a;
	const float* b;
	float* result;
	char padding[64 - 3 * sizeof(float*) - sizeof(float)];
};
typedef struct matrix_group matrix_group;

struct submatrix {
	const float* matrix_a;
	const float* matrix_b;
	float* result;
	int working_i;
	int working_j;
	int start_i;
	int start_j;
};
typedef struct submatrix submatrix;

struct sub_lu_matrix {
	const float* a;
	double* u;
	double* l;
	int working_i;
	int working_j;
	int start_i;
	int start_j;
};
typedef struct sub_lu_matrix sub_lu_matrix;

struct double_matrix_group {
	double value;
	const double* a;
	const double* b;
	double* result;
	char padding[64 - 3 * sizeof(double*) - sizeof(double)];
};
typedef struct double_matrix_group double_matrix_group;

struct Merge_args {
	float* array;
	int length;
	int depth;
};
typedef struct Merge_args Merge_args;
#endif
