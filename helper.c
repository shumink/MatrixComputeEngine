#include <xmmintrin.h>
#include "helper.h"

extern int g_seed;

extern ssize_t g_width;
extern ssize_t g_height;
extern ssize_t g_elements;

extern ssize_t g_nthreads;

void* para_sum(void* matrix) {
	float* arr = (float*)matrix;

	ssize_t work_per_thread = (int)floor((float)g_elements / (float)g_nthreads);
	float* local_sum = (float*)malloc(sizeof(float));
	register float temp_sum = 0;
	float partialsum[4] = {0};
	__m128 resultvector, working;
	resultvector = _mm_load_ps(&partialsum[0]);
	int k = 0;
	for (k = 0; k < work_per_thread; k += 4) {
		working = _mm_load_ps(&arr[k]);
		resultvector = _mm_add_ps(working, resultvector);
	}
	_mm_store_ps(&partialsum[0], resultvector);
	for (int i = 0; i < 4; i++) {
		temp_sum += partialsum[i];
	}

	for (int i = k - 4; k < work_per_thread % 4; k++) {
		temp_sum += arr[i];
	}
	*local_sum = temp_sum;
	pthread_exit((void*)local_sum);
}

void* para_matrix_add(void* args) {
	matrix_group* group = (matrix_group*)args;
	ssize_t work_per_thread = g_elements / g_nthreads;
	float* result = group->result;
	const float* a = group->a;
	const float* b = group->b;
	__m128 matrixa, matrixb, resultvector;
	int k = 0;
	for (k = 0; k < work_per_thread; k += 4) {
		matrixa = _mm_load_ps(&a[k]);
		matrixb = _mm_load_ps(&b[k]);
		resultvector = _mm_add_ps(matrixa, matrixb);
		_mm_store_ps(&result[k], resultvector);
	}
	for (int j = k - 4; j < work_per_thread % 4; j++) {
		result[j] = a[j] + b[j];
	}
	pthread_exit(NULL);
}
void* para_matrix_mul(void* args) {
	int threads_available = (g_width < 4)? 1: g_nthreads;
	const float* X = ((matrix_group*)args)->a;
	const float* Y = ((matrix_group*)args)->b;
	float* result = ((matrix_group*)args)->result;
	int start_i = ((matrix_group*)args)->value;
	ssize_t work_per_thread = g_elements / threads_available;
	int row = 0;
	int column = 0;
	register float element = 0;

	for (int i = start_i; i < start_i + work_per_thread; i++) {
		row = floor(i / g_width);
		column = i % g_width;
		element = 0;
		__m128 hvector, vvector, mul, D;
		float d[4] = {0};
		int k = 0;
		for (k = 0; k < g_width; k += 4) {
			hvector = _mm_load_ps(&X[IDX(row, k)]);
	        vvector = _mm_load_ps(&Y[IDX(column, k)]);
	        mul = _mm_mul_ps(hvector, vvector);
	        D = _mm_load_ps(&d[0]);
	        D = _mm_add_ps(mul, D);
	        _mm_store_ps(&d[0], D);
		}
		for (int j = 0; j < 4; j++) {
			element += d[j];
		}
		for (int j = k - 4; j < work_per_thread % 4; j++) {
			element += X[IDX(row,j)] * Y[IDX(column,j)];
		}
		result[i] = element;
	}
	pthread_exit(NULL);
}

void* para_scalar_add(void* args) {
	matrix_group* matrices = (matrix_group*)args;
	const float* old = matrices->a;
	float* new = matrices->result;
	float value = matrices->value;
	ssize_t work_per_thread = g_elements / g_nthreads;
	__m128 oldvector, D, newvector;
	float d[4] = {value, value, value, value};
	int k = 0;
	for (k = 0; k < work_per_thread; k += 4) {
		oldvector = _mm_load_ps(&old[k]);
		D = _mm_load_ps(&d[0]);
		newvector = _mm_add_ps(oldvector, D);
		_mm_store_ps(&new[k], newvector);
	}
	for (int j = k - 4; j < work_per_thread % 4; j++) {
		new[j] = old[j] + value;
	}

	pthread_exit(NULL);
}

void* para_identity(void* args) {
	matrix_group* matrices = (matrix_group*)args;
	float* result = matrices->result;
	int start_i = (int)matrices->value;
	ssize_t work_per_thread = g_width / g_nthreads;
	for (int i = start_i; i < start_i + work_per_thread; i++) {
		result[IDX(i, i)] = 1;
	}
	pthread_exit(NULL);
}

void* para_uniform(void* args) {
	matrix_group* matrices = (matrix_group*)args;
	float* result = matrices->result;
	float value = matrices->value;
	ssize_t work_per_thread = g_elements / g_nthreads;
	for (int i = 0; i < work_per_thread; i++) {
		result[i] = value;
	}
	pthread_exit(NULL);
}

void* para_clone(void* args) {
	const float* source = ((matrix_group*)args)->a;
	float* result = ((matrix_group*)args)->result;
	int work_per_thread = g_elements / g_nthreads;
	memmove(result, source, sizeof(float) * work_per_thread);
	pthread_exit(NULL);
}

void* para_qsort(void* args) {
	qsort(((Merge_args*)args)->array,
	((Merge_args*)args)->length, sizeof(float), compare);
	pthread_exit(NULL);
}

void* para_reversed(void* args) {
	float* result = ((matrix_group*)args)->result;
	int start = ((matrix_group*)args)->value;
	const float* source = ((matrix_group*)args)->a;
	int work_per_thread = g_elements / g_nthreads;
	for (int i = start; i < start + work_per_thread; i++) {
		result[i] = source[g_elements - 1 - i];
	}
	pthread_exit(NULL);
}

void* para_transposed(void* args) {
	int threads_available = (g_width < 4)? 1: g_nthreads;
	const float* matrix_a = ((submatrix*)args)->matrix_a;
	float* result = ((submatrix*)args)->result;
	int start_i = ((submatrix*)args)->start_i;
	int start_j = ((submatrix*)args)->start_j;
	ssize_t work_per_thread = g_width / (int)sqrt(threads_available);



	for (int j = start_j; j < start_j + work_per_thread; j++) {
		for (int i = start_i; i < start_i + work_per_thread; i++) {
			result[IDX(j, i)] = matrix_a[IDX(i, j)];
		}

	}
	pthread_exit(NULL);
}

void* para_scalar_mul(void* args) {
	matrix_group* matrices = (matrix_group*)args;
	const float* old = matrices->a;
	float* new = matrices->result;
	float value = matrices->value;
	ssize_t work_per_thread = g_elements / g_nthreads;
	__m128 oldvector, D, newvector;
	float d[4] = {value, value, value, value};
	int k = 0;
	for (k = 0; k < work_per_thread; k += 4) {
		oldvector = _mm_load_ps(&old[k]);
		D = _mm_load_ps(&d[0]);
		newvector = _mm_mul_ps(oldvector, D);
		_mm_store_ps(&new[k], newvector);
	}
	for (int j = k - 4; j < work_per_thread % 4; j++) {
		new[j] = old[j] * value;
	}
	pthread_exit(NULL);
}

void* para_trace(void* args) {
	matrix_group* received = (matrix_group*)args;
	const float* matrix= received->a;
	int start_index = received->value;
	ssize_t work_per_thread = start_index + g_width / g_nthreads;
	float* local_trace = (float*)malloc(sizeof(float));
	register float temp_trace = 0;

	for (int i = start_index; i < work_per_thread; i++) {
		temp_trace += matrix[IDX(i, i)];
	}
	*local_trace = temp_trace;
	pthread_exit((void*)local_trace);
}

void* para_min(void* matrix) {
	float* arr = (float*)matrix;

	ssize_t work_per_thread = g_elements / g_nthreads;
	float* local_min = (float*)malloc(sizeof(float));
	register float min = arr[0];

	__m128 minvector, working;
	float d[4] = {arr[0], arr[0], arr[0], arr[0]};
	minvector = _mm_load_ps(&d[0]);
	int k = 0;
	for (k = 0; k < work_per_thread; k += 4) {
		working = _mm_load_ps(&arr[k]);
		minvector = _mm_min_ps(working, minvector);
	}
	_mm_store_ps(&d[0], minvector);
	for (int i = 0; i < 4; i++) {
		if (min > d[i]) {
			min = d[i];
		}
	}
	for (int i = k - 4; i < work_per_thread; i++) {
		if (min > arr[i]) {
			min = arr[i];
		}
	}

	*local_min = min;
	pthread_exit((void*)local_min);
}

void* para_max(void* matrix) {
	float* arr = (float*)matrix;

	ssize_t work_per_thread = g_elements / g_nthreads;
	float* local_max = (float*)malloc(sizeof(float));
	register float max = arr[0];
	__m128 maxvector, working;
	float d[4] = {arr[0], arr[0], arr[0], arr[0]};
	maxvector = _mm_load_ps(&d[0]);
	int k = 0;
	for (k = 0; k < work_per_thread; k += 4) {
		working = _mm_load_ps(&arr[k]);
		maxvector = _mm_max_ps(working, maxvector);
	}
	_mm_store_ps(&d[0], maxvector);
	for (int i = 0; i < 4; i++) {
		if (max < d[i]) {
			max = d[i];
		}
	}
	for (int i = k - 4; i < work_per_thread; i++) {
		if (max < arr[i]) {
			max = arr[i];
		}
	}
	*local_max = max;
	pthread_exit((void*)local_max);
}

void* para_determinant_triangle(void* args) {
	const double* matrix = ((double_matrix_group*)args)->a;
	int start_i = ((double_matrix_group*)args)->value;
	double* local_result = calloc(1, sizeof(double));
	*local_result = 1;
	ssize_t work_per_thread = g_width / g_nthreads;
	for (int i = start_i; i < start_i + work_per_thread; i++) {
		*local_result *= matrix[IDX(i, i)];
	}
	pthread_exit((void*)local_result);
}

float dete_triangle(float* triangle) {
	float result = 1;
	for (int i = 0; i < g_width; i++) {
		result *= triangle[IDX(i, i)];
	}
	return result;
}

void* para_freq(void* args) {
	const float* matrix = ((matrix_group*)args)->a;
	const float value = ((matrix_group*)args)->value;
	ssize_t work_per_thread = g_elements / g_nthreads;
	ssize_t* local_freq = (ssize_t*)malloc(sizeof(ssize_t));
	register int frequency = 0;
	__m128 matrixvector, vvector, resultvector, TR, temp, tempresult;
	float resultarr[4] = {0, 0, 0, 0};
	float valuearr[4] = {value, value, value, value};
	float trim[4] = {1, 1, 1, 1};


	vvector = _mm_load_ps(&valuearr[0]);
	TR = _mm_load_ps(&trim[0]);

	int k = 0;
	for (k = 0; k < work_per_thread; k += 4) {
		resultvector = _mm_load_ps(&resultarr[0]);
		matrixvector = _mm_load_ps(&matrix[k]);
	    temp = _mm_cmpeq_ps(matrixvector, vvector);
	    tempresult = _mm_and_ps(temp, TR);
	    resultvector = _mm_add_ps(resultvector, tempresult);
		_mm_store_ps(&resultarr[0], resultvector);
	}

	for (int i = 0; i < 4; i++) {
		frequency += resultarr[i];
	}

	for (int i = k - 4; i < work_per_thread % 4; i++) {
		if (matrix[i] == value) {
			frequency++;
		}
	}
	*local_freq = frequency;
	pthread_exit((void*)local_freq);
}

//Laplace expansion
//Adapted from: http://www.vikparuchuri.com/blog/find-the-determinant-of-a-matrix/
float laplace(const float* matrix, int length) {
  if (length == 2) {
    return matrix[0] * matrix[3] - matrix[1] * matrix[2];
  }
  float result = 0;
  for (int x = 0; x < length; x++) {
    float* smaller = calloc((length - 1) * (length - 1), sizeof(float));
    load_matrix(smaller, matrix, x, length - 1);
    result += pow(-1, x) * matrix[x] * laplace (smaller, length - 1);
    free(smaller);
  }
  return result;
}

void load_matrix(float* small, const float* big, int deletedcolumn, int length) {
  for (int i = 0; i < length; i++) {
    for (int j = 0; j < length; j++) {
      if (j  >= deletedcolumn) {
        small[i * length + j] = big[(i + 1) * (length + 1) + j + 1];
      } else {
        small[i * length + j] = big[(i + 1) * (length + 1) + j];
      }
    }
  }
}

void Merge(float* left, float* right, float* dst, int leftLength, int rightLength) {
	int arr1 = 0;
	int arr2 = 0;
	while (arr1 < leftLength && arr2 < rightLength) {
		if (left[arr1] >= right[arr2]) {
			dst[arr1 + arr2] = right[arr2];
			arr2++;
		} else {
			dst[arr1 + arr2] = left[arr1];
			arr1++;
		}
	}
	if (arr1 < leftLength) {
		memcpy(dst + arr1 + arr2, left + arr1, (leftLength - arr1) * sizeof(float));
	} else if (arr2 < rightLength) {
		memcpy(dst + arr1 + arr2, right + arr2, (rightLength - arr2) * sizeof(float));
	}
}

int compare(const void* a, const void* b) {
	return *((const float*)a) - *((const float*)b);
}

matrix_group* form_group(const float* a, const float* b, float* result,float value) {
	matrix_group* new_group = (matrix_group*)malloc(sizeof(matrix_group) * g_nthreads);
	new_group->value = value;
	new_group->a = a;
	new_group->b = b;
	new_group->result = result;
	return new_group;
}
