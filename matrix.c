#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <stdbool.h>
#include <pthread.h>
#include <string.h>
#include <math.h>
#include <errno.h>
#include "matrix.h"
#include "helper.c"

int g_seed = 0;

ssize_t g_width = 0;
ssize_t g_height = 0;
ssize_t g_elements = 0;

ssize_t g_nthreads = 1;

////////////////////////////////
///     UTILITY FUNCTIONS    ///
////////////////////////////////


/**
* Returns pseudorandom number determined by the seed.
*/
int fast_rand(void) {

	g_seed = (214013 * g_seed + 2531011);
	return (g_seed >> 16) & 0x7FFF;
}

/**
* Sets the seed used when generating pseudorandom numbers.
*/
void set_seed(int seed) {

	g_seed = seed;
}

/**
* Sets the number of threads available.
*/
void set_nthreads(ssize_t count) {

	g_nthreads = count;
}

/**
* Sets the dimensions of the matrix.
*/
void set_dimensions(ssize_t order) {

	g_width = order;
	g_height = order;

	g_elements = g_width * g_height;
}

/**
* Displays given matrix.
*/
void display(const float* matrix) {

	for (ssize_t y = 0; y < g_height; y++) {
		for (ssize_t x = 0; x < g_width; x++) {
			if (x > 0) printf(" ");
			printf("%.2f", matrix[y * g_width + x]);
		}

		printf("\n");
	}
}

/**
* Displays given matrix row.
*/
void display_row(const float* matrix, ssize_t row) {

	for (ssize_t x = 0; x < g_width; x++) {
		if (x > 0) printf(" ");
		printf("%.2f", matrix[row * g_width + x]);
	}

	printf("\n");
}

/**
* Displays given matrix column.
*/
void display_column(const float* matrix, ssize_t column) {

	for (ssize_t i = 0; i < g_height; i++) {
		printf("%.2f\n", matrix[i * g_width + column]);
	}
}

/**
* Displays the value stored at the given element index.
*/
void display_element(const float* matrix, ssize_t row, ssize_t column) {

	printf("%.2f\n", matrix[row * g_width + column]);
}

////////////////////////////////
///   MATRIX INITALISATIONS  ///
////////////////////////////////

/**
* Returns new matrix with all elements set to zero.
*/
float* new_matrix(void) {

	return calloc(g_elements, sizeof(float));
}

/**
* Returns new identity matrix.
*/
float* identity_matrix(void) {

	float* result = new_matrix();

	/*
	1 0
	0 1
	*/
	pthread_t id[g_nthreads];
	ssize_t work_per_thread = g_width / g_nthreads;
	matrix_group* sequence = (matrix_group*)malloc(sizeof(matrix_group) * g_nthreads);
	for (int i = 0; i < g_nthreads; i++) {
		sequence[i].value = i * work_per_thread;
		sequence[i].result = result;
		pthread_create(&id[i], NULL, para_identity, (void*)(sequence + i));
	}
	for (int i = 0; i < g_nthreads; i++) {
		pthread_join(id[i], NULL);
	}

	int remaining = work_per_thread * g_nthreads;
	for (int i = remaining; i < g_width; i++) {
		result[i] = 1;
	}

	free(sequence);

	return result;
}


/**
* Returns new matrix with elements generated at random using given seed.
*/
float* random_matrix(int seed) {

	float* matrix = new_matrix();

	set_seed(seed);
	for (ssize_t i = 0; i < g_elements; i++) {
		matrix[i] = fast_rand();
	}

	return matrix;
}

/**
* Returns new matrix with all elements set to given value.
*/
float* uniform_matrix(float value) {

	float* result = new_matrix();

	/*
	1 1
	1 => 1 1
	*/
	if (value == 0.0) {
		return result;
	}

	pthread_t id[g_nthreads];
	ssize_t work_per_thread = g_elements / g_nthreads;
	matrix_group* sequence = (matrix_group*)malloc(sizeof(matrix_group) * g_nthreads);
	for (int i = 0; i < g_nthreads; i++) {
		sequence[i].value = value;
		sequence[i].result = result + i * work_per_thread;
		pthread_create(&id[i], NULL, para_uniform, (void*)(sequence + i));
	}
	for (int i = 0; i < g_nthreads; i++) {
		pthread_join(id[i], NULL);
	}

	int remaining = work_per_thread * g_nthreads;
	for (int i = remaining; i < g_elements; i++) {
		result[i] = value;
	}

	free(sequence);
	return result;
}

/**
* Returns new matrix with elements in sequence from given start and step
*/
float* sequence_matrix(float start, float step) {

	float* result = new_matrix();

	/*

	1 2
	1 1 => 3 4
	*/
	if (step != 0.0) {
		float increment = start;
		for (int i = 0; i < g_elements; i ++) {
			result[i] = increment;
			increment += step;
		}
	} else {
		return uniform_matrix(start);
	}
	return result;
}

////////////////////////////////
///     MATRIX OPERATIONS    ///
////////////////////////////////

/**
* Returns new matrix with elements cloned from given matrix.
*/
float* cloned(const float* matrix) {

	float* result = new_matrix();

	pthread_t id[g_nthreads];

	ssize_t work_per_thread = (int)floor((float)g_elements / (float)g_nthreads);
	matrix_group sequence[g_nthreads];
	for (int i = 0; i < g_nthreads; i++) {
		sequence[i].a = matrix + i * work_per_thread;
		sequence[i].result = result + i * work_per_thread;
		pthread_create(&id[i], NULL, para_clone, (void*)(&sequence[i]));
	}

	for (int i = 0; i < g_nthreads; i++) {
		pthread_join(id[i], NULL);
	}
	int remaining = work_per_thread * g_nthreads;
	if (remaining != g_elements) {
		memcpy(result + remaining, matrix + remaining, sizeof(float) * (g_elements - remaining));
	}
	return result;
}

/**
* Returns new matrix with elements in ascending order.
*/
float* sorted(const float* matrix) {

	// float* result = new_matrix();

	/*

	3 4    1 2
	2 1 => 3 4

	*/
	// memcpy(result, matrix, sizeof(float) * g_elements);
	float* result = cloned(matrix);
	pthread_t id[4];
	Merge_args sequence[4];
	int work_per_thread = g_elements / 4;
	int work_last_thread = g_elements - 3 * work_per_thread;
	for (int i = 0; i < 4; i++) {
		sequence[i].array = result + i * work_per_thread;
		sequence[i].length = (i == 3)? work_last_thread: work_per_thread;
		sequence[i].depth = 0;
		pthread_create(&id[i], NULL, para_qsort, (void*)(sequence + i));
	}
	for (int i = 0; i < 4; i++) {
		pthread_join(id[i], NULL);
	}
	float* merge2 = calloc(2 * work_per_thread, sizeof(float));
	Merge(result, result + work_per_thread, merge2, work_per_thread, work_per_thread);

	float* merge3 = calloc(3 * work_per_thread, sizeof(float));
	Merge(merge2, result + 2 * work_per_thread, merge3, 2 * work_per_thread, work_per_thread);

	float* sort = calloc(g_elements, sizeof(float));
	Merge(merge3, result + 3 * work_per_thread, sort, 3 * work_per_thread, work_last_thread);
	free(merge2);
	free(merge3);
	free(result);
	return sort;
}


/**
* Returns new matrix with elements rotated 90 degrees clockwise.
*/
float* rotated(const float* matrix) {

	float* result = transposed(matrix);

	/*
	1 2    3 1
	3 4 => 4 2
	*/
	for (int i = 0; i < g_height; i++) {
		for (int j = 0; j < g_width / 2; j++) {
			float temp = result[i * g_width + j];
			result[i * g_width + j] = result[i * g_width + g_width - j - 1];
			result[i * g_width + g_width - j - 1] = temp;
		}
	}
	return result;
}

/**
* Returns new matrix with elements ordered in reverse.
*/
float* reversed(const float* matrix) {

	float* result = new_matrix();

	/*


	1 2    4 3
	3 4 => 2 1
	*/
	pthread_t id[g_nthreads];

	ssize_t work_per_thread = (int)floor((float)g_elements / (float)g_nthreads);
	matrix_group sequence[g_nthreads];
	for (int i = 0; i < g_nthreads; i++) {
		sequence[i].a = matrix;
		sequence[i].result = result;
		sequence[i].value = i * work_per_thread;
		pthread_create(&id[i], NULL, para_reversed, (void*)(&sequence[i]));
	}

	for (int i = 0; i < g_nthreads; i++) {
		pthread_join(id[i], NULL);
	}
	int remaining = work_per_thread * g_nthreads;
	for (int i = remaining; i < g_elements; i++) {
		result[i] = matrix[g_elements - 1 - i];
	}
	return result;
}


/**
* Returns new transposed matrix.
*/
float* transposed(const float* matrix) {

	float* result = new_matrix();
	/*


	1 2    1 3
	3 4 => 2 4
	*/
	int threads_available = (g_width < 4)? 1: g_nthreads;
	pthread_t id[threads_available];
	ssize_t work_per_thread = g_width / (int)sqrt(threads_available);
	submatrix* args = malloc(g_elements / threads_available * sizeof(submatrix));
	int num_args = 0;
	for (int i = 0; i < (int)sqrt(threads_available) ; i++) {
		for (int j = 0; j < (int)sqrt(threads_available) ; j++) {
			args[num_args].matrix_a = matrix;
			args[num_args].result = result;
			args[num_args].start_i = i * work_per_thread;
			args[num_args].start_j = j * work_per_thread;
			pthread_create(&id[num_args], NULL, para_transposed, (void*)(args + num_args));
			num_args++;
		}
	}


	num_args = 0;
	for (int i = 0; i < (int)sqrt(threads_available); i++) {
		for (int j = 0; j < (int)sqrt(threads_available); j++) {
			pthread_join(id[num_args], NULL);
			num_args++;
		}
	}
free(args);
	int remaining_index = (int)sqrt(threads_available) * work_per_thread;
	for (int j = 0; j < g_width; j++) {
		for (int i = remaining_index; i < g_width; i++) {
			result[IDX(j, i)] = matrix[IDX(i, j)];
		}
	}

	for (int j = remaining_index; j < g_width; j++) {
		for (int i = 0; i < g_width; i++) {
			result[IDX(j, i)] = matrix[IDX(i, j)];
		}
	}
	if (remaining_index != g_width) {
		result[IDX(g_height - 1, g_width - 1)] /= 2;
	}
	return result;
}


/**
* Returns new matrix with scalar added to each element.
*/
float* scalar_add(const float* matrix, float scalar) {

	float* result = new_matrix();

	/*
	1 0        2 1
	0 1 + 1 => 1 2

	1 2        5 6
	3 4 + 4 => 7 8
	*/
	if (scalar == 0) {
		return cloned(matrix);
	}
	__m128 oldvector, D, newvector;
	float d[4] = {scalar, scalar, scalar, scalar};
	D = _mm_load_ps(&d[0]);
	int k = 0;
	for (k = 0; k < g_elements; k += 4) {
		oldvector = _mm_load_ps(&matrix[k]);
		newvector = _mm_add_ps(oldvector, D);
		_mm_store_ps(&result[k], newvector);
	}
	for (int j = k - 4; j < g_elements % 4; j++) {
		result[j] = matrix[j] + scalar;
	}
	return result;
}


/**
* Returns new matrix with scalar multiplied to each element.
*/
float* scalar_mul(const float* matrix, float scalar) {

	float* result = new_matrix();

	/*


	1 0        2 0
	0 1 x 2 => 0 2

	1 2        2 4
	3 4 x 2 => 6 8
	*/
	if (scalar == 0) {
		return result;
	}

	__m128 oldvector, D, newvector;
	float d[4] = {scalar, scalar, scalar, scalar};
	D = _mm_load_ps(&d[0]);
	int k = 0;
	for (k = 0; k < g_elements; k += 4) {
		oldvector = _mm_load_ps(&matrix[k]);
		newvector = _mm_mul_ps(oldvector, D);
		_mm_store_ps(&result[k], newvector);
	}
	for (int j = k - 4; j < g_elements % 4; j++) {
		result[j] = matrix[j] * scalar;
	}
	return result;
}


/**
* Returns new matrix that is the result of
* adding the two given matrices together.
*/
float* matrix_add(const float* matrix_a, const float* matrix_b) {

	float* result = new_matrix();

	/*


	1 0   0 1    1 1
	0 1 + 1 0 => 1 1

	1 2   4 4    5 6
	3 4 + 4 4 => 7 8
	*/

	ssize_t work_per_thread = g_elements / g_nthreads;

	pthread_t* id = (pthread_t*)malloc(sizeof(pthread_t) * g_nthreads);

	matrix_group* sequence = (matrix_group*)malloc(sizeof(matrix_group) * g_nthreads);
	for (int i = 0; i < g_nthreads; i++) {
		sequence[i].a = matrix_a + i * work_per_thread;
		sequence[i].b = matrix_b + i * work_per_thread;//* sizeof(float) ??
		sequence[i].result = result + i * work_per_thread;
		pthread_create(&id[i], NULL, para_matrix_add, (void*)(sequence + i));
	}
	for (int i = 0; i < g_nthreads; i++) {
		pthread_join(id[i], NULL);
	}

	int remaining = work_per_thread * g_nthreads;
	for (int i = remaining; i < g_elements; i++) {
		result[i] = matrix_a[i] + matrix_b[i];
	}

	free(sequence);
	free(id);
	// free(group);
	return result;
}


/**
* Returns new matrix that is the result of
* multiplying the two matrices together.
*/
float* matrix_mul(const float* matrix_a, const float* matrix_b) {

	float* result = new_matrix();

	/*

	1 2   1 0    1 2
	3 4 x 0 1 => 3 4

	1 2   5 6    19 22
	3 4 x 7 8 => 43 50
	*/
	float* transposed_b = transposed(matrix_b);
	int threads_available = (g_width < 4)? 1: g_nthreads;
	pthread_t id[threads_available];
	ssize_t work_per_thread = g_elements / threads_available;
	matrix_group sequence[threads_available];
	for (int i = 0; i < threads_available; i++) {
		sequence[i].a = matrix_a;
		sequence[i].b = transposed_b;
		sequence[i].result = result;
		sequence[i].value = i * work_per_thread;
		pthread_create(&id[i], NULL, para_matrix_mul, (void*)(sequence + i));
	}
	for (int i = 0; i < threads_available; i++) {
		pthread_join(id[i], NULL);
	}
	int remaining_index = threads_available * work_per_thread;
	for (int i = remaining_index; i < g_elements; i++) {
		int row = floor(i / g_width);
		int column = i % g_width;
		register float element = 0;
		__m128 hvector, vvector, mul, D;
		float d[4] = {0};
		int k = 0;
		for (k = 0; k < g_width; k += 4) {
			hvector = _mm_load_ps(&matrix_a[IDX(row, k)]);
	        vvector = _mm_load_ps(&matrix_b[IDX(column, k)]);
	        mul = _mm_mul_ps(hvector, vvector);
	        D = _mm_load_ps(&d[0]);
	        D = _mm_add_ps(mul, D);
	        _mm_store_ps(&d[0], D);
		}
		for (int j = 0; j < 4; j++) {
			element += d[j];
		}
		for (int j = k - 4; j < work_per_thread % 4; j++) {
			element += matrix_a[IDX(row,j)] * matrix_b[IDX(column,j)];
		}
		result[i] = element;
	}

	free(transposed_b);
	return result;
}

/**
* Returns new matrix that is the result of
* powering the given matrix to the exponent.
*/
float* matrix_pow(const float* matrix, float exponent) {

	// float* result = new_matrix();

	/*


	1 2        1 0
	3 4 ^ 0 => 0 1

	1 2        1 2
	3 4 ^ 1 => 3 4

	1 2        199 290
	3 4 ^ 4 => 435 634
	*/
	if (exponent == 0) {
		return identity_matrix();
	} else if (exponent == 1) {
		float* result = cloned(matrix);
		return result;
	} else {
		return recursive_power((float*)matrix, (int)exponent);
	}
}

//Algorithm adapted from :
//http://cs.stackexchange.com/questions/4998/matrix-powering-in-o-log-n-time
//With some modifications
float* recursive_power(float* matrix, int exponent) {
	if (exponent == 0) {
		free(matrix);
		return identity_matrix();
	} else if (exponent == 1) {
		float* result = cloned(matrix);
		return result;
	} else {
		float* P = recursive_power(matrix, floor(exponent / 2));
		if (exponent % 2 == 0) {
			float* result = matrix_mul(P, P);
			free(P);
			return result;
		} else {
			float* temp = matrix_mul(P, P);
			float* result = matrix_mul(temp, matrix);
			free(temp);
			free(P);
			return result;
		}
	}
}

/**
* Returns new matrix that is the result of
* convolving given matrix with a 3x3 kernel matrix.
*/
float* matrix_conv(const float* matrix, const float* kernel) {

	float* result = new_matrix();

	/*
	TODO

	Convolution is the process in which the values of a matrix are
	computed according to the weighted sum of each value and it's
	neighbours, where the weights are given by the kernel matrix.
	*/
	for (ssize_t i = 0; i < g_width; i++) {
		for (ssize_t j = 0; j < g_width; j++) {
			double element_result = 0;
			int cnt = 0;
			for (int y = i - 1; y < i + 2; y++) {

				for (int x = j - 1; x < j + 2; x++) {
					if ((x < 0 || y < 0) && x < g_width && y < g_height) {
						//upper left corner
						if (x < 0 && y < 0) {
							element_result += matrix[IDX(y + 1, x + 1)] * kernel[cnt];
							//left
						} else if (x < 0 && y >= 0) {
							element_result += matrix[IDX(y, x + 1)] * kernel[cnt];
							//up
						} else if (x >= 0 && y < 0){
							element_result += matrix[IDX(y + 1, x)] * kernel[cnt];
						}
					} else if ((x >= g_width || y >= g_height) && x >= 0 && y >= 0) {
						//bottom right coner
						if (x >= g_width && y >= g_height) {
							element_result += matrix[IDX(y - 1, x - 1)] * kernel[cnt];
							//right
						} else if (x >= g_width && y < g_height) {
							element_result += matrix[IDX(y, x - 1)] * kernel[cnt];
							//bottom
						} else if (x < g_width && y >= g_height) {
							element_result += matrix[IDX(y - 1, x)] * kernel[cnt];
						}
						//bottom left
						// printf("bottom left: %d, %d\n", y, x);
					} else if (x < 0 && y >= g_height) {
						element_result += matrix[IDX(y - 1, x + 1)] * kernel[cnt];
						//upper right coner
					} else if (x >= g_width && y < 0) {
						element_result += matrix[IDX(y + 1, x - 1)] * kernel[cnt];
						//rest
					} else {
						element_result += matrix[IDX(y, x)] * kernel[cnt];
					}
					cnt++;
				}
			}

			result[IDX(i, j)] = element_result;
		}
	}
	return result;
}

////////////////////////////////
///       COMPUTATIONS       ///
////////////////////////////////

/**
* Returns the sum of all elements.
*/
float get_sum(const float* matrix) {

	/*

	2 1
	1 2 => 6

	1 1
	1 1 => 4
	*/
	pthread_t id[g_nthreads];

	ssize_t work_per_thread = (int)floor((float)g_elements / (float)g_nthreads);
	int* sequence = malloc(sizeof(int) * g_nthreads);
	for (int i = 0; i < g_nthreads; i++) {
		sequence[i] = i;
		pthread_create(&id[i], NULL, para_sum, (void*)(matrix + sequence[i] * work_per_thread));
	}
	void** result = malloc(sizeof(float*));
	float sum = 0;
	for (int i = 0; i < g_nthreads; i++) {
		//Segment fault
		pthread_join(id[i], result);
		float returned = *((float*)*result);
		sum = sum + returned;
		free(*result);
	}
	free(result);


	int remaining = work_per_thread * g_nthreads;
	for (int i = remaining; i < g_elements; i++) {
		sum += matrix[i];
	}

	free(sequence);
	return sum;
}


/**
* Returns the trace of the matrix.
*/
float get_trace(const float* matrix) {

	/*
	1 0
	0 1 => 2

	2 1
	1 2 => 4
	*/
	//Parralleling it makes it fail all tests
	//BUGGY!!
	pthread_t id[g_nthreads];

	ssize_t work_per_thread = g_width / g_nthreads;
	matrix_group* sequence = (matrix_group*)malloc(sizeof(matrix_group) * g_nthreads);
	for (int i = 0; i < g_nthreads; i++) {
		sequence[i].a = matrix;
		sequence[i].value = i * work_per_thread;
		pthread_create(&id[i], NULL, para_trace, (void*)(sequence + i));
	}

	void** result = malloc(sizeof(double*));
	float trace = 0;
	for (int i = 0; i < g_nthreads; i++) {
		pthread_join(id[i], result);
		trace += *((float*)*result);
		free(((void*)*result));
	}
	free(result);

	int remaining = work_per_thread * g_nthreads;
	for (int i = remaining; i < g_width; i++) {
		trace += matrix[IDX(i, i)];
	}

	free(sequence);
	return trace;
}


/**
* Returns the smallest value in the matrix.
*/
float get_minimum(const float* matrix) {

	/*
	1 2
	3 4 => 1

	4 3
	2 1 => 1
	*/
	pthread_t id[g_nthreads];
	float min = matrix[0];
	ssize_t work_per_thread = g_elements / g_nthreads;
	int* sequence = malloc(sizeof(int) * g_nthreads);
	for (int i = 0; i < g_nthreads; i++) {
		sequence[i] = i;
		pthread_create(&id[i], NULL, para_min, (void*)(matrix + sequence[i] * work_per_thread));
	}
	void** result = malloc(sizeof(float*));
	for (int i = 0; i < g_nthreads; i++) {
		pthread_join(id[i], result);
		if (min > *((float*)*result)) {
			min = *((float*)*result);
		}
		free(((void*)*result));
	}
	free(result);
	int remaining = work_per_thread * g_nthreads;
	for (int i = remaining; i < g_elements; i++) {
		if (min > matrix[i]) {
			min = matrix[i];
		}
	}

	free(sequence);

	return min;
}

/**
* Returns the largest value in the matrix.
*/
float get_maximum(const float* matrix) {

	/*

	1 2
	3 4 => 4

	4 3
	2 1 => 4
	*/

	pthread_t id[g_nthreads];
	float max = matrix[0];
	ssize_t work_per_thread = g_elements / g_nthreads;
	int* sequence = malloc(sizeof(int) * g_nthreads);
	for (int i = 0; i < g_nthreads; i++) {
		sequence[i] = i;
		pthread_create(&id[i], NULL, para_max, (void*)(matrix + sequence[i] * work_per_thread));
	}
	void** result = malloc(sizeof(float*));
	for (int i = 0; i < g_nthreads; i++) {
		pthread_join(id[i], result);
		if (max < *((float*)*result)) {
			max = *((float*)*result);
		}
		free(((void*)*result));
	}
	free(result);
	int remaining = work_per_thread * g_nthreads;
	for (int i = remaining; i < g_elements; i++) {
		if (max < matrix[i]) {
			max = matrix[i];
		}
	}

	free(sequence);

	return max;
}

/**
* Returns the determinant of the matrix.
*/
float get_determinant(const float* matrix) {

	/*
	TODO
	Permutaion matrix???
	1 0
	0 1 => 1

	1 2
	3 4 => -2

	8 0 2
	0 4 0
	2 0 8 => 240
	*/
	return laplace(matrix, g_width);
}


/**
* Returns the frequency of the given value in the matrix.
*/
ssize_t get_frequency(const float* matrix, float value) {

	/*

	1 1
	1 1 :: 1 => 4

	1 0
	0 1 :: 2 => 0
	*/

	ssize_t work_per_thread = g_elements / g_nthreads;
	pthread_t id[g_nthreads];
	matrix_group sequence[g_nthreads];
	for (int i = 0; i < g_nthreads; i++) {
		sequence[i].a = matrix + i * work_per_thread;
		sequence[i].value = value;
		pthread_create(&id[i], NULL, para_freq, (void*)(sequence + i));
	}

	void** result = malloc(sizeof(ssize_t*));
	ssize_t freq = 0;

	for (int i = 0; i < g_nthreads; i++) {
		pthread_join(id[i], result);
		freq += *((ssize_t*)*result);
		free(((void*)*result));
	}
	free(result);
	int remaining = work_per_thread * g_nthreads;
	for (int i = remaining; i < g_elements; i++) {
		if (matrix[i] == value) {
			freq ++;
		}
	}
	return freq;
}
