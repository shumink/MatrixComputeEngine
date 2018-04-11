#ifndef MATRIX_H
#define MATRIX_H

#include <unistd.h>

// utility functions

int fast_rand(void);

void set_seed(int value);
void set_nthreads(ssize_t count);
void set_dimensions(ssize_t width);

void display(const float* matrix);
void display_row(const float* matrix, ssize_t row);
void display_column(const float* matrix, ssize_t column);
void display_element(const float* matrix, ssize_t row, ssize_t column);

// matrix operations

float* new_matrix(void);
float* identity_matrix(void);

float* random_matrix(int seed);
float* uniform_matrix(float value);
float* sequence_matrix(float start, float step);

float* cloned(const float* matrix);
float* sorted(const float* matrix);
float* rotated(const float* matrix);
float* reversed(const float* matrix);
float* transposed(const float* matrix);

float* scalar_add(const float* matrix, float scalar);
float* scalar_mul(const float* matrix, float scalar);
float* matrix_pow(const float* matrix, float exponent);
float* matrix_conv(const float* matrix, const float* kernel);
float* matrix_add(const float* matrix_a, const float* matrix_b);
float* matrix_mul(const float* matrix_a, const float* matrix_b);

// compute operations

float get_sum(const float* matrix);
float get_trace(const float* matrix);
float get_minimum(const float* matrix);
float get_maximum(const float* matrix);
float get_determinant(const float* matrix);
ssize_t get_frequency(const float* matrix, float value);

#endif
