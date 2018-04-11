/*
######################################
### DO NOT MODIFY THIS SOURCE FILE ###
######################################
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <stdbool.h>

#include "matrix.h"

#define MAX_BUFFER 256
#define MAX_ENTRIES 512

#define MATRIX_GUARD(x) \
	float* m = find_matrix(x); \
	if (m == NULL) { \
		puts("no such matrix"); \
		return; \
	}

#define MATRIX_GUARD_PAIR(x, y) \
	float* m1 = find_matrix(x); \
	float* m2 = find_matrix(y); \
	if (m1 == NULL || m2 == NULL) { \
		puts("no such matrix"); \
		return; \
	}

typedef struct entry {
	char key[MAX_BUFFER];
	float* matrix;
} entry;

static ssize_t g_order	= 0;
static ssize_t g_nthreads = 0;
static ssize_t g_nentries = 0;

static entry** g_entries = NULL;

/**
 * Adds entry to list of entries.
 */
entry* add_entry(char* key) {

	entry* e = calloc(1, sizeof(entry));
	strcpy(e->key, key);
	g_entries[g_nentries] = e;
	g_nentries += 1;

	return e;
}

/**
 * Returns entry with given key.
 */
entry* find_entry(char* key) {

	for (ssize_t i = 0; i < g_nentries; i++) {
		if (strcmp(key, g_entries[i]->key) == 0) {
			return g_entries[i];
		}
	}

	return NULL;
}

/**
 * Returns matrix with given key.
 */
float* find_matrix(char* key) {

	entry* e = find_entry(key);
	if (e == NULL) {
		return NULL;
	}

	return e->matrix;
}

/**
 * Returns matrix kernel with given name.
 */
const float* find_kernel(char* name) {

	static const float blur[9] = {
		0.0625, 0.1250, 0.0625,
		0.1250, 0.2500, 0.1250,
		0.0625, 0.1250, 0.0625,
	};

	static const float edge[9] = {
		 0,  1,  0,
		 1, -4,  1,
		 0,  1,  0,
	};

	static const float emboss[9] = {
		-2, -1,  0,
		-1,  1,  1,
		 0,  1,  2,
	};

	static const float sharpen[9] = {
	     0, -1,  0,
		-1,  5, -1,
		 0, -1,  0,
	};

	static const float outline[9] = {
		-1, -1, -1,
		-1,  8, -1,
		-1, -1, -1,
	};

	static const float identity[9] = {
		 0,  0,  0,
		 0,  1,  0,
		 0,  0,  0,
	};

	if (strcasecmp(name, "blur") == 0) {
		return blur;
	} else if (strcasecmp(name, "edge") == 0) {
		return edge;
	} else if (strcasecmp(name, "emboss") == 0) {
		return emboss;
	} else if (strcasecmp(name, "sharpen") == 0) {
		return sharpen;
	} else if (strcasecmp(name, "outline") == 0) {
		return outline;
	} else if (strcasecmp(name, "identity") == 0) {
		return identity;
	}

	return NULL;
}

/**
 * Releases all dynamically allocated memory used by the entries.
 */
void release(void) {

	if (g_entries == NULL) {
		return;
	}

	for (ssize_t i = 0; i < g_nentries; i++) {
		free(g_entries[i]->matrix);
		free(g_entries[i]);
	}

	free(g_entries);
}

/**
 * Initialises settings based on command line arguments.
 */
void init(int argc, char** argv) {

	if (argc != 3) {
		goto invalid;
	}

	g_order = atol(argv[1]);
	g_nthreads = atol(argv[2]);

	if (g_order < 1 || g_nthreads < 1) {
		goto invalid;
	}

	set_nthreads(g_nthreads);
	set_dimensions(g_order);
	return;

invalid:
	puts("Usage: matrix <width> <threads>");
	exit(1);
}

/**
 * Bye command.
 */
void command_bye(void) {

	printf("bye\n");
	release();
	exit(0);
}

/**
 * Help command.
 */
void command_help(void) {

	const char* HELP =
		"BYE\n"
		"HELP\n"
		"\n"
		"SET <key> = identity\n"
		"SET <key> = random <seed>\n"
		"SET <key> = uniform <value>\n"
		"SET <key> = sequence <start> <step>\n"
		"\n"
		"SET <key> = cloned <matrix>\n"
		"SET <key> = sorted <matrix>\n"
		"SET <key> = rotated <matrix>\n"
		"SET <key> = reversed <matrix>\n"
		"SET <key> = transposed <matrix>\n"
		"\n"
		"SET <key> = scalar.add <matrix> <scalar>\n"
		"SET <key> = scalar.mul <matrix> <scalar>\n"
		"SET <key> = matrix.add <matrix a> <matrix b>\n"
		"SET <key> = matrix.mul <matrix a> <matrix b>\n"
		"SET <key> = matrix.pow <matrix> <exponent>\n"
		"SET <key> = matrix.conv <matrix> <kernel>\n"
		"\n"
		"SHOW <key>\n"
		"SHOW <key> row <number>\n"
		"SHOW <key> column <number>\n"
		"SHOW <key> element <row> <column>\n"
		"\n"
		"COMPUTE sum <key>\n"
		"COMPUTE trace <key>\n"
		"COMPUTE minimum <key>\n"
		"COMPUTE maximum <key>\n"
		"COMPUTE determinant <key>\n"
		"COMPUTE frequency <key> <value>\n";

	printf("%s", HELP);
}

/**
 * Set command.
 */
void command_set(const char* line) {

	char cmd[MAX_BUFFER];
	char key[MAX_BUFFER];
	char func[MAX_BUFFER];
	char arg1[MAX_BUFFER];
	char arg2[MAX_BUFFER];

	int argc = sscanf(line, "%s %s = %s %s %s", cmd, key, func, arg1, arg2);
	if (argc < 3) {
		puts("invalid arguments");
		return;
	}

	float* matrix = NULL;

	switch (argc) {
		case 3:
			if (strcasecmp(func, "identity") == 0) {
				matrix = identity_matrix();
			} else {
				goto invalid;
			}
			break;

		case 4:
			if (strcasecmp(func, "random") == 0) {
				int seed = atol(arg1);
				matrix = random_matrix(seed);
			} else if (strcasecmp(func, "uniform") == 0) {
				float value = atof(arg1);
				matrix = uniform_matrix(value);
			} else if (strcasecmp(func, "cloned") == 0) {
				MATRIX_GUARD(arg1);
				matrix = cloned(m);
			} else if (strcasecmp(func, "sorted") == 0) {
				MATRIX_GUARD(arg1);
				matrix = sorted(m);
			} else if (strcasecmp(func, "rotated") == 0) {
				MATRIX_GUARD(arg1);
				matrix = rotated(m);
			} else if (strcasecmp(func, "reversed") == 0) {
				MATRIX_GUARD(arg1);
				matrix = reversed(m);
			} else if (strcasecmp(func, "transposed") == 0) {
				MATRIX_GUARD(arg1);
				matrix = transposed(m);
			} else {
				goto invalid;
			}
			break;

		case 5:
			if (strcasecmp(func, "sequence") == 0) {
				float start = atof(arg1);
				float step = atof(arg2);
				matrix = sequence_matrix(start, step);
			} else if (strcasecmp(func, "scalar.add") == 0) {
				MATRIX_GUARD(arg1);
				float value = atof(arg2);
				matrix = scalar_add(m, value);
			} else if (strcasecmp(func, "scalar.mul") == 0) {
				MATRIX_GUARD(arg1);
				float value = atof(arg2);
				matrix = scalar_mul(m, value);
			} else if (strcasecmp(func, "matrix.add") == 0) {
				MATRIX_GUARD_PAIR(arg1, arg2);
				matrix = matrix_add(m1, m2);
			} else if (strcasecmp(func, "matrix.mul") == 0) {
				MATRIX_GUARD_PAIR(arg1, arg2);
				matrix = matrix_mul(m1, m2);
			} else if (strcasecmp(func, "matrix.pow") == 0) {
				MATRIX_GUARD(arg1);
				float exponent = atof(arg2);
				matrix = matrix_pow(m, exponent);
			} else if (strcasecmp(func, "matrix.conv") == 0) {
				MATRIX_GUARD(arg1);
				const float* kernel = find_kernel(arg2);
				if (kernel == NULL) {
					puts("no such kernel");
				}
				matrix = matrix_conv(m, kernel);
			} else {
				goto invalid;
			}
			break;
	}

	entry* e = find_entry(key);
	if (e == NULL) {
		e = add_entry(key);
	} else {
		free(e->matrix);
	}

	e->matrix = matrix;

	puts("ok");
	return;

invalid:
	puts("invalid arguments");
}

/**
 * Show command.
 */
void command_show(char* line) {

	char cmd[MAX_BUFFER];
	char key[MAX_BUFFER];
	char func[MAX_BUFFER];
	char arg1[MAX_BUFFER];
	char arg2[MAX_BUFFER];

	int argc = sscanf(line, "%s %s %s %s %s", cmd, key, func, arg1, arg2);
	if (argc < 2) {
		goto invalid;
	}

	MATRIX_GUARD(key);
	if (argc == 2) {
		display(m);
		return;
	}

	const float v1 = atof(arg1) - 1;
	if (v1 >= g_order) {
		goto invalid;
	}

	if (argc == 4 && strcasecmp(func, "row") == 0) {
		display_row(m, v1);
	} else if (argc == 4 && strcasecmp(func, "column") == 0) {
		display_column(m, v1);
	} else if (argc == 5 && strcasecmp(func, "element") == 0) {
		const float v2 = atof(arg2) - 1;
		if (v2 >= g_order) {
			goto invalid;
		}
		display_element(m, v1, v2);
	}

	return;

invalid:
	puts("invalid arguments");
}

/**
 * Compute command.
 */
void command_compute(char* line) {

	char cmd[MAX_BUFFER];
	char key[MAX_BUFFER];
	char func[MAX_BUFFER];
	char arg1[MAX_BUFFER];

	int argc = sscanf(line, "%s %s %s %s", cmd, func, key, arg1);
	if (argc < 3) {
		goto invalid;
	}

	MATRIX_GUARD(key);
	float result = 0;

	if (strcasecmp(func, "sum") == 0) {
		result = get_sum(m);
	} else if (strcasecmp(func, "trace") == 0) {
		result = get_trace(m);
	} else if (strcasecmp(func, "minimum") == 0) {
		result = get_minimum(m);
	} else if (strcasecmp(func, "maximum") == 0) {
		result = get_maximum(m);
	} else if (strcasecmp(func, "determinant") == 0) {
		result = get_determinant(m);
	} else if (strcasecmp(func, "frequency") == 0) {
		ssize_t count = get_frequency(m, atof(arg1));
		printf("%zu\n", count);
		return;
	} else {
		goto invalid;
	}

	printf("%.2f\n", result);
	return;

invalid:
	puts("invalid arguments");
}

/**
 * Runs computations and stores matrices based on given input.
 */
void compute_engine(void) {

	g_entries = calloc(MAX_ENTRIES, sizeof(entry*));

	while (true) {
		printf("> ");

		char line[MAX_BUFFER];
		if (fgets(line, MAX_BUFFER, stdin) == NULL) {
			printf("\n");
			command_bye();
		}

		char command[MAX_BUFFER];
		if (sscanf(line, "%s", command) != 1) {
			printf("\n");
			continue;
		}

		if (strcasecmp(command, "bye") == 0) {
			command_bye();
		} else if (strcasecmp(command, "help") == 0) {
			command_help();
		} else if (strcasecmp(command, "set") == 0) {
			command_set(line);
		} else if (strcasecmp(command, "show") == 0) {
			command_show(line);
		} else if (strcasecmp(command, "compute") == 0) {
			command_compute(line);
		} else {
			puts("no such command");
		}

		printf("\n");
	}
}

/**
 * Main function.
 */
int main(int argc, char** argv)
{
	init(argc, argv);
	compute_engine();

	return 0;
}
