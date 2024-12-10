#ifndef BENCHMARK_H
#define BENCHMARK_H

#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "sequence.h"


struct sw_implementation
{
	const char* name;
	int (*function)(const struct sequence_t*, const struct sequence_t*, const struct scores_t*);
	size_t (*flops)(size_t, size_t);
};


void benchmark(struct sw_implementation* implementation, int n_implementations, int n_runs,  struct sequence_t A, struct sequence_t B, const struct scores_t *scores_param);
double avg(double* array, size_t array_len);
double std_dev(double* array, size_t array_len, double avg);
int are_scores_equal(int* array, size_t array_len);

#endif
