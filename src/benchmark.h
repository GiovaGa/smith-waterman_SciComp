#ifndef BENCHMARK_H
#define BENCHMARK_H

#include <omp.h>
#include <stdio.h>
#include <stdlib.h>

#include "sequence.h"

struct sw_implementation
{
	const char* name;
	int (*function)(int, int, const char*, const char*, int, int, int, int, int*);
	long (*flops)(int, int);
};

void benchmark(struct sw_implementation* implementation, int n_implementations, struct sequence_t A, struct sequence_t B)
{
	int* H = (int*)calloc((A.length + 1) * (B.length + 1), sizeof(int));
	if (!H)
	{
		printf("Cannot allocate memory for score matrix\n");
		return;
	}
	double start_time, elapsed_time, gflops;
	int score;

	printf("\n\n%-20s %-10s %-15s %-15s\n", "Function", "Score", "Elapsed Time (s)", "Performance (GFLOPS/s)");
	printf("-------------------------------------------------------------\n");

	for (int i = 0; i < n_implementations; i++)
	{
		start_time = omp_get_wtime();
		score = implementation[i].function(A.length, B.length, A.data, B.data, 4, -3, -1, -2, H); // TODO: allow for different scores
		elapsed_time = omp_get_wtime() - start_time;

		gflops = implementation[i].flops(A.length, B.length) / (elapsed_time * 1e9);
		printf("%-20s %-10d %-15f %-15f\n", implementation[i].name, score, elapsed_time, gflops);
	}
	free(H);
}


#endif
