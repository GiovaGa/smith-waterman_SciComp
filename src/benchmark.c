#include<stdlib.h>
#include<stdio.h>
#include<omp.h>
#include"benchmark.h"

void benchmark(struct sw_implementation* implementation, int n_implementations, struct sequence_t A, struct sequence_t B, const struct scores_t *scores_param)
{
	int* H = (int*)malloc((A.length + 1) * (B.length + 1) * sizeof(int));
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
		memset(H, 0, (A.length + 1) * (B.length + 1) * sizeof(int));
		start_time = omp_get_wtime();
		score = implementation[i].function(&A, &B, scores_param, H);
		elapsed_time = omp_get_wtime() - start_time;

		gflops = implementation[i].flops(A.length, B.length) / (elapsed_time * 1e9);
		printf("%-20s %-10d %-15f %-15f\n", implementation[i].name, score, elapsed_time, gflops);
	}
	free(H);
}
