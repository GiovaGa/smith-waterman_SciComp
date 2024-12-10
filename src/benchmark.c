#include<stdlib.h>
#include<stdio.h>
#include<omp.h>
#include<math.h>

#include"benchmark.h"

double avg(double* array, size_t array_len) {
	double acc = 0;
	for (size_t i=0; i<array_len; i++)
		acc += array[i];
	return acc/array_len;
}
double std_dev(double* array, size_t array_len, double avg) {
	double acc = 0;
	for (size_t i=0; i<array_len; i++)
		acc += pow(array[i] - avg, 2);
	return sqrt(acc/array_len);
}

int are_scores_equal(int* array, size_t array_len) {
	for (size_t i=0; i<array_len-1; i++) {
		if (array[i] != array[i+1]) return 0;
	}
	return 1;
}

void benchmark(struct sw_implementation* implementation, int n_implementations, int n_runs, struct sequence_t A, struct sequence_t B, const struct scores_t *scores_param)
{
	double start_time, elapsed_time, gflops;
	int score = 0;
	int* score_container = (int*)malloc((n_runs) * sizeof(int));
	double* time_container = (double*)malloc((n_runs) * sizeof(double));
	double* gflops_container = (double*)malloc((n_runs) * sizeof(double));

	printf("%-30s %-10s %-10s %-22s %-10s %-27s %s\n", "Function", "Score", "Consistent", "AVG Elapsed Time (s)", "STD DEV", "AVG Performance (GFLOPS/s)", "STD DEV");
	printf("-----------------------------------------------------------------------------------------------------------------\n");

	for (int i = 0; i < n_implementations; i++) {
		for (int j=0; j<n_runs; j++) {
			start_time = omp_get_wtime();
			score = implementation[i].function(&A, &B, scores_param);
			score_container[j] = score;

			elapsed_time = omp_get_wtime() - start_time;
			time_container[j] = elapsed_time;

			gflops = implementation[i].flops(A.length, B.length) / (elapsed_time * 1e9);
			gflops_container[j] = gflops;
		}
		double time_avg = avg(time_container, n_runs);
		double time_stdev = std_dev(time_container, n_runs, time_avg);
		char *is_score_same = are_scores_equal(score_container, n_runs) ? "Yes" : "No";

		double gflops_avg = avg(gflops_container, n_runs);
		double gflops_stdev = std_dev(gflops_container, n_runs, gflops_avg);
		
		printf("%-30s %-10d %-10s %-22f %-10f %-27f %f\n", implementation[i].name, score, is_score_same, time_avg, time_stdev, gflops_avg, gflops_stdev);
		
		memset(time_container, 0, n_runs * sizeof(double));
		memset(gflops_container, 0, n_runs * sizeof(double));
	}
}
