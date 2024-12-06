#include<stdlib.h>
#include<string.h>
#include<omp.h>
#include"benchmark.h"

static inline int max(const int a, const int b) {
	if (a >= b) return a;
	else return b;
}

int smith_waterman_quadratic(const int N, const int M, const char* restrict A, const char* restrict B, const struct scores_t*scores, int* restrict H) {
	int ans = 0;
    // int* W; // Gap weight
	int *Mi = malloc((N+1)*sizeof(int));
	for(int k = 0; k < (N+1); ++k)
		Mi[k] = scores->gap_opening;
    int *Mj = malloc((M+1)*sizeof(int));
	for(int k = 0; k < (M+1); ++k)
		Mj[k] = scores->gap_opening;

	for (int j = 1; j <= M; ++j) {
	    for (int i = 1; i <= N; ++i) {
			const int score = scores->match * (A[i-1] == B[j-1]) + scores->mismatch * (A[i-1] != B[j-1]);
			H[i * (M + 1) + j] = max(H[i * (M + 1) + j], H[(i - 1) * (M + 1) + j - 1] + score);

			H[i * (M + 1) + j] = max(H[i * (M + 1) + j], Mj[j]);
			H[i * (M + 1) + j] = max(H[i * (M + 1) + j], Mi[i]);

			ans = max(ans, H[i * (M + 1) + j]);

            Mj[j] = max(Mj[j] + scores->gap_extension, H[i * (M + 1) + j] + scores->gap_opening);
            Mi[i] = max(Mi[i] + scores->gap_extension, H[i * (M + 1) + j] + scores->gap_opening);
		}
	}
	
	free(Mi);
	free(Mj);
	return ans;
}

long smith_waterman_flops_quadratic(int N, int M)
{
	return (long)10*(N)*M;
}
