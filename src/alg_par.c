#include<omp.h>
#include"benchmark.h"

static inline int max(const int a, const int b) {
	if (a >= b) return a;
	else return b;
}

int smith_waterman_parallel(const int N, const int M, const char* restrict A, const char* restrict B, const struct scores_t*scores, int* restrict H) {

	int ans = 0;
	// int* W; // Gap weight

	for (int j = 1; j <= M; ++j) {
#pragma omp parallel for reduction (max:ans) firstprivate(j)
	    for (int i = 1; i <= N; ++i) {
			const int score = scores->match * (A[i-1] == B[j-1]) + scores->mismatch * (A[i-1] != B[j-1]);
			H[i * (M + 1) + j] = max(H[i * (M + 1) + j], H[(i - 1) * (M + 1) + j - 1] + score);
			for (int k = 1; k < i; ++k) {
				H[i * (M + 1) + j] = max(H[i * (M + 1) + j], H[(i - k) * (M + 1) + j] + scores->gap_opening + scores->gap_extension * (k-1));
			}
			for (int k = 1; k < j; ++k) {
				H[i * (M + 1) + j] = max(H[i * (M + 1) + j], H[i * (M + 1) + j - k] + scores->gap_opening + scores->gap_extension * (k-1));
			}
			ans = max(ans, H[i * (M + 1) + j]);
		}
	}
	return ans;
}

