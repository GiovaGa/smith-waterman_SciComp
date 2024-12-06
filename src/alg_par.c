#include<omp.h>

static inline int max(const int a, const int b) {
	if (a >= b) return a;
	else return b;
}

int smith_waterman_parallel(const int N, const int M, const char* restrict A, const char* restrict B, int score_match, int score_open_gap, int score_continue_gap, int score_mismatch, int* restrict H) {

	int ans = 0;
	// int* W; // Gap weight

	for (int j = 1; j <= M; ++j) {
#pragma omp parallel for reduction (max:ans) firstprivate(j)
	    for (int i = 1; i <= N; ++i) {
			const int score = score_match * (A[i-1] == B[j-1]) + score_mismatch * (A[i-1] != B[j-1]);
			H[i * (M + 1) + j] = max(H[i * (M + 1) + j], H[(i - 1) * (M + 1) + j - 1] + score);
			for (int k = 1; k < i; ++k) {
				H[i * (M + 1) + j] = max(H[i * (M + 1) + j], H[(i - k) * (M + 1) + j] + score_open_gap + score_continue_gap * (k-1));
			}
			for (int k = 1; k < j; ++k) {
				H[i * (M + 1) + j] = max(H[i * (M + 1) + j], H[i * (M + 1) + j - k] + score_open_gap + score_continue_gap * (k-1));
			}
			ans = max(ans, H[i * (M + 1) + j]);
		}
	}
	return ans;
}

