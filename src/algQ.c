#include<stdlib.h>
#include<string.h>
#include<omp.h>

static inline int max(const int a, const int b) {
	if (a >= b) return a;
	else return b;
}

int smith_waterman_quadratic(const int N, const int M, const char* restrict A, const char* restrict B, int score_match, int score_open_gap, int score_continue_gap, int score_mismatch, int* restrict H) {

	int ans = 0;
    // int* W; // Gap weight
	int *Mi = malloc((N+1)*sizeof(int)); memset(Mi,score_open_gap,(N+1)*sizeof(int));
    int *Mj = malloc((M+1)*sizeof(int)); memset(Mj,score_open_gap,(M+1)*sizeof(int));

	for (int j = 1; j <= M; ++j) {
	    for (int i = 1; i <= N; ++i) {
			const int score = score_match * (A[i-1] == B[j-1]) + score_mismatch * (A[i-1] != B[j-1]);
			H[i * (M + 1) + j] = max(H[i * (M + 1) + j], H[(i - 1) * (M + 1) + j - 1] + score);

			H[i * (M + 1) + j] = max(H[i * (M + 1) + j], Mj[j]);
			H[i * (M + 1) + j] = max(H[i * (M + 1) + j], Mi[i]);

			ans = max(ans, H[i * (M + 1) + j]);

            Mj[j] = max(Mj[j] + score_continue_gap, H[i * (M + 1) + j] + score_open_gap);
            Mi[i] = max(Mi[i] + score_continue_gap, H[i * (M + 1) + j] + score_open_gap);
		}
	}
	return ans;
}

long smith_waterman_flops_quadratic(int N, int M)
{
	return (long)10*(N)*M;
}
