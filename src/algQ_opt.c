#include<stdlib.h>
#include<string.h>
#include<omp.h>

static inline int max(const int a, const int b) {
	if (a > b) return a;
	return b;
}

int smith_waterman_quadratic_opt(const int N, const int M, const char* restrict A, const char* restrict B, int score_match, int score_open_gap, int score_continue_gap, int score_mismatch, int* restrict H) {

	//if(M < N) return smith_waterman_quadratic_opt(M, N, B, A, score_match, score_open_gap, score_continue_gap, score_mismatch, H);
	int ans = 0;
    int *Mj = malloc((M+1)*sizeof(int));
	for(int k = 0; k < (M+1); ++k)
		Mj[k] = score_open_gap;

	for (int i = 1; i <= N; ++i)
	{
		int Mi = score_open_gap;
		const char a = A[i - 1];
	    for (int j = 1; j <= M; ++j)
		{
			const int score = score_match * (a == B[j-1]) + score_mismatch * (a != B[j-1]);
			
			const int h = max(0, H[(i - 1) * (M + 1) + j - 1] + score);
			h = max(h, Mj[j]);
			h = max(h, Mi);

			ans = max(ans, h);

            Mj[j] = max(Mj[j] + score_continue_gap, h + score_open_gap);
            Mi = max(Mi + score_continue_gap, h + score_open_gap);
			H[i * (M + 1) + j] = h;
		}
	}
	free(Mj);
	return ans;
}

long smith_waterman_flops_quadratic_opt(int N, int M)
{
	return (long)10*(N)*M;
}
