#include<stdlib.h>
#include<string.h>
#include<omp.h>
#include"benchmark.h"

static inline int max(const int a, const int b) {
	if (a > b) return a;
	return b;
}

int smith_waterman_quadratic_opt(const int N, const int M, const char* restrict A, const char* restrict B, const struct scores_t*scores, int* restrict H) {

	//if(M < N) return smith_waterman_quadratic_opt(M, N, B, A, scores->match, scores->gap_opening, scores->gap_extension, scores->mismatch, H);
	int ans = 0;
    int *Mj = malloc((M+1)*sizeof(int));
	for(int k = 0; k < (M+1); ++k)
		Mj[k] = scores->gap_opening;

	for (int i = 1; i <= N; ++i)
	{
		int Mi = scores->gap_opening;
		const char a = A[i - 1];
	    for (int j = 1; j <= M; ++j)
		{
			const int score = scores->match * (a == B[j-1]) + scores->mismatch * (a != B[j-1]);
			
			int h = max(0, H[(i - 1) * (M + 1) + j - 1] + score);
			h = max(h, Mj[j]);
			h = max(h, Mi);

			ans = max(ans, h);

            Mj[j] = max(Mj[j] + scores->gap_extension, h + scores->gap_opening);
            Mi = max(Mi + scores->gap_extension, h + scores->gap_opening);
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
