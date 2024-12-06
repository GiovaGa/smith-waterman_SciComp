#include<stdlib.h>
#include<string.h>
#include<omp.h>
#include <stdbool.h>

static inline int min(const int a, const int b) {
	if (a < b) return a;
	return b;
}

static inline int max(const int a, const int b) {
	if (a > b) return a;
	return b;
}

static inline void fill(int v, int* p, size_t n)
{
	for(size_t k = 0; k < n; ++k)
		p[k] = v;
}

#define I_BLOCK_SIZE 8
#define J_BLOCK_SIZE (64/sizeof(int))


int smith_waterman_quadratic_parallel(size_t N, size_t M, const char* restrict A, const char* restrict B, int score_match, int score_open_gap, int score_continue_gap, int score_mismatch, int* restrict H) {

	//if(M < N) return smith_waterman_quadratic_opt(M, N, B, A, score_match, score_open_gap, score_continue_gap, score_mismatch, H);
	int ans = 0;
    int *Mj = malloc((M+1)*sizeof(int));
	fill(score_open_gap, Mj, M+1);
	
	for (size_t ii = 0; ii <= N; ii += I_BLOCK_SIZE)
	{
		int Mi[I_BLOCK_SIZE];
		fill(score_open_gap, Mi, I_BLOCK_SIZE);
		
		for (size_t jj = 0; jj <= M; jj += J_BLOCK_SIZE)
		{
			for (size_t i = ii + (ii==0); i < ii+I_BLOCK_SIZE && i < N+1; ++i)
			{
				char a = A[i - 1];
				for (size_t j = jj + (jj==0); j < jj+J_BLOCK_SIZE && j < M+1; ++j)
				{
					int h = H[i * (M + 1) + j];
					int h_ne = H[(i - 1) * ( M + 1) + j - 1];

					h = max(h, h_ne + score_match*(a==B[j-1]) + score_mismatch*(a!=B[j-1]));
					h = max(h, Mj[j]);
					h = max(h, Mi[i - ii]);

					ans = max(ans, h);

					Mj[j] = max(Mj[j] + score_continue_gap, h + score_open_gap);
					Mi[i - ii] = max(Mi[i - ii] + score_continue_gap, h + score_open_gap);
					H[i * (M + 1) + j] = h;
				}
			}
		}
	}
	
	free(Mj);
	return ans;
}

long smith_waterman_flops_quadratic_parallel(int N, int M)
{
	return (long)10*(N)*M;
}
