#include<stdlib.h>
#include<string.h>
#include<omp.h>
#include <stdbool.h>
#include"benchmark.h"

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


int smith_waterman_quadratic_parallel(const struct sequence_t* A, const struct sequence_t* B, const struct scores_t*scores, int* restrict H) {
	//if(B->length < A->length) return smith_waterman_quadratic_opt(B->length, N, B, A, scores->match, scores->gap_opening, scores->gap_extension, scores->mismatch, H);
	int ans = 0;
    int *Mj = malloc((B->length+1)*sizeof(int));
	fill(scores->gap_opening, Mj, B->length+1);
	
	for (size_t ii = 0; ii <= A->length; ii += I_BLOCK_SIZE)
	{
		int Mi[I_BLOCK_SIZE];
		fill(scores->gap_opening, Mi, I_BLOCK_SIZE);
		
		for (size_t jj = 0; jj <= B->length; jj += J_BLOCK_SIZE)
		{
			for (size_t i = ii + (ii==0); i < ii+I_BLOCK_SIZE && i < A->length+1; ++i)
			{
				const char a = A->data[i - 1];
				for (size_t j = jj + (jj==0); j < jj+J_BLOCK_SIZE && j < B->length+1; ++j)
				{
					int h = H[i * (B->length + 1) + j];
					int h_ne = H[(i - 1) * ( B->length + 1) + j - 1];

					h = max(h, h_ne + scores->match*(a==B->data[j-1]) + scores->mismatch*(a!=B->data[j-1]));
					h = max(h, Mj[j]);
					h = max(h, Mi[i - ii]);

					ans = max(ans, h);

					Mj[j] = max(Mj[j] + scores->gap_extension, h + scores->gap_opening);
					Mi[i - ii] = max(Mi[i - ii] + scores->gap_extension, h + scores->gap_opening);
					H[i * (B->length + 1) + j] = h;
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
