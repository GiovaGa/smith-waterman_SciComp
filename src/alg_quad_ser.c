
#include<stdlib.h>
#include <stdio.h>
#include<omp.h>

#include "sequence.h"

static inline int max(const int a, const int b)
{
	if (a > b) return a;
	return b;
}

int sw_quad_ser(const struct sequence_t* A, const struct sequence_t* B, const struct scores_t*scores)
{
	//if(B->length < N) return smith_waterman_quadratic_opt(B->length, N, B, A, scores->match, scores->gap_opening, scores->gap_extension, scores->mismatch, H);
	
	int* H = (int*)malloc((A->length + 1) * (B->length + 1) * sizeof(int));
	if (!H)
	{
		fprintf(stderr, "Cannot allocate memory for score matrix\n");
		abort();
	}
	int ans = 0;
    int *Mj = malloc((B->length+1)*sizeof(int));
	for(size_t k = 0; k < (B->length+1); ++k)
		Mj[k] = scores->gap_opening;

	for (size_t i = 1; i <= A->length; ++i)
	{
		int Mi = scores->gap_opening;
		const char a = A->data[i - 1];
		for (size_t j = 1; j <= B->length; ++j)
		{
			const int score = scores->match * (a == B->data[j-1]) + scores->mismatch * (a != B->data[j-1]);
			
			const int h_ne = H[(i - 1) * ( B->length + 1) + j - 1] * (i!=1) * (j!=1);
			int h = max(0, h_ne + score);
			h = max(h, Mj[j]);
			h = max(h, Mi);

			ans = max(ans, h);

            Mj[j] = max(Mj[j] + scores->gap_extension, h + scores->gap_opening);
            Mi = max(Mi + scores->gap_extension, h + scores->gap_opening);
			H[i * (B->length + 1) + j] = h;
		}
	}
	free(Mj);
	free(H);
	return ans;
}

size_t sw_quad_ser_flops(size_t N, size_t M)
{
	return 10*N*M;
}
