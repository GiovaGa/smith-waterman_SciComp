
#include<stdlib.h>
#include <stdio.h>
#include<omp.h>

#include "sequence.h"

static inline int max(const int a, const int b)
{
	if (a > b) return a;
	return b;
}

int sw_cub_par(const struct sequence_t* A, const struct sequence_t* B, const struct scores_t*scores, int* restrict H)
{
	// int* H = (int*)malloc((A->length + 1) * (B->length + 1) * sizeof(int));
	if (!H)
	{
		fprintf(stderr, "Cannot allocate memory for score matrix\n");
		abort();
	}
	
	int ans = 0;
	for (size_t j = 1; j <= B->length; ++j)
	{
		const char b = B->data[j-1];
#pragma omp parallel for reduction (max:ans) firstprivate(j)
	    for (size_t i = 1; i <= A->length; ++i)
		{
			const int score = scores->match * (A->data[i-1] == b) + scores->mismatch * (A->data[i-1] != b);
			const int h_ne = H[(i - 1) * ( A->length + 1) + j - 1] * (i!=1) * (j!=1);
			
			int h = max(0, h_ne + score);
			
			for (size_t k = 1; k < i; ++k)
				h = max(h, H[(i - k) * (A->length + 1) + j] + scores->gap_opening + scores->gap_extension * (k-1));
			
			for (size_t k = 1; k < j; ++k)
				h = max(h, H[i * (A->length + 1) + j - k] + scores->gap_opening + scores->gap_extension * (k-1));
			
			ans = max(ans, h);
			H[i * (A->length + 1) + j] = h;
		}
	}
	return ans;
}

size_t sw_cub_par_flops(size_t N, size_t M)
{
	return N*M;
}
