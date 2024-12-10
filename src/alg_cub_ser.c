
#include<stdlib.h>
#include <stdio.h>
#include<omp.h>

#include "sequence.h"

static inline int max(const int a, const int b)
{
	if (a > b) return a;
	return b;
}

int sw_cub_ser(const struct sequence_t* A, const struct sequence_t* B, const struct scores_t* scores)
{
	int* H = (int*)malloc((A->length + 1) * (B->length + 1) * sizeof(int));
	if (!H)
	{
		fprintf(stderr, "Cannot allocate memory for score matrix\n");
		abort();
	}
	
	int ans = 0;
	for (size_t i = 1; i <= A->length; ++i)
	{
		const char a = A->data[i-1];
		H[i * (A->length + 1)] = 0;
		for (size_t j = 1; j <= B->length; ++j)
		{
			const int score = scores->match * (a == B->data[j-1]) + scores->mismatch * (a != B->data[j-1]);
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
	free(H);
	return ans;
}

size_t sw_cub_ser_flops(size_t N, size_t M)
{
	return N*M;
}
