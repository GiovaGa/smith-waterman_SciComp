#include<stdlib.h>
#include<string.h>
#include<omp.h>
#include"benchmark.h"

static inline int max(const int a, const int b) {
	if (a >= b) return a;
	else return b;
}

int smith_waterman_quadratic(const struct sequence_t* A, const struct sequence_t* B, const struct scores_t*scores, int* restrict H) {
	int ans = 0;
    // int* W; // Gap weight
	int *Mi = malloc((A->length+1)*sizeof(int));
	for(size_t k = 0; k < (A->length+1); ++k)
		Mi[k] = scores->gap_opening;
    int *Mj = malloc((B->length+1)*sizeof(int));
	for(size_t k = 0; k < (B->length+1); ++k)
		Mj[k] = scores->gap_opening;

	for (size_t i = 1; i <= A->length; ++i) {
		for (size_t j = 1; j <= B->length; ++j) {
			const int score = scores->match * (A->data[i-1] == B->data[j-1]) + scores->mismatch * (A->data[i-1] != B->data[j-1]);
			H[i * (B->length + 1) + j] = max(H[i * (B->length + 1) + j], H[(i - 1) * (B->length + 1) + j - 1] + score);

			H[i * (B->length + 1) + j] = max(H[i * (B->length + 1) + j], Mj[j]);
			H[i * (B->length + 1) + j] = max(H[i * (B->length + 1) + j], Mi[i]);

			ans = max(ans, H[i * (B->length + 1) + j]);

            Mj[j] = max(Mj[j] + scores->gap_extension, H[i * (B->length + 1) + j] + scores->gap_opening);
            Mi[i] = max(Mi[i] + scores->gap_extension, H[i * (B->length + 1) + j] + scores->gap_opening);
		}
	}
	
	free(Mi);
	free(Mj);
	return ans;
}

long smith_waterman_flops_quadratic(int N, int M)
{
	return (long)10*(N)*M;
}
