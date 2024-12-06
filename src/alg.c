#include<omp.h>
#include"benchmark.h"
#include"sequence.h"

static inline int max(const int a, int b) {
    return (a >= b)? a:b;
}


int smith_waterman(const struct sequence_t* A, const struct sequence_t* B, const struct scores_t*scores, int* restrict H) {
	int ans = 0;
	// int* W; // Gap weight

	for (size_t i = 1; i <= A->length; ++i) {
		for (size_t j = 1; j <= B->length; ++j) {
			const int score = scores->match * (A->data[i-1] == B->data[j-1]) + scores->mismatch * (A->data[i-1] != B->data[j-1]);
			H[i * (A->length + 1) + j] = max(H[i * (A->length + 1) + j], H[(i - 1) * (A->length + 1) + j - 1] + score);
			for (size_t k = 1; k < i; ++k) {
				H[i * (A->length + 1) + j] = max(H[i * (A->length + 1) + j], H[(i - k) * (A->length + 1) + j] + scores->gap_opening + scores->gap_extension * (k-1));
			}
			for (size_t k = 1; k < j; ++k) {
				H[i * (A->length + 1) + j] = max(H[i * (A->length + 1) + j], H[i * (A->length + 1) + j - k] + scores->gap_opening + scores->gap_extension * (k-1));
			}
			ans = max(ans, H[i * (A->length + 1) + j]);
		}
	}
	return ans;
}

long smith_waterman_flops(int N, int M)
{
	return (long)(N)*M;
}
