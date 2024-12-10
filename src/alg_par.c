#include<omp.h>
#include"benchmark.h"

static inline int max(const int a, const int b) {
	if (a >= b) return a;
	else return b;
}

int smith_waterman_parallel(const struct sequence_t* A, const struct sequence_t* B, const struct scores_t*scores, int* restrict H) {

	int ans = 0;
	// int* W; // Gap weight

	for (size_t j = 1; j <= B->length; ++j) {
#pragma omp parallel for reduction (max:ans) firstprivate(j)
	    for (size_t i = 1; i <= A->length; ++i) {
			const int score = scores->match * (A->data[i-1] == B->data[j-1]) + scores->mismatch * (A->data[i-1] != B->data[j-1]);
			H[i * (B->length + 1) + j] = max(H[i * (B->length + 1) + j], H[(i - 1) * (B->length + 1) + j - 1] + score);
			for (size_t k = 1; k < i; ++k) {
				H[i * (B->length + 1) + j] = max(H[i * (B->length + 1) + j], H[(i - k) * (B->length + 1) + j] + scores->gap_opening + scores->gap_extension * (k-1));
			}
			for (size_t k = 1; k < j; ++k) {
				H[i * (B->length + 1) + j] = max(H[i * (B->length + 1) + j], H[i * (B->length + 1) + j - k] + scores->gap_opening + scores->gap_extension * (k-1));
			}
			ans = max(ans, H[i * (B->length + 1) + j]);
		}
	}
	return ans;
}

