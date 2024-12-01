
static inline int min_between(int a, int b) {
	if (a <= b) return a;
	else return b;
}

static inline int max_between(int a, int b) {
	if (a >= b) return a;
	else return b;
}

int smith_waterman(int N, int M, const char* restrict A, const char* restrict B, int score_match, int score_skip, int score_mismatch, int* restrict H) {

	int ans = 0;
	// int* W; // Gap weight

	for (int i = 1; i <= N; ++i) {
		for (int j = 1; j <= M; ++j) {
			const int score = score_match * (A[i] == B[j]) + score_mismatch * (A[i] != B[j]);
			H[i * (M + 1) + j] = max_between(H[i * (M + 1) + j], H[(i - 1) * (M + 1) + j - 1] + score);
			for (int k = 1; k < i; ++k) {
				H[i * (M + 1) + j] = max_between(H[i * (M + 1) + j], H[(i - k) * (M + 1) + j] - score_skip * k);
			}
			for (int k = 1; k < j; ++k) {
				H[i * (M + 1) + j] = max_between(H[i * (M + 1) + j], H[i * (M + 1) + j - k] - score_skip * k);
			}
			ans = max_between(ans, H[i * (M + 1) + j]);
		}
	}
	return ans;
}

long smith_waterman_flops(int N, int M)
{
	return (long)(N)*M;
}