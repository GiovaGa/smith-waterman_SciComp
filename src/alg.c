#include<stdlib.h>

inline int min(int a,int b){
    if(a <= b) return a;
    else return b;
}

inline int max(int a,int b){
    if(a >= b) return a;
    else return b;
}

int smith_waterman(int N, int M, char*A, char*B, int score_match, int score_skip, int score_mismatch){

    int ans = 0;
    // int* W; // Gap weight
    int* H = malloc((N+1)*(M+1)*sizeof(int)); // Score matrix
    for(int i = 0;i <= N;++i){
        for(int j = 0;j <= M;++j){
            H[i*M+j] = 0;
        }
    }
    for(int i = 1;i <= N;++i){
        for(int j = 1;j <= M;++j){
            const int score = score_match*(A[i] == B[j]) + score_mismatch*(A[i] != B[j]);
            H[i*(M+1)+j] = max(H[i*(M+1)+j], H[(i-1)*(M+1)+j-1] + score);
            for(int k = 1;k < i;++k){
                H[i*(M+1)+j] = max(H[i*(M+1)+j], H[(i-k)*(M+1)+j] - score_skip*k);
            }
            for(int k = 1;k < j;++k){
                H[i*(M+1)+j] = max(H[i*(M+1)+j], H[i*(M+1)+j-k] - score_skip*k);
            }
            ans = max(ans, H[i*(M+1)+j]);
        }
    }
    return ans;
}
