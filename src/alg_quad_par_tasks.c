#include<stdlib.h>
#include <stdio.h>
#include<omp.h>

#include "sequence.h"

static inline int max(const int a, const int b)
{
	if (a > b) return a;
	return b;
}

static const int BLOCK_SIZE = 200;

static int recursive_solve(const size_t i0, const size_t i1, const size_t j0, const size_t j1,
                           const struct sequence_t*A,const struct sequence_t*B, const struct scores_t*scores,int* H,int* Mj,int* Mi){
    // printf("%lu -- %lu, %lu -- %lu",i0,i1,j0,j1); 
    int ans = 0,ans1,ans2;
    if((i1 > i0+BLOCK_SIZE || j1 > j0+BLOCK_SIZE)&&(i1 > i0+BLOCK_SIZE && j1 > j0+BLOCK_SIZE)){
        const size_t im = (i0+i1)/2, jm = (j0+j1)/2;
        ans  = recursive_solve(i0, im, j0, jm, A, B, scores, H, Mj, Mi);
    #pragma omp task shared(ans1, im, i1, j0, jm, A, B, scores, H, Mj, Mi)
        ans1 = recursive_solve(im, i1, j0, jm, A, B, scores, H, Mj, Mi);
    #pragma omp task shared(ans2, i0, im, jm, j1, A, B, scores, H, Mj, Mi)
        ans2 = recursive_solve(i0, im, jm, j1, A, B, scores, H, Mj, Mi);
    #pragma omp taskwait
        ans  = max(max(ans,ans1),max(ans2, recursive_solve(im, i1, jm,j1, A, B, scores, H, Mj, Mi)));
    }else{
        // if(i0 >= i1 || j0 >= j1) return 0;
        // actually solve
        for (size_t i = i0; i < i1; ++i)
        {
            const char a = A->data[i - 1];
            for (size_t j = j0; j < j1; ++j)
            {
                const int score = scores->match * (a == B->data[j-1]) + scores->mismatch * (a != B->data[j-1]);
                const int h_ne = H[(i - 1) * ( B->length + 1) + j - 1] * (i!=1) * (j!=1);
                int h = max(0, h_ne + score);
                h = max(h, Mj[j]); h = max(h, Mi[i]);

                ans = max(ans, h);

                Mj[j] = max(Mj[j] + scores->gap_extension, h + scores->gap_opening);
                Mi[i] = max(Mi[i] + scores->gap_extension, h + scores->gap_opening);
                H[i * (B->length + 1) + j] = h;
            }
        }
    }
    return ans;
}


int sw_quad_par_tasks(const struct sequence_t* A, const struct sequence_t* B, const struct scores_t*scores, int* restrict H)
{
	//if(B->length < N) return smith_waterman_quadratic_opt(B->length, N, B, A, scores->match, scores->gap_opening, scores->gap_extension, scores->mismatch, H);
	
	int ans = 0;
    int *Mi = malloc((A->length+1)*sizeof(int));
	for(size_t k = 0; k < (A->length+1); ++k) Mi[k] = scores->gap_opening;
    int *Mj = malloc((B->length+1)*sizeof(int));
	for(size_t k = 0; k < (B->length+1); ++k) Mj[k] = scores->gap_opening;

#pragma omp parallel
    {
#pragma omp single
        {
        ans = recursive_solve(1, A->length+1, 1, B->length+1, A, B, scores, H, Mj, Mi);
        }
    }
	free(Mi);
	free(Mj);
	return ans;
}
