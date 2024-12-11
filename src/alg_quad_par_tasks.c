#include<stdlib.h>
#include <stdio.h>
#include<omp.h>

#include "sequence.h"

static inline int max(const int a, const int b)
{
	if (a > b) return a;
	return b;
}

static const int BLOCK_SIZE = 10;  // 50  

// static inline void process_block(size_t i0, size_t j0,const struct sequence_t*A,const struct sequence_t*B, const struct scores_t*scores,int* H,int* Mj, int*ans){ }


int sw_quad_par_tasks(const struct sequence_t* A, const struct sequence_t* B, const struct scores_t*scores)
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
	for(size_t k = 0; k < (B->length+1); ++k) Mj[k] = scores->gap_opening;

#pragma omp parallel reduction(max:ans)
    {
#pragma omp single
        {
	    for (size_t i0 = 1; i0 <= A->length; i0 += BLOCK_SIZE){
            const size_t j0 = 1;
            const size_t v0 = (i0-1)*(B->length + 1)+ j0-1,
                         v1 = (i0)*(B->length + 1), // j0 + BLOCK_SIZE-1,
                         v2 = (i0+BLOCK_SIZE-1)*(B->length + 1)+ j0-1,
                         v3 = (i0+BLOCK_SIZE)*(B->length + 1); // + j0 + BLOCK_SIZE-1 // vertices of the block
#pragma omp task shared(ans,H,Mj,A,B,scores) firstprivate(i0,j0) default(none) depend(in:H[v0:v1]) depend(out:H[v2:v3]) // depend(in:H[v0:v2])depend(out:H[v2:v3])
        {
            // printf("\n%i: block: %lu,%lu\n", omp_get_thread_num(),i0,j0);
            for (size_t i = i0; i < i0+BLOCK_SIZE; ++i)
            {
                int Mi = scores->gap_opening;
                const char a = A->data[i - 1];
            for (size_t j0 = 1; j0 <= B->length; j0 += BLOCK_SIZE)
                for (size_t j = j0; j < j0+BLOCK_SIZE; ++j)
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
          }
        }
      }
    }
	free(Mj);
	free(H);
	return ans;
}
