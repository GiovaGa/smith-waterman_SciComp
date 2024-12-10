
#include<stdlib.h>
#include <stdio.h>
#include<omp.h>

#include "sequence.h"

static inline int min(const int a, const int b)
{
	if (a < b) return a;
	return b;
}

static inline int max(const int a, const int b)
{
	if (a > b) return a;
	return b;
}

static inline void fill(int v, int* p, size_t n)
{
	for(size_t k = 0; k < n; ++k)
		p[k] = v;
}

#define I_BLOCK_SIZE 8
#define J_BLOCK_SIZE (64/sizeof(int))


int sw_quad_par_atomic(const struct sequence_t* A, const struct sequence_t* B, const struct scores_t*scores) {
	//if(B->length < A->length) return smith_waterman_quadratic_opt(B->length, N, B, A, scores->match, scores->gap_opening, scores->gap_extension, scores->mismatch, H);
	int* H = (int*)malloc((A->length + 1) * (B->length + 1) * sizeof(int));
	if (!H)
	{
		fprintf(stderr, "Cannot allocate memory for score matrix\n");
		abort();
	}
	
	int ans = 0;
    int *Mj = malloc((B->length+1)*sizeof(int));
	fill(scores->gap_opening, Mj, B->length+1);

	size_t rows = (A->length+1) / I_BLOCK_SIZE + (((A->length+1) % I_BLOCK_SIZE) != 0);
	size_t cols = (B->length+1) / J_BLOCK_SIZE + (((B->length+1) % J_BLOCK_SIZE) != 0);
	
	size_t** counter_ptrs;

#pragma omp parallel reduction(max:ans)
	{
		int num_threads = omp_get_num_threads();
		int thread_id = omp_get_thread_num();
		int prev_thread_id = (thread_id == 0) ? (num_threads - 1) : (thread_id - 1);
		#pragma omp single
		{
			counter_ptrs = malloc(num_threads * sizeof(size_t*));
		}
		size_t counter = 0;
		counter_ptrs[thread_id] = &counter;
		#pragma omp barrier
		size_t* prev_counter = counter_ptrs[prev_thread_id];
		
		
		for (size_t ii = thread_id; ii < rows; ii += num_threads)
		{
			int Mi[I_BLOCK_SIZE];
			fill(scores->gap_opening, Mi, I_BLOCK_SIZE);
			
			for (size_t jj = 0; jj < cols; ++jj)
			{
				if (thread_id != 0)
				{
					size_t x;
					do
					{
						#pragma omp atomic read
						x = *prev_counter;
					}
					while(x <= counter); 
				}
				else if (ii != 0)
				{
					size_t x;
					do
					{
						#pragma omp atomic read
						x = *prev_counter;
					}
					while(x <= counter - cols); 	
				}
				
				
				for (size_t i = ii*I_BLOCK_SIZE + (ii==0); i < (ii+1)*I_BLOCK_SIZE && i < A->length+1; ++i)
				{
					const char a = A->data[i - 1];
					for (size_t j = jj*J_BLOCK_SIZE + (jj==0); j < (jj+1)*J_BLOCK_SIZE && j < B->length+1; ++j)
					{
						// TODOOOOOOO::: fix
						int h_ne = H[(i - 1) * ( B->length + 1) + j - 1] * (i!=1) * (j!=1);

						int h = max(0, h_ne + scores->match*(a==B->data[j-1]) + scores->mismatch*(a!=B->data[j-1]));
						h = max(h, Mj[j]);
						h = max(h, Mi[i - ii*I_BLOCK_SIZE]);

						ans = max(ans, h);

						Mj[j] = max(Mj[j] + scores->gap_extension, h + scores->gap_opening);
						Mi[i - ii*I_BLOCK_SIZE] = max(Mi[i - ii*I_BLOCK_SIZE] + scores->gap_extension, h + scores->gap_opening);
						H[i * (B->length + 1) + j] = h;
					}
				}
				++counter;
			}
		}
	}
	
	free(Mj);
	free(counter_ptrs);
	free(H);
	return ans;
}

size_t sw_quad_par_atomic_flops(size_t N, size_t M)
{
	return 10*N*M;
}
