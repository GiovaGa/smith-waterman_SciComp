int solve(int i0,  int i1,  int j0,  int j1, char*A, char*B, scores_t*scores,int* H, int*Mj, int*Mi){
   int ans0,ans1,ans2,ans3;
   if(/*block is big enough*/){
       int im = (i0+i1)/2, jm = (j0+j1)/2;
      ans0  = solve(i0, im, j0, jm, A, B, scores, H, Mj, Mi);
#pragma omp task
      ans1 = solve(im, i1, j0, jm, A, B, scores, H, Mj, Mi);
#pragma omp task
      ans2 = solve(i0, im, jm, j1, A, B, scores, H, Mj, Mi);
#pragma omp taskwait
      ans3 = solve(im, i1, jm,j1, A, B, scores, H, Mj, Mi);
      return max(ans0,ans1,ans2,ans3);
   }else ; // solve normally...