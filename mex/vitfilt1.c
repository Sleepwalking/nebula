#include <mex.h>
#include <math.h>
#include <stdio.h>

#define FP_TYPE float

#include "../external/libgvps/gvps.h"

/* [s L] = vitfilt1(smap, pout, ptran)

  Inputs:
   smap:        2d mat    real
   pout:        2d mat    real
   ptran:       1d vec    real
  
  Outputs:
   s:           1d vec    real
*/

#define max(a, b) ((a) > (b) ? (a) : (b))
#define min(a, b) ((a) < (b) ? (a) : (b))
#define mxGetL(pr) max(mxGetN(pr), mxGetM(pr))

int ntrans = 0;
int fntrans(void* task, int t) {
  return ntrans;
}

double* ptrans = NULL;
FP_TYPE ptransition(void* task, int ds, int t) {
  return ptrans[min(ntrans - 1, ds)];
}

void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[]) {
  if(nrhs != 3)
    mexErrMsgTxt("Incorrect number of inputs.");
  if(nlhs > 2 || nlhs < 1)
    mexErrMsgTxt("Incorrect number of outputs.");

  double* smap = mxGetPr(prhs[0]);
  double* pout = mxGetPr(prhs[1]);
  mwSize nfrm = mxGetM(prhs[1]);
  mwSize nstate = mxGetN(prhs[1]);
  ptrans = mxGetPr(prhs[2]);
  ntrans = mxGetL(prhs[2]);
  mwSize nvirtstate = 0;

  plhs[0] = mxCreateDoubleMatrix(nfrm, 1, mxREAL);
  double* s = mxGetPr(plhs[0]);
  
  gvps_obsrv* obsrv = gvps_obsrv_create(nfrm);
  for(int i = 0; i < nfrm; i ++) {
    obsrv -> slice[i] = gvps_obsrv_slice_create(nstate);
    for(int j = 0; j < nstate; j ++) {
      obsrv -> slice[i] -> pair[j].state = smap[i + j * nfrm];
      obsrv -> slice[i] -> pair[j].p = pout[i + j * nfrm];
      nvirtstate = max(nvirtstate, smap[i + j * nfrm]);
    }
  }
  
  int* path = mxCalloc(nfrm, sizeof(int));
  double L = gvps_sparse_sampled_static(path, NULL, nvirtstate, obsrv,
    ptransition, fntrans);
  for(int i = 0; i < nfrm; i ++)
    s[i] = path[i];
  
  if(nlhs == 2) {
    // output likelihood
    plhs[1] = mxCreateDoubleScalar(L);
  }
  mxFree(path);
  gvps_obsrv_free(obsrv);
}

