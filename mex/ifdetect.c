#include <math.h>
#include <stdio.h>

#define FP_TYPE double
#include "../external/ciglet/ciglet.c"

#include <mex.h>

/* [f] = ifdetect(x, center, fc, fres)

  Inputs:
   x            1d vec    real
   center       1d vec    real
   fc           scalar    real
   fres         scalar    real
   
  Outputs:
   f            1d vec    real
*/

#define max(a, b) ((a) > (b) ? (a) : (b))
#define min(a, b) ((a) < (b) ? (a) : (b))
#define mxGetL(pr) max(mxGetN(pr), mxGetM(pr))

static double fftbuff[65536];
void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[]) {
  if(nrhs != 4)
    mexErrMsgTxt("Incorrect number of inputs.");
  if(nlhs != 1)
    mexErrMsgTxt("Incorrect number of outputs.");

  double* x = mxGetPr(prhs[0]);
  mwSize nx = mxGetL(prhs[0]);
  double* center = mxGetPr(prhs[1]);
  mwSize nf = mxGetL(prhs[1]);
  double fc = mxGetScalar(prhs[2]);
  double fres = mxGetScalar(prhs[3]);

  plhs[0] = mxCreateDoubleMatrix(nf, 1, mxREAL);
  double* f = mxGetPr(plhs[0]);
  
  ifdetector* ifd = create_ifdetector(fc, fres);
  for(int i = 0; i < nf; i ++) {
    double* xi = fetch_frame(x, nx, (int)center[i], ifd -> nh);
    f[i] = ifdetector_estimate(ifd, xi, ifd -> nh);
    free(xi);
  }
  delete_ifdetector(ifd);
}

