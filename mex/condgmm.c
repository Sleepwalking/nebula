#include <mex.h>
#include <math.h>
#include <stdio.h>

#define max(a, b) ((a) > (b) ? (a) : (b))
#define min(a, b) ((a) < (b) ? (a) : (b))
#define mxGetL(pr) max(mxGetN(pr), mxGetM(pr))

/*
  L = condgmm(G, f, x, stdfloor)
*/

inline static double lse(double a, double b) {
  double maxab = a > b ? a : b;
  double minab = a > b ? b : a;
  if(maxab - minab > 15)
    return maxab;
  return log(1.0 + exp(minab - maxab)) + maxab;
}

void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[]) {
  if(nrhs != 3 && nrhs != 4)
    mexErrMsgTxt("Incorrect number of inputs.");
  if(nlhs != 1)
    mexErrMsgTxt("Incorrect number of outputs.");
  
  const mxArray* G = prhs[0];
  double* f = mxGetPr(prhs[1]);
  int nf = mxGetL(prhs[1]);
  double* x = mxGetPr(prhs[2]);
  int nx = mxGetL(prhs[2]);
  double stdfloor = 0;
  if(nrhs == 4)
    stdfloor = mxGetScalar(prhs[3]);
  
  int nmix = mxGetL(G);
  double* f0_mu    = mxCalloc(nmix, sizeof(double));
  double* f0_sigma = mxCalloc(nmix, sizeof(double));
  double* f0_lL    = mxCalloc(nmix, sizeof(double));
  double* f0_w     = mxCalloc(nmix, sizeof(double));
  double* dx       = mxCalloc(nx, sizeof(double));
  double lwsum = -1e10;
  
  // convert to conditional GMM
  for(int i = 0; i < nmix; i ++) {
    mxArray* imix = mxGetCell(G, i);
    mxArray* joint_mu_    = mxGetField(imix, 0, "mu");
    mxArray* joint_w_     = mxGetField(imix, 0, "w");
    mxArray* inv11_       = mxGetField(imix, 0, "inv11");
    mxArray* det11_       = mxGetField(imix, 0, "det11");
    mxArray* S21_S11inv_  = mxGetField(imix, 0, "s21s11i");
    mxArray* cond_sigma_  = mxGetField(imix, 0, "condsigma");
    double joint_w = mxGetScalar(joint_w_);
    double* joint_mu = mxGetPr(joint_mu_);
    double cond_sigma = mxGetScalar(cond_sigma_);
    double det11 = fabs(mxGetScalar(det11_));
    double* inv11 = mxGetPr(inv11_);
    double* S21_S11inv = mxGetPr(S21_S11inv_);

    f0_mu[i] = joint_mu[nx];
    f0_sigma[i] = max(stdfloor * stdfloor, cond_sigma);
    for(int k = 0; k < nx; k ++) {
      dx[k] = x[k] - joint_mu[k];
      f0_mu[i] += S21_S11inv[k] * dx[k];
    }
    double dxS11invdx = 0;
    for(int j = 0; j < nx; j ++) {
      double dx_S11inv = 0;
      for(int k = 0; k < nx; k ++)
        dx_S11inv += dx[k] * inv11[k + j * nx];
      dxS11invdx += dx_S11inv * dx[j];
    }

    f0_lL[i] = - 0.5 * dxS11invdx - (nx - 1) * log(2.0 * M_PI) -
      0.5 * log(det11) + log(joint_w);
    lwsum = lse(lwsum, f0_lL[i]);
  }
  
  // reduce dynamic range of mixture weights and normalize
  double wsum = 0;
  for(int i = 0; i < nmix; i ++) {
    double lw = f0_lL[i] - lwsum;
    if(lw < -10) lw = -10;
    f0_w[i] = exp(lw);
    wsum += f0_w[i];
  }
  for(int i = 0; i < nmix; i ++) {
    f0_w[i] /= wsum;
  }
  
  // evaluate conditional likelihood
  plhs[0] = mxCreateDoubleMatrix(nf, 1, mxREAL);
  double* L = mxGetPr(plhs[0]);
  for(int i = 0; i < nf; i ++) {
    double p = 0;
    for(int j = 0; j < nmix; j ++) {
      double d = (f[i] - f0_mu[j]) * (f[i] - f0_mu[j]) / f0_sigma[j];
      p += f0_w[j] * exp(- 0.5 * d) / sqrt(2.0 * M_PI * f0_sigma[j]);
    }
    L[i] = log(p);
  }

  mxFree(f0_mu);
  mxFree(f0_sigma);
  mxFree(f0_lL);
  mxFree(f0_w);
  mxFree(dx);
}

