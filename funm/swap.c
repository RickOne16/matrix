/*
 *  C mex file for MATLAB that implements LAPACK ctrexc_ for
 *  reordering the complex Schur decomposition A = UTU^* such that
 *  T is in block form. This program combines Cswap and Crealswap
 *  which deal with the complex and the real case respectively.
 *
 *  Input matrices 'U' and 'T' are those from Schur decompositon
 *
 *  Called by [U,T] = swap(U,T,M)
 *
 *  The matrix M is produced using the m-files
 *  >> m = blocking(T,noplot,delta); where delta is some tolerance
 *                                   default: delta = 0.1;
 *  >> [n,M,ind,totalswaps] = swapping5(m);
 *
 *  Output matrices 'U' and 'T' are the updated Schur decomposition
 *
 */

#include "mex.h"
#include "matrix.h"
#include "fort.h"      /* defines mat2fort and fort2mat */

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{

  /* compq=V then the matrix of Schur vectors is updated */

  char *compq = "V", msg[80];
  int n, ldt, ldq, i, info, ifst, ilst, lengthm;
  double *t, *q, *m, *work;
  mxArray  *mm;

  /* expect 3 inputs and 2 outputs */
  if ((nrhs != 3) || (nlhs != 2)){
     mexErrMsgTxt("Expected 3 inputs and 2 outputs");
  }

  /* Sizes of stuff */

  n = mxGetN( prhs[0] );
  ldt = n;
  ldq = n;
  lengthm = mxGetM( prhs[2] );

  /* Copy input and output matrices to stuff */

  mm = mxDuplicateArray( prhs[2] );
  m = mxGetPr(mm);

  if (mxIsComplex( prhs[1] )) {

      /* Convert complex input to complex FORTRAN format */

      q = mat2fort( prhs[0], ldq, n );
      t = mat2fort( prhs[1], ldt, n ); }

  else {

      plhs[0] = mxCreateDoubleMatrix(n,n,mxREAL);
      plhs[1] = mxCreateDoubleMatrix(n,n,mxREAL);

      plhs[1] = mxDuplicateArray( prhs[1] );
      t = mxGetPr(plhs[1]);
      plhs[0] = mxDuplicateArray( prhs[0] );
      q = mxGetPr(plhs[0]);
  }

  /* Allocate workspace */

  work = (double *)mxCalloc(n,sizeof(double));

  /* Do the loop */

  for ( i = 0; i < lengthm; i++) {
    info = 0;
    ifst = m[lengthm + i];
    ilst = m[i];
    if (mxIsComplex( prhs[1] )) {
        ztrexc(compq,&n,t,&ldt,q,&ldq,&ifst,&ilst,&info); }
    else {
        dtrexc(compq,&n,t,&ldt,q,&ldq,&ifst,&ilst,work,&info);
    }
    if (info < 0){
      sprintf(msg, "The routine DTREXC has detected an error");
      mexErrMsgTxt(msg);
    }
  }

  /* Convert output to MATLAB format */

  if (mxIsComplex( prhs[1] )) {
      plhs[0] = fort2mat( q, ldq, ldq, n );
      plhs[1] = fort2mat( t, ldt, ldt, n );
  }

  /* Free up memory */

  mxFree(work);
}

