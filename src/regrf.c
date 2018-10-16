/*******************************************************************
 Copyright (C) 2001-2015 Leo Breiman, Adele Cutler and Merck & Co., Inc.
 
 This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
*******************************************************************/

#include <R.h>
#include "rf.h"

double poissondev(double y, double y_pred){
  return y_pred - y + y * log(y + (y==0)) - y * log(y_pred+(y_pred==0));
}


void regRF(double *x,double *offset, double *y, int *xdim, int *sampsize,
           int *nthsize, int *nrnodes, int *nTree, int *mtry, int *imp,
           int *cat, int *maxcat, int *jprint, 
           double *yptr, double *errimp, 
           int *treeSize, int *nodestatus,
           int *lDaughter, int *rDaughter, double *avnode, int *mbest,
           double *upper, double *dev, int *keepf, int *replace,
           int *testdat, double *xts, int *nts, double *yts, double *offsetts, int *labelts,
           double *yTestPred, 
           double *devts, double *coef,
           int *nout, int *inbag) {
  /*************************************************************************
  Input:
  mdim=number of variables in data set
  nsample=number of cases
  
  nthsize=number of cases in a node below which the tree will not split,
  setting nthsize=5 generally gives good results.
  
  nTree=number of trees in run.  200-500 gives pretty good results
  
  mtry=number of variables to pick to split on at each node.  mdim/3
  seems to give genrally good performance, but it can be
  altered up or down
  
  imp=1 turns on variable importance.  This is computed for the
  mth variable as the percent rise in the test set mean sum-of-
  squared errors when the mth variable is randomly permuted.
  
  *************************************************************************/
  double errts = 0.0, averrb, meanY, meanYts, varY, varYts, r, xrand,
    errb = 0.0, resid=0.0, ooberr, ooberrperm, delta, *resOOB;
  
  double *yb,*offsetb, *xtmp, *xb, *ytr, *ytree, *tgini;
  
  int k, m, mr, n, nOOB, j, jout, idx, ntest, last, ktmp, nPerm,
  nsample, mdim, keepF, keepInbag;
  int *oobpair, varImp, *varUsed; 
  
  int *in, *nind, *nodex, *nodexts;
  
  nsample = xdim[0];
  mdim = xdim[1];
  ntest = *nts;
  varImp = imp[0];
  nPerm = imp[1];
  keepF = keepf[0];
  keepInbag = keepf[1];
  
  if (*jprint == 0) *jprint = *nTree + 1;
  
  yb         = (double *) S_alloc(*sampsize, sizeof(double));
  offsetb    = (double *) S_alloc(*sampsize, sizeof(double));
  xb         = (double *) S_alloc(mdim * *sampsize, sizeof(double));
  ytr        = (double *) S_alloc(nsample, sizeof(double));
  xtmp       = (double *) S_alloc(nsample, sizeof(double));
  resOOB     = (double *) S_alloc(nsample, sizeof(double));
  
  in        = (int *) S_alloc(nsample, sizeof(int));
  nodex      = (int *) S_alloc(nsample, sizeof(int));
  varUsed    = (int *) S_alloc(mdim, sizeof(int));
  nind = *replace ? NULL : (int *) S_alloc(nsample, sizeof(int));
  
  if (*testdat) {
    ytree      = (double *) S_alloc(ntest, sizeof(double));
    nodexts    = (int *) S_alloc(ntest, sizeof(int));
  }
  
  /* If variable importance is requested, tgini points to the second
  "column" of errimp, otherwise it's just the same as errimp. */
  tgini = varImp ? errimp + mdim : errimp;
  
  averrb = 0.0;
  meanY = 0.0;
  varY = 0.0;
  
  zeroDouble(yptr, nsample);
  zeroInt(nout, nsample);
  for (n = 0; n < nsample; ++n) {
    varY += n * (y[n] - meanY)*(y[n] - meanY) / (n + 1);
    meanY = (n * meanY + y[n]) / (n + 1);
  }
  varY /= nsample;
  
  varYts = 0.0;
  meanYts = 0.0;
  if (*testdat) {
    for (n = 0; n < ntest; ++n) {
      varYts += n * (yts[n] - meanYts)*(yts[n] - meanYts) / (n + 1);
      meanYts = (n * meanYts + yts[n]) / (n + 1);
    }
    varYts /= ntest;
  }

  if (varImp) {
    zeroDouble(errimp, mdim * 2);
  } else {
    zeroDouble(errimp, mdim);
  }
  if (*labelts) zeroDouble(yTestPred, ntest);
  
  /* print header for running output */
  if (*jprint <= *nTree) {
             Rprintf("     | OOB (1 tree)");
             Rprintf(" | OOB (forest) ");
    if (*testdat) Rprintf("|   Test set   ");
    Rprintf("|\n");
    Rprintf("Tree |Loss Function ");
             Rprintf("|Loss Function ");
    if (*testdat) Rprintf("|Loss Function ");
    Rprintf("|\n");
  }
  GetRNGstate();
  /*************************************
  * Start the loop over trees.
  *************************************/
  for (j = 0; j < *nTree; ++j) {
    idx = keepF ? j * *nrnodes : 0;
    zeroInt(in, nsample);
    zeroInt(varUsed, mdim);
    /* Draw a random sample for growing a tree. */
    if (*replace) { /* sampling with replacement */
    for (n = 0; n < *sampsize; ++n) {
      xrand = unif_rand();
      k = xrand * nsample;
      in[k] += 1;
      yb[n] = y[k];
      offsetb[n] = offset[k];
      for(m = 0; m < mdim; ++m) {
        xb[m + n * mdim] = x[m + k * mdim];
      }
    }
    } else { /* sampling w/o replacement */
    for (n = 0; n < nsample; ++n) nind[n] = n;
      last = nsample - 1;
      for (n = 0; n < *sampsize; ++n) {
        ktmp = (int) (unif_rand() * (last+1));
        k = nind[ktmp];
        swapInt(nind[ktmp], nind[last]);
        last--;
        in[k] += 1;
        yb[n] = y[k];
        offsetb[n] = offset[k];
        for(m = 0; m < mdim; ++m) {
          xb[m + n * mdim] = x[m + k * mdim];
        }
      }
    }
    if (keepInbag) {
      for (n = 0; n < nsample; ++n) inbag[n + j * nsample] = in[n];
    }
    /* grow the regression tree */
    regTree(xb, offsetb, yb, mdim, *sampsize, lDaughter + idx, rDaughter + idx,
            upper + idx, avnode + idx, nodestatus + idx, *nrnodes,
            treeSize + j, *nthsize, *mtry, mbest + idx, cat, tgini,
            varUsed);
    /* predict the OOB data with the current tree */
    /* ytr is the prediction on OOB data by the current tree */
    predictRegTree(x, offset, nsample, mdim, lDaughter + idx,
                   rDaughter + idx, nodestatus + idx, ytr, upper + idx,
                   avnode + idx, mbest + idx, treeSize[j], cat, *maxcat,
                   nodex);
    /* yptr is the aggregated prediction by all trees grown so far */
    errb = 0.0;
    ooberr = 0.0;
    jout = 0; /* jout is the number of cases that has been OOB so far */
    nOOB = 0; /* nOOB is the number of OOB samples for this tree */
    for (n = 0; n < nsample; ++n) {
      if (in[n] == 0) {
        nout[n]++;
        nOOB++;
        yptr[n] = ((nout[n]-1) * yptr[n] + ytr[n]) / nout[n];
        resOOB[n] = poissondev(y[n], ytr[n]); 
        ooberr += 2*resOOB[n];
      }
      if (nout[n]) {
        jout++;
        errb +=   poissondev(y[n], yptr[n]);
      }
    }
    ooberr /= nOOB;
    errb /= jout;
    errb *= 2;
    
    /* predict testset data with the current tree */
    if (*testdat) {
      predictRegTree(xts,offsetts, ntest, mdim, lDaughter + idx,
                     rDaughter + idx, nodestatus + idx, ytree,
                     upper + idx, avnode + idx,
                     mbest + idx, treeSize[j], cat, *maxcat, nodexts);
      /* ytree is the prediction for test data by the current tree */
      /* yTestPred is the average prediction by all trees grown so far */
      errts = 0.0;
      for (n = 0; n < ntest; ++n) {
        yTestPred[n] = (j * yTestPred[n] + ytree[n]) / (j + 1);
      }
      /* compute testset Deviance */
      if (*labelts) {
        for (n = 0; n < ntest; ++n) {
          resid = poissondev(yts[n], yTestPred[n]); 
          errts += 2*resid;
        }
        errts /= ntest;
      }
    }
    /* Print running output. */
    if ((j + 1) % *jprint == 0) {
      Rprintf("%4d |", j + 1);
      Rprintf(" %12.6g ", ooberr);
      Rprintf("| %12.6g ", errb);
      if(*labelts == 1) Rprintf("| %12.6g ",
         errts);
      Rprintf("|\n");
    }
    dev[j] = errb;
    if (*labelts) devts[j] = errts;
    /* Variable importance */
    if (varImp) {
      for (mr = 0; mr < mdim; ++mr) {
        if (varUsed[mr]) { /* Go ahead if the variable is used */
    /* make a copy of the m-th variable into xtmp */
    for (n = 0; n < nsample; ++n)
      xtmp[n] = x[mr + n * mdim];
          ooberrperm = 0.0;
          for (k = 0; k < nPerm; ++k) {
            permuteOOB(mr, x, in, nsample, mdim);
            predictRegTree(x, offset, nsample, mdim, lDaughter + idx,
                           rDaughter + idx, nodestatus + idx, ytr,
                           upper + idx, avnode + idx, mbest + idx,
                           treeSize[j], cat, *maxcat, nodex);
            for (n = 0; n < nsample; ++n) {
              if (in[n] == 0) {
                r = poissondev(y[n], ytr[n]);
                ooberrperm += 2*r;
              }
            }
          }
          ooberrperm /=nPerm*nOOB;
          delta = ooberrperm - ooberr;
          errimp[mr] += delta;
          /* copy original data back */
          for (n = 0; n < nsample; ++n)
            x[mr + n * mdim] = xtmp[n];
        }
      }
    }
  }
  PutRNGstate();
  /* end of tree iterations=======================================*/

  if (varImp) {
    for (m = 0; m < mdim; ++m) {
      errimp[m] = errimp[m] / *nTree;
    }
  }
  for (m = 0; m < mdim; ++m) tgini[m] /= *nTree;
}