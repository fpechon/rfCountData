/*******************************************************************
 /*******************************************************************
  Copyright (C) 2001-7 Leo Breiman, Adele Cutler and Merck & Co., Inc.
  
  This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
*******************************************************************/
 
 /******************************************************************
 * buildtree and findbestsplit routines translated from Leo's
 * original Fortran code.
 *
 *      copyright 1999 by leo Breiman
 *      this is free software and can be used for any purpose.
 *      It comes with no guarantee.
 *
 ******************************************************************/
#include <Rmath.h>
#include <R.h>
#include "rf.h"
 
 
 void regTree(double *x, double *offset, double *y, int mdim, int nsample, int *lDaughter,
              int *rDaughter,
              double *upper, double *avnode, int *nodestatus, int nrnodes,
              int *treeSize, int nthsize, int mtry, int *mbest, int *cat,
              double *tgini, int *varUsed) {
   int i, j, k, m, ncur, *jdex, *nodestart, *nodepop;
   int ndstart, ndend, ndendl, jstat, msplit;
   double d, d_offset, ss, av, decsplit, ubest, sumnode, nodecnt;
   
   nodestart = (int *) Calloc(nrnodes, int);
   nodepop   = (int *) Calloc(nrnodes, int);
   
   /* initialize some arrays for the tree */
   
   zeroInt(nodestatus, nrnodes);
   zeroInt(nodestart, nrnodes);
   zeroInt(nodepop, nrnodes);
   zeroDouble(avnode, nrnodes);
   
   jdex = (int *) Calloc(nsample, int);
   for (i = 1; i <= nsample; ++i) jdex[i-1] = i;
   
   ncur = 0;
   nodestart[0] = 0;
   nodepop[0] = nsample;
   nodestatus[0] = NODE_TOSPLIT;
   
   /* compute mean and sum of squares for Y 
   * av = sum x_i / sum offset_i
   * ss = - Poisson LL
   */
   av = 0.0;
   d = 0.0;
   d_offset = 0.0;
   ss = 0.0;
   for (i = 0; i < nsample; ++i) {
     d += y[jdex[i] - 1];
     d_offset += exp(offset[jdex[i]-1]);
     ss += exp(offset[jdex[i]-1])  - y[jdex[i]-1] *offset[jdex[i]-1] ;
   }
   av = d / d_offset ;
   avnode[0] = av;
   
   /* start main loop */
   
   
   for (k = 0; k < nrnodes - 2; ++k) {
     if (k > ncur || ncur >= nrnodes - 2) break;
     /* skip if the node is not to be split */
     if (nodestatus[k] != NODE_TOSPLIT) continue;
     
#ifdef RF_DEBUG
     Rprintf("regTree: k=%d, av=%f, ss=%f\n", k, av, ss);
#endif
     /* initialize for next call to findbestsplit */
     ndstart = nodestart[k];
     ndend = ndstart + nodepop[k] - 1;
     d_offset = 0.0;
     for (j = 0; j < nsample; ++j) d_offset += exp(offset[jdex[j]-1]);
     nodecnt = d_offset;
     sumnode = nodecnt * avnode[k];
     jstat = 0;
     decsplit = 0.0;
     
#ifdef RF_DEBUG
     Rprintf("before findBestSplit: ndstart=%d, ndend=%d, jstat=%d, decsplit=%f\n",
             ndstart, ndend, jstat, decsplit);
#endif
     findBestSplit(x, offset, jdex, y, mdim, nsample, ndstart, ndend, &msplit,
                   &decsplit, &ubest, &ndendl, &jstat, mtry, sumnode,
                   nodecnt, cat); 
#ifdef RF_DEBUG
     Rprintf(" after findBestSplit: ndstart=%d, ndend=%d, jstat=%d, decsplit=%f, msplit=%d\n",
             ndstart, ndend, jstat, decsplit, msplit);
#endif
     if (jstat == 1) {
       /* Node is terminal: Mark it as such and move on to the next. */
       nodestatus[k] = NODE_TERMINAL;
       continue;
     }
     /* Found the best split. */
     mbest[k] = msplit;
     varUsed[msplit - 1] = 1;
     upper[k] = ubest;
     tgini[msplit - 1] += decsplit;
     nodestatus[k] = NODE_INTERIOR;
     
     /* leftnode no.= ncur+1, rightnode no. = ncur+2. */
     nodepop[ncur + 1] = ndendl - ndstart + 1;
     nodepop[ncur + 2] = ndend - ndendl;
     nodestart[ncur + 1] = ndstart;
     nodestart[ncur + 2] = ndendl + 1;
     
     /* compute mean and sum of squares for the left daughter node 
     * av = sum y_i / sum offset_i
     * ss = - Poisson LL
     */
     
     av = 0.0;
     d = 0.0;
     d_offset = 0.0;
     ss = 0.0;
     for (j = ndstart; j <= ndendl; ++j) {
       d += y[jdex[j]-1];
       d_offset += exp(offset[jdex[j]-1]);
       ss += exp(offset[jdex[j]-1])  - y[jdex[j]-1] *offset[jdex[j]-1] ;
     }
     av = d / d_offset;
     avnode[ncur+1] = av;
     
     nodestatus[ncur+1] = NODE_TOSPLIT;
     if (nodepop[ncur + 1] <= nthsize) {
       nodestatus[ncur + 1] = NODE_TERMINAL;
     }
     /* compute mean and sum of squares for the right daughter node
     * av = sum y_i / sum offset_i
     * ss = - Poisson LL
     */
     av = 0.0;
     ss = 0.0;
     for (j = ndendl + 1; j <= ndend; ++j) {
       d += y[jdex[j]-1];
       d_offset += exp(offset[jdex[j]-1]);
       ss += exp(offset[jdex[j]-1])  - y[jdex[j]-1] *offset[jdex[j]-1] ;
     }
     av = d / d_offset;
     avnode[ncur + 2] = av;
     
     nodestatus[ncur + 2] = NODE_TOSPLIT;
     if (nodepop[ncur + 2] <= nthsize) {
       nodestatus[ncur + 2] = NODE_TERMINAL;
     }
     
     /* map the daughter nodes */
     lDaughter[k] = ncur + 1 + 1;
     rDaughter[k] = ncur + 2 + 1;
     /* Augment the tree by two nodes. */
     ncur += 2;
#ifdef RF_DEBUG
     Rprintf(" after split: ldaughter=%d, rdaughter=%d, ncur=%d\n",
             lDaughter[k], rDaughter[k], ncur);
#endif
     
   }
   *treeSize = nrnodes;
   for (k = nrnodes - 1; k >= 0; --k) {
     if (nodestatus[k] == 0) (*treeSize)--;
     if (nodestatus[k] == NODE_TOSPLIT) {
       nodestatus[k] = NODE_TERMINAL;
     }
   }
   
   Free(nodestart);
   Free(jdex);
   Free(nodepop);
 }
 
 /*--------------------------------------------------------------*/
 void findBestSplit(double *x, double *offset, int *jdex, double *y, int mdim, int nsample,
                    int ndstart, int ndend, int *msplit, double *decsplit,
                    double *ubest, int *ndendl, int *jstat, int mtry,
                    double sumnode, double nodecnt, int *cat) {
   int last, ncat[MAX_CAT], icat[MAX_CAT], lc, nl, nr;
   int i, j, kv, l, *mind, *ncase;
   double *xt, *ut, *v, *yl, *offsetl, sumcat[MAX_CAT], sumcat_offset[MAX_CAT], avcat[MAX_CAT], tavcat[MAX_CAT], ubestt;
   double crit, critmax, critvar, suml, sumr, d, critParent, npopl, npopr;
   
   ut = (double *) Calloc(nsample, double);
   xt = (double *) Calloc(nsample, double);
   v  = (double *) Calloc(nsample, double);
   yl = (double *) Calloc(nsample, double);
   offsetl = (double *) Calloc(nsample, double);
   mind  = (int *) Calloc(mdim, int);
   ncase = (int *) Calloc(nsample, int);
   zeroDouble(avcat, MAX_CAT);
   zeroDouble(tavcat, MAX_CAT);
   
   /* START BIG LOOP */
   *msplit = -1;
   *decsplit = 0.0;
   critmax = 0.0;
   ubestt = 0.0;
   for (i=0; i < mdim; ++i) mind[i] = i;
   
   last = mdim - 1;
   for (i = 0; i < mtry; ++i) {
     critvar = 0.0;
     j = (int) (unif_rand() * (last+1));
     kv = mind[j];
     swapInt(mind[j], mind[last]);
     last--;
     
     lc = cat[kv];
     if (lc == 1) {
       /* numeric variable */
       for (j = ndstart; j <= ndend; ++j) {
         xt[j] = x[kv + (jdex[j] - 1) * mdim];
         yl[j] = y[jdex[j] - 1];
         offsetl[j] = exp(offset[jdex[j] - 1]);
       }
     } else {
       /* categorical variable */
       zeroInt(ncat, MAX_CAT);
       zeroDouble(sumcat, MAX_CAT);
       zeroDouble(sumcat_offset, MAX_CAT);
       for (j = ndstart; j <= ndend; ++j) {
         l = (int) x[kv + (jdex[j] - 1) * mdim];
         sumcat[l - 1] += y[jdex[j] - 1];
         sumcat_offset[l-1] += exp(offset[jdex[j]-1]);
         ncat[l - 1] ++;
       }
       /* Compute means of Y by category. */
       for (j = 0; j < lc; ++j) {
         avcat[j] = sumcat_offset[j] ? sumcat[j] / sumcat_offset[j] : 0.0;
       }
       /* Make the category mean the `pseudo' X data. */
       for (j = 0; j < nsample; ++j) {
         xt[j] = avcat[(int) x[kv + (jdex[j] - 1) * mdim] - 1];
         yl[j] = y[jdex[j] - 1];
         offsetl[j] = exp(offset[jdex[j] - 1]);
       }
     }
     /* copy the x data in this node. */
     for (j = ndstart; j <= ndend; ++j) v[j] = xt[j];
     for (j = 1; j <= nsample; ++j) ncase[j - 1] = j;
     R_qsort_I(v, ncase, ndstart + 1, ndend + 1);
     if (v[ndstart] >= v[ndend]) continue;
     /* ncase(n)=case number of v nth from bottom */
     /* Start from the right and search to the left. */
     critParent = sumnode * sumnode / nodecnt;
     suml = 0.0;
     sumr = sumnode;
     npopl = 0.0;
     npopr = nodecnt;
     crit = 0.0;
     /* Search through the "gaps" in the x-variable. */
     for (j = ndstart; j <= ndend - 1; ++j) {
       d = yl[ncase[j] - 1];
       suml += d;
       sumr -= d;
       npopl += offsetl[ncase[j]-1];
       npopr -= offsetl[ncase[j]-1];
       if (v[j] < v[j+1]) {
         /* Poisson deviance in both nodes
         * in left node: suml claims with npopl exposure -> suml / npopl claim frequency
         * sum_i lambda_i - n_i * log(lambda_i);
         * but lambda_i = d_i lambda
         * then we have
         *  lambda_1 sum_i d_i - sum_i n_i * log(d_i) + sum_i n_i * log(lambda_1)
         * = sum_i n_i         - sum_i n_i * log(d_i) + sum_i n_i * log(lambda_1)
         * in right node:
         * lambda_2 sum_i d_i - sum_i n_i * log(d_i) + sum_i n_i * log(lambda_2)
         * 
         * lambda. - suml * log(lambda.) where lambda. is the aggregated claim frequency
         * i.e. lambda. = sum_l / popl * sum_i offsetl_i  = sum_l
         *  log-likelihood in left node is  : log(lambda_i)  sum_i y_i
         *  log-likelihood in right node is : sum_r - sum_r * log(sum_r)
         */
         crit = -sumnode * log(sumnode / nodecnt) 
         +suml * log((suml + (suml ==0)) / npopl) + sumr * log((sumr +(sumr==0)) / npopr);
         if (crit > critvar) {
           ubestt = (v[j] + v[j+1]) / 2.0;
           critvar = crit;
         }
       }
     }
     if (critvar > critmax) {
       *ubest = ubestt;
       *msplit = kv + 1;
       critmax = critvar;
       for (j = ndstart; j <= ndend; ++j) {
         ut[j] = xt[j];
       }
       if (cat[kv] > 1) {
         for (j = 0; j < cat[kv]; ++j) tavcat[j] = avcat[j];
       }
     }
   }
   *decsplit = critmax;
   
   /* If best split can not be found, set to terminal node and return. */
   if (*msplit != -1) {
     nl = ndstart;
     for (j = ndstart; j <= ndend; ++j) {
       if (ut[j] <= *ubest) {
         nl++;
         ncase[nl-1] = jdex[j];
       }
     }
     *ndendl = imax2(nl - 1, ndstart);
     nr = *ndendl + 1;
     for (j = ndstart; j <= ndend; ++j) {
       if (ut[j] > *ubest) {
         if (nr >= nsample) break;
         nr++;
         ncase[nr - 1] = jdex[j];
       }
     }
     if (*ndendl >= ndend) *ndendl = ndend - 1;
     for (j = ndstart; j <= ndend; ++j) jdex[j] = ncase[j];
     
     lc = cat[*msplit - 1];
     if (lc > 1) {
       for (j = 0; j < lc; ++j) {
         icat[j] = (tavcat[j] < *ubest) ? 1 : 0;
       }
       *ubest = pack(lc, icat);
     }
   } else *jstat = 1;
   
   Free(ncase);
   Free(mind);
   Free(v);
   Free(yl);
   Free(offsetl);
   Free(xt);
   Free(ut);
 }
 
 /*====================================================================*/
 void predictRegTree(double *x,double *offset, int nsample, int mdim,
                     int *lDaughter, int *rDaughter, int *nodestatus,
                     double *ypred, double *split, double *nodepred,
                     int *splitVar, int treeSize, int *cat, int maxcat,
                     int *nodex) {
   int i, j, k, m, *cbestsplit;
   double dpack;
   
   /* decode the categorical splits */
   if (maxcat > 1) {
     cbestsplit = (int *) Calloc(maxcat * treeSize, int);
     zeroInt(cbestsplit, maxcat * treeSize);
     for (i = 0; i < treeSize; ++i) {
       if (nodestatus[i] != NODE_TERMINAL && cat[splitVar[i] - 1] > 1) {
         dpack = split[i];
         /* unpack `npack' into bits */
         /* unpack(dpack, maxcat, cbestsplit + i * maxcat); */
         for (j = 0; j < cat[splitVar[i] - 1]; ++j) {
           cbestsplit[j + i*maxcat] = ((unsigned long) dpack & 1) ? 1 : 0;
           dpack = dpack / 2.0 ;
           /* cbestsplit[j + i*maxcat] = npack & 1; */
         }
       }
     }
   }
   
   for (i = 0; i < nsample; ++i) {
     k = 0;
     while (nodestatus[k] != NODE_TERMINAL) { /* go down the tree */
           m = splitVar[k] - 1;
       if (cat[m] == 1) {
         k = (x[m + i*mdim] <= split[k]) ?
         lDaughter[k] - 1 : rDaughter[k] - 1;
       } else {
         /* Split by a categorical predictor */
         k = cbestsplit[(int) x[m + i * mdim] - 1 + k * maxcat] ?
         lDaughter[k] - 1 : rDaughter[k] - 1;
       }
     }
     /* terminal node: assign prediction and move on to next */
     ypred[i] = nodepred[k] * exp(offset[i]) ;
     nodex[i] = k + 1;
   }
   if (maxcat > 1) Free(cbestsplit);
 }
 