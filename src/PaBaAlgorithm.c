/******************************************************************************
  *
  * calcPaBaMat.c
  *
  * Calculate matrix of all slope pairs for Passing-Bablok algorithm.
  * Slopes are returned as radian angles.
  *
  * Copyright (C) 2011 Roche Diagnostics GmbH
  *
  * This program is free software: you can redistribute it and/or modify
  * it under the terms of the GNU General Public License as published by
  * the Free Software Foundation, either version 3 of the License, or
  * any later version.
  *
  * This program is distributed in the hope that it will be useful,
  * but WITHOUT ANY WARRANTY; without even the implied warranty of
  * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  * GNU General Public License for more details.
  *
  * You should have received a copy of the GNU General Public License
  * along with this program.  If not, see <http://www.gnu.org/licenses/>.
  *
  *****************************************************************************/
#include <R.h>
#include <Rdefines.h>
#include <math.h>
#include <assert.h>
#include <stdlib.h>
typedef enum {false, true} bool;        /* false=0, true=1 */
const double INF = 1020304.050607;      /* evaluated on R-side as Inf */
const double EPS = 1.e-12;
const double PI2 = PI/2;  
const double PI4 = PI/4;
/* function prototypes to avoid implict declaration of functions */
void FillBins(  int * pnSlots, int * posCor, double * pXVals, double * pYVals,
                int * nData, int * nPos, int * nPos2, int * nNeg, int * nNeg2, 
                int * nRegular, int * nAllItems, int * nSlots);
int IndexOf(int nItems, int* pnSlots, int nSlots);
double Tan(double x);
int IndexOf(int nItems, int* pnSlots, int nSlots);
double getMedianSlope(int index1, int index2, int NBins);
double calcDiff(double x,double y) 
{
	double dRes = x-y;
	if(dRes == 0.) return 0.;
	if(fabs(dRes) < EPS*(fabs(x)+fabs(y))/2.) return 0.;
	return dRes;
}
SEXP calcAngleMat(SEXP X, SEXP Y, SEXP posCor) 
{
	// Check length
	assert(length(X)>1);
	assert(length(X)==length(Y));
	// Get X and Y data
	PROTECT(X=AS_NUMERIC(X));
	PROTECT(Y=AS_NUMERIC(Y));
	double* pX = NUMERIC_POINTER(X);
	double* pY = NUMERIC_POINTER(Y);
	unsigned int nData = length(X);
	// Get correlation
	int pCor = INTEGER_VALUE(posCor);
	// Create result matrix
	SEXP ans;
	PROTECT(ans=allocMatrix(REALSXP,nData,nData));
	double* pans = REAL(ans);
	// Calculate results
	double dx,dy;
	int j, k;
	for(j = 0; j < nData; j++) {
		for(k = 0; k < nData; k++)
		{
			if(k > j) {
				dx = calcDiff(pX[k],pX[j]);
				dy = calcDiff(pY[k],pY[j]);
				if(dx!=.0) pans[j+nData*k] = atan(dy/dx);
				else if( dy != 0. ) {
					// x==0, y!=0
					// only positive infinity for pos correlated
					if(pCor) pans[j+nData*k] = M_PI_2;
					else pans[j+nData*k] = -M_PI_2;
				}
				else pans[j+nData*k] = NA_REAL; // dx==0, dy==0, set to NA
			}
			else pans[j+nData*k] = NA_REAL; // not upper triangle, set to NA
		}
	}
	// Wrap up and return results
	UNPROTECT(3);
	return(ans);
}
/*
  Function serves as interface between the DLL and R.
  All arguments are passed "by-reference", i.e. only pointers to variables created on R-side are passed. 
  pX            ... (double)  pointer indicating the 1st element of a double-array (vector) with X-coordinates of points
  pY            ... (double)  pointer indicating the 1st element of a double-array (vector) with Y-coordinates of points
  pNData        ... (int)     pointer to an integer storing the number of data-points
  pPosCor       ... (int)     pointer to an integer indicating whether X and Y is positively correlated , 1 = TRUE, 0 = FALSE
  pNBins        ... (int)     pointer to an integer specifying the number of bins used to characterize slope-values
  pSlope        ... (double)  pointer to a double, where the computed slope-value will be stored
  pQuantile     ... (double)  pointer to a double representing the 1-alpha/2 quantile of the standard normal distribution
  pSlopeLower   ... (double)  pointer to a double, where the computed lower confidence bound of the slope CI will be stored
  pSlopeUpper   ... (double)  pointer to a double, where the computed upper confidence bound of the slope CI will be stored
  pCIundefined   ... (int)    pointer to an integer indicating errors computing the confidence interval for the slope
*/
void PaBaLargeData( double * pX, double * pY, int * pNData, int * pPosCor, int * pNBins, double * pSlope,
                    double * pQuantile, double * pSlopeLower, double * pSlopeUpper, int * pCIundefined)
{
  bool LCLundef=false, UCLundef=false;                              /* indicate problems if set to TRUE */
  int *pBins;
  int nInd, nItems;
  int NPos, NPos2, NNeg, NNeg2, NRegular;                           /* these variables are only needed on C-side, not on R-side */
  int NAllItems, Index, Offset, nValIndex2, half;               
  double dConf;
  int Index1, Index2, LCLindex, LCLindex1, LCLindex2, UCLindex, UCLindex1, UCLindex2;     /* needed if Indices are even numbers */
  NPos = NPos2 = NNeg = NNeg2 = NRegular = NAllItems = Index = Offset = nValIndex2 = 0;
  pBins = (int *) calloc((*pNBins+1), sizeof(int));                 /* allocate memory for the binning-array */
  for(int i=0; i <=*pNBins; i++)
    pBins[i] = 0;
  FillBins( pBins, pPosCor, pX, pY, pNData, &NPos, &NPos2,          /* call workhorse-function */
            &NNeg, &NNeg2, &NRegular, &NAllItems, pNBins);
 
  if(*pPosCor == 1)                                                 /* compute Bin-index of the offsetted median */
    Offset = NNeg + NNeg2;
  else
    Offset = (-1) * (NPos + NPos2);
  /* determine slope of the regression line */
  nValIndex2 = NAllItems + Offset;
  half = (int)(nValIndex2 + 1)/2;                                   /* integer part, i.e. floor() */
  if(nValIndex2 % 2 == 0)                                           /* nValIndex is an even number */  
  {
    Index1 = IndexOf( half,     pBins, *pNBins);                    
    Index2 = IndexOf( half + 1, pBins, *pNBins);
    *pSlope = getMedianSlope(Index1, Index2, *pNBins);
  }
  else                                                              /* nValIndex2 is an odd number */ 
  {
    Index = IndexOf( half, pBins, *pNBins);                         
    *pSlope = Tan( ((double)(Index) / (*pNBins)) * PI - PI2);
  }
  /* calculate confidence interval for the slope, the CI for the intercept will be calculated on R-side */
  /* Lower CI-Bound */
	dConf = (*pQuantile) * sqrt( ((double)*pNData) *(((double)*pNData)-1.) * (2. * ((double)*pNData) + 5) / 18.);
  dConf = round(dConf);                                         /* round to the nearest integer */
  nInd     = (int) (NAllItems - dConf + Offset);                /* as in the exact algo */
  nItems   = (int)((nInd + 1)/2);                               /* !!! added + 1 as in the exact algo (nInd+1L)%/%2) */
  if(nItems >= 0)                                               /* otherwise the -INF init-value coming from R remains unchanged */
  {
    if( nInd % 2 == 0 )                                         /* nInd is an even number */
    {
      LCLindex1 = IndexOf(nItems,     pBins, *pNBins);
      LCLindex2 = IndexOf(nItems + 1, pBins, *pNBins);
		  if( (LCLindex1 >= 0) && (LCLindex2 >= 0) )                /* IndexOf returns -1 if the index found is equal to NBins + 1 */
      {
        *pSlopeLower = getMedianSlope(LCLindex1, LCLindex2, *pNBins);
      }
      else
		    LCLundef = true; 
    }
    else                                                        /* nInd is an odd number */
    {
	    LCLindex = IndexOf(nItems, pBins, *pNBins);
		  if(LCLindex >= 0)
        *pSlopeLower = Tan(((double)LCLindex/(*pNBins))*PI - PI2);
      else
		    LCLundef = true; 
    }
	}
  /* Upper CI-Bound */
  nInd   = (int) (NAllItems + dConf + Offset);
  nItems = (int) ((nInd + 1)/2);
	if(nItems <= NAllItems)                                 /* otherwise the INF init-value coming from R remains unchanged */
	{
    if( nInd % 2 == 0 )                                   /* nInd is an even number */
    {
      UCLindex1 = IndexOf(nItems,     pBins, *pNBins);
      UCLindex2 = IndexOf(nItems + 1, pBins, *pNBins);
	    if( (UCLindex1 >= 0) && (UCLindex2 >= 0) )          /* IndexOf returns -1 if the index found is equal to NBins + 1 */
      {
        *pSlopeUpper = getMedianSlope(UCLindex1, UCLindex2, *pNBins);
      }
		  else
        UCLundef = true;
    }
    else                                                  /* nInd is an odd number */
    {
	    UCLindex = IndexOf(nItems, pBins, *pNBins);
	    if(UCLindex >= 0)
        *pSlopeUpper = Tan((double)UCLindex/(*pNBins) * PI - PI2);
		  else
        UCLundef = true;
    }
	}
  if(LCLundef || UCLundef)        /* signal error(s) computing CI-bounds for slope */
  {
    *pCIundefined = 1;
  }
  free(pBins);                    
}
/*
  Iterates over all pairs of points, computes the slope and intercept of the
  emerging line and counts the number of occurences of a distinct slope falling
  into one of nSlope bins.
  pnSlots     ... (int)    pointer indicating the 1st element of a integer-array (vector) where the number of occurrences will be stored ... 
  posCor      ... (int)    pointer to an integer indicating whether X and Y is positively correlated (posCor == 1) or not (posCor != 1)
  pXvals      ... (double) pointer indicating the 1st element of a double-array (vector) with X-coordinates of points
  pYvals      ... (double) pointer indicating the 1st element of a double-array (vector) with Y-coordinates of points
  nData       ... (int)    pointer to an integer corresponding to the number of data pairs 
  nPos        ... (int)     pointer to an integer storing the number of slopes being >= 1  (pi/4)
  nPos2       ... (int)     pointer to an integer storing the number of slopes being > 1   (pi/4)
  nNeg        ... (int)     pointer to an integer storing the number of slopes being <= -1 (-pi/4)
  nNeg2       ... (int)     pointer to an integer storing the number of slopes being <  -1 (-pi/4)     
  nRegular    ... (int)     pointer to an integer counting regular slopes   
  nAllItems   ... (int)     pointer to an integer counting all items
  nSlots      ... (int)     pointer to an integer specifying the number of bins used to characterize slope-values
     
*/
void FillBins(  int * pnSlots, int * posCor, double * pXVals, double * pYVals,
                int * nData, int * nPos, int * nPos2, int * nNeg, int * nNeg2, 
                int * nRegular, int * nAllItems, int * nSlots)
{
	double dx,dy,phi;
  int Index;
  *nPos = *nPos2 = *nNeg = *nNeg2 = *nAllItems = *nRegular = 0;
	for(int j = 0; j < *nData; j++)
	{
		for(int k = j+1; k < *nData; k++)
		{
			dx = calcDiff(pXVals[k], pXVals[j]);
			dy = calcDiff(pYVals[k], pYVals[j]);
			if(dx != .0)                                        // avoiding division by zero              
			{
				phi = atan(dy/dx);
				Index = (int) (.5 + (*nSlots)*(phi+PI2)/PI);           // indices are computed for all angles within 1st and 2nd quadrant of a cartesian coordinate system
				pnSlots[Index]++;                                 // ... casting cuts decimal digits
				(*nAllItems)++;
        if( phi >= PI4 )
        {
          (*nPos)++;
          if( phi > PI4 )
          {
            (*nPos2)++;
          }
        }
        else if( phi <= -PI4 )
        {
          (*nNeg)++;
          if( phi < -PI4 )
          {
            (*nNeg2)++;
          }
        }
        else
        {
          (*nRegular)++;
        }
   		}
			else if( dy != 0. )                               // either -Inf or Inf (division by Zero)
			{
				if(*posCor == 1)                                // positively correlated
				{
					pnSlots[*nSlots]++;                           // slope equal to Infinity
					(*nAllItems)++;
					(*nPos)++;
					(*nPos2)++;
				}
				else                                            // negatively correlated
				{
					(*pnSlots)++;                                 // slope equal to -Infinity
					(*nAllItems)++;     
					(*nNeg)++;
					(*nNeg2)++;
				}
			}
      //			else dy == 0 
      //		-> ommit this value, because the atan (or slope) isn't defined !!!
		}
	}
}
/*
  compute index of that bin containing 'nItems' if all items of preceeding bins
  were cumulated, i.e. in which bin is the 'nItems'-th slope
  (taken from C.Kuhn C++-code)
  Note: argument 'nSlots' was added since there are no class attributes
        available in C.
*/
int IndexOf(int nItems, int * pnSlots, int nSlots)
{
	int nInd = 0;
  int Index = 0;
	for(int j = 0; j < nSlots + 1; j++)          // replaced < by <= nSlots, since there are N+1 slots actually
	{
		nInd += pnSlots[j];
		if(nInd >=  nItems)
    {
      Index = j;
			break;
    }
	}
	assert(Index < nSlots + 1);
	if(Index == nSlots + 1)
		return -1;
	return Index;
}
/* 
  Function computes tangent and returns INF or -INF constants in case
  of infinite or minus infinite value. 
*/
double Tan(double x)
{
	if(calcDiff(fabs(x), PI2) == 0.)
	{
		if(x > 0)
    {
			return INF;
    }
		else
    {
			return -INF;
    }
	}
	return tan(x);
}
/*
  Function computes slope values for two indices according to the number of bins
  and averages both values. This is used for computing the median slope out of an even
  number of slopes. Both slopes are avaraged on radians-scale, not until then, 
  the slope is computed back from this circle-measure!
*/
double getMedianSlope(int index1, int index2, int NBins)
{
  double Slope1 = ((double)index1 / (NBins)) * PI - PI2;
  double Slope2 = ((double)index2 / (NBins)) * PI - PI2;
  double Slope = (Slope1 + Slope2) / 2;
  Slope = Tan(Slope);
  return Slope;         
}
