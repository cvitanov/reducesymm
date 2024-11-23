#include <math.h>
#include "mex.h"

void mexFunction(int nlhs, mxArray *plhs[],
    int nrhs, const mxArray *prhs[]) 
{
    /* There must be 4 arguments.  The first is pointer to x, the second is pointer to deltax, the third and fourth are the arrays of doubles 
       to be interpolated.  It is assumed that arrays have at least x/deltax + 1 elements. */
    
    double x,dx,r,cr;
    double *yinput, *youtput;
    int i;
    
    x = *(mxGetPr(prhs[0]));  /* get the input x value */
    dx = *(mxGetPr(prhs[1]));  /* get the input deltax value */
    r = x/dx;
    i = (int)r; /* this is floor of x/dx, that is, the integer part */
    r = r - (double)i;  /* this is fractional part */
    cr = 1.0 - r;
    plhs[0] = mxCreateDoubleMatrix(1, 2, mxREAL);/* Create an mxArray for the output data, which is two interpolated values. */    
    youtput = mxGetPr(plhs[0]);    
    yinput = mxGetPr(prhs[2]); /* pointer to first array to interpolate */
    youtput[0] = cr*yinput[i] + r*yinput[i+1];  /* first interpolated value */
    yinput = mxGetPr(prhs[3]); /* pointer to second array to interpolate */
    youtput[1] = cr*yinput[i] + r*yinput[i+1];  /* second interpolated value */
    
} 

