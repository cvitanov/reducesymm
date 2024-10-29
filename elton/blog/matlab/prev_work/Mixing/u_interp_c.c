#include "mex.h"
#include <math.h>
void mexFunction(int nlhs, mxArray *plhs[], 
    int nrhs, const mxArray *prhs[])
{
    /* There are supposed to be 6 input items (nrhs = 6), which we do not check, assume it is called correctly */
    /* Inputs are pointers to: dummy t-value,(x,y,z) vector, Lx,Lz,(nx,ny,nz) vector, and v_grid (a 4-d array) */
    /* There is one ouput item (nlhs = 1).  It is (u,v,w) vector, will show up as a column vector in matlab */
    double *ptr_x_vect,*ptr_nx_vect,*ptr_vgrid, *ptr_u_output;
    double x,y,z,Lx,Lz,dx,dy,dz,gx,gy,gz,cgx,cgy,cgz,div;
    double a000,a001,a010,a011,a100,a101,a110,a111;
    int i,j,k,j_mult,k_mult,nx,ny,nz,offst,offst001,offst010,offst011,offst100,offst101,offst110,offst111;
    int n;
    
    /* get x,y,z that was input */
    ptr_x_vect = mxGetPr(prhs[1]); /* prhs[0] points to dummy t value which is ignored */
    x = ptr_x_vect[0];
    y = ptr_x_vect[1] + 1; /* change y to range [0,2] rather than [-1,1] for our purposes here */
    z = ptr_x_vect[2];
    //mexPrintf("x = %lf y = %lf z = %lf\n",x,y,z);
    
    /* create an output column vector and a pointer to it */
    plhs[0] = mxCreateDoubleMatrix(3, 1, mxREAL); /* three rows, one column */
    ptr_u_output = mxGetPr(plhs[0]);
       
    /* Check y. If y outside the box, just return with a fixed u-vector that runs along the bdry */
    if (y > 2 || y < 0)
    {
        if(y > 2)
            ptr_u_output[0] = 1.0; /* run forward at top of box */
        else
            ptr_u_output[0] = -1.0; /* run backwards at bottom of box */            
        ptr_u_output[1] = 0.0;
        ptr_u_output[2] = 0.0;
        return;
    }
    /* Get Lx and Lz values */
    Lx = (mxGetPr(prhs[2]))[0];
    Lz = (mxGetPr(prhs[3]))[0];    
    
    /* Check x and z. If x or z is outside of the periodic box, shift it back inside to calculate the velocity field.
       Probably x is not greater than 2*Lx, so we will try that first, then in unlikely case it is > Lx, do something more.
        Similarly for x < 0, and for z. */
    if( x > Lx )
    {
        x = x - Lx;
        if (x > Lx)
        {           
            div = (double)((int)(x/Lx));  /* get the floor of x/Lx */
            x = x - div*Lx;
        }
    }
    else if (x < 0 )
    {
        x = x + Lx;
        if (x < 0)
        {
            div = (double)((int)abs(x/Lx));  /* get the floor of abs(x/Lx) */
            x = x + (1.0 + div)*Lx;
        }
    }
    if( z > Lz )
    {
        z = z - Lz;
        if (z > Lz)
        {           
            div = (double)((int)(z/Lz));  /* get the floor of z/Lz */
            z = z - div*Lz;
        }
    }
    else if (z < 0 )
    {
        z = z + Lz;
        if (z < 0)
        {
            div = (double)((int)abs(z/Lz));  /* get the floor of abs(z/Lz) */
            z = z + (1.0 + div)*Lz;
        }
    }
    
    
    //mexPrintf("Lx = %lf Lz = %lf\n",Lx,Lz);
    ptr_nx_vect = mxGetPr(prhs[4]);
    nx = (int)(ptr_nx_vect[0]);
    ny = (int)(ptr_nx_vect[1]);
    nz = (int)(ptr_nx_vect[2]);
    //mexPrintf("nx = %d ny = %d nz = %d\n",nx,ny,nz);
    ptr_vgrid = mxGetPr(prhs[5]);
    
    /* The following things (dx, dy,dz, j_mult, k_mult, and offstxxx's, would be computed just once at beginning in the real version */
    dx = Lx/(double)nx;
    dy = 2.0/(double)ny;
    dz = Lz/(double)nz;
    //mexPrintf("dx = %lf dy = %lf dz = %lf\n",dx,dy,dz);
    j_mult = 3*(nx+1);
    k_mult = j_mult*(ny+1);
    /* compute offsets from bottom-left-lower corner of grid box to other corners */
    offst001 = k_mult;
    offst010 = j_mult;
    offst011 = j_mult + k_mult;
    offst100 = 3;
    offst101 = 3 + k_mult;
    offst110 = 3 + j_mult;
    offst111 = 3 + j_mult + k_mult;
    
    gx = x/dx;
    gy = y/dy;
    gz = z/dz;
    i = (int)gx; 
    j = (int)gy;
    k = (int)gz; 
    //mexPrintf("i = %d j = %d k = %d\n",i,j,k);
    gx = gx-(double)i; /* fractional part */
    gy = gy-(double)j; /* fractional part */
    gz = gz-(double)k; /* fractional part */
    cgx = 1.0 - gx;
    cgy = 1.0 - gy;
    cgz = 1.0 - gz;
    
    /* get the 8 coefficents for interpolating.  These are the same for all three components of output.
       Do this with 12 multiplications. */
    a000 = cgx*cgy;
    a001 = a000*gz; /* so this is cgx*cgy*gz */
    a000 = a000*cgz; /* so this is cgx*cgy*cgz.  We managed to get two coeffs with three multiplications instead of four. */
    a010 = cgx*gy;
    a011 = a010*gz;
    a010 = a010*cgz;
    a100 = gx*cgy;
    a101 = a100*gz;
    a100 = a100*cgz;
    a110 = gx*gy;
    a111 = a110*gz;
    a110 = a110*cgz;
    
    /* It is at least as efficient for us to compute the offset into the 4-d array ourselves, than use the 4-d indexing that c would do itself,
       because we have precomputed the multipliers. */
    offst = i*3 + j*j_mult + k*k_mult; /* bottom-left-lower corner of grid_box */
    
    ptr_u_output[0] = a000*ptr_vgrid[offst] + a001*ptr_vgrid[offst+offst001] + a010*ptr_vgrid[offst+offst010] + a011*ptr_vgrid[offst+offst011] +
                      a100*ptr_vgrid[offst+offst100] + a101*ptr_vgrid[offst+offst101] + a110*ptr_vgrid[offst+offst110] + a111*ptr_vgrid[offst+offst111];
    offst = offst + 1; /* move to next output component. */
    ptr_u_output[1] = a000*ptr_vgrid[offst] + a001*ptr_vgrid[offst+offst001] + a010*ptr_vgrid[offst+offst010] + a011*ptr_vgrid[offst+offst011] +
                      a100*ptr_vgrid[offst+offst100] + a101*ptr_vgrid[offst+offst101] + a110*ptr_vgrid[offst+offst110] + a111*ptr_vgrid[offst+offst111];;
    offst = offst + 1; /* move to next output component. */
    ptr_u_output[2] = a000*ptr_vgrid[offst] + a001*ptr_vgrid[offst+offst001] + a010*ptr_vgrid[offst+offst010] + a011*ptr_vgrid[offst+offst011] +
                      a100*ptr_vgrid[offst+offst100] + a101*ptr_vgrid[offst+offst101] + a110*ptr_vgrid[offst+offst110] + a111*ptr_vgrid[offst+offst111];;
/* Note this has taken a total of 36 double-precision multiplications. */                      
    
}

    
    
    
    