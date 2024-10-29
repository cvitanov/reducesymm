#include "mex.h"
#include <math.h>
typedef struct 
{
    double dr[3];
    double Lx;
    double Lz;
    double *ptr_vgrid;
    int offsts[8];
    int j_mult;
    int k_mult;
} INTERP_PARMS;
void interp_c(double *kk, double *r, INTERP_PARMS *p);

void mexFunction(int nlhs, mxArray *plhs[], 
    int nrhs, const mxArray *prhs[])
{
    /* There are supposed to be 3 input items (nrhs = 3), which we do not check, assume it is called correctly */
    /* Inputs are pointers to: initial(rx,ry,rz) vector, parameters, and v_grid (a 4-d array).  The parameters are:
       h,np,dx,dy,dz,Lx,Lz,j_mult,k_mult where
       h is step-size, np = num.pts, not counting initial pt; the other parameters are used in interpolating within a grid box. */
    /* There is one ouput item (nlhs = 1).  It is a 3 by (np+1) matrix */
    INTERP_PARMS ip;
    double *ptr_r0_vect,*ptr_parms_vect,*ptr_r_matrix;
    double r[3],k1[3],k2[3],k3[3],k4[3],h,h2,h6;
    int n,np;
    
    /* get r0 ptr that was input */
    ptr_r0_vect = mxGetPr(prhs[0]);
    /* get input parms ptr */
    ptr_parms_vect = mxGetPr(prhs[1]);
    h = ptr_parms_vect[0];
    h2 = h*.5;
    h6 = h/6.0;
    np = (int)(ptr_parms_vect[1]);
    ip.dr[0] = ptr_parms_vect[2];
    ip.dr[1] = ptr_parms_vect[3];
    ip.dr[2] = ptr_parms_vect[4];
    ip.Lx = ptr_parms_vect[5];
    ip.Lz = ptr_parms_vect[6];
    ip.j_mult = (int)(ptr_parms_vect[7]);
    ip.k_mult = (int)(ptr_parms_vect[8]);
    ip.offsts[0] = 0; //000 corner
    ip.offsts[1] = ip.k_mult; //001 corner
    ip.offsts[2] = ip.j_mult; //010 corner
    ip.offsts[3] = ip.j_mult + ip.k_mult; //011 corner
    ip.offsts[4] = 3; //100 corner
    ip.offsts[5] = 3 + ip.k_mult; //101 corner
    ip.offsts[6] = 3 + ip.j_mult; //110 corner
    ip.offsts[7] = 3 + ip.j_mult + ip.k_mult; //111 corner
    /* get v_grid ptr */
    ip.ptr_vgrid = mxGetPr(prhs[2]);
       
    /* create an output matrix and a pointer to it */
    plhs[0] = mxCreateDoubleMatrix(3, np+1, mxREAL); /* three rows, np+1 columns */
    ptr_r_matrix = mxGetPr(plhs[0]);
    /* put first column to be the initial point */
    ptr_r_matrix[0] = ptr_r0_vect[0];
    ptr_r_matrix[1] = ptr_r0_vect[1];
    ptr_r_matrix[2] = ptr_r0_vect[2];
    //mexPrintf("j_mult = %d dr[0] = %g\n",ip.j_mult,ip.dr[0]);
    for(n = 0;n < np;++n)  //np
    {
        interp_c(k1,ptr_r_matrix,&ip);  /* returns answer in k1 */
        r[0] = ptr_r_matrix[0] + h2*k1[0];
        r[1] = ptr_r_matrix[1] + h2*k1[1];
        r[2] = ptr_r_matrix[2] + h2*k1[2];
        
        interp_c(k2,r,&ip);  /* returns answer in k2 */
        r[0] = ptr_r_matrix[0] + h2*k2[0];
        r[1] = ptr_r_matrix[1] + h2*k2[1];
        r[2] = ptr_r_matrix[2] + h2*k2[2];
        
        interp_c(k3,r,&ip);  /* returns answer in k3 */
        r[0] = ptr_r_matrix[0] + h*k3[0];
        r[1] = ptr_r_matrix[1] + h*k3[1];
        r[2] = ptr_r_matrix[2] + h*k3[2];
    
        interp_c(k4,r,&ip);  /* returns answer in k4 */
     
        ptr_r_matrix[3] = ptr_r_matrix[0] + h6*(k1[0] + 2.0*k2[0] + 2.0*k3[0] + k4[0]);
        ptr_r_matrix[4] = ptr_r_matrix[1] + h6*(k1[1] + 2.0*k2[1] + 2.0*k3[1] + k4[1]);
        ptr_r_matrix[5] = ptr_r_matrix[2] + h6*(k1[2] + 2.0*k2[2] + 2.0*k3[2] + k4[2]);
        ptr_r_matrix += 3; /* move to next column of output matrix */
    }
}

void interp_c(double *kk, double *r, INTERP_PARMS *p)
{   /* kk is size-3 array of outputs, r is size-three array of inputs, p is ptr to structure containing needed parameters for interpolation */
    
    double x,y,z,div,gx,gy,gz,cgx,cgy,cgz;
    double a[8];
    double *pr_vgrid;
    int i,j,k;
    
    x = r[0];
    y = r[1]+1;
    z = r[2];
    
    /* Check y. If y outside the box, just return with a fixed vector that runs along the bdry */
    if (y > 2 || y < 0)
    {
        if(y > 2)
            kk[0] = 1.0; /* run forward at top of box */
        else
            kk[0] = -1.0; /* run backwards at bottom of box */            
        kk[1] = 0.0;
        kk[2] = 0.0;
        return;
    }
    /* Check x and z. If x or z is outside of the periodic box, shift it back inside to calculate the velocity field.
       Probably x is not greater than 2*Lx, so we will try that first, then in unlikely case it is > Lx, do something more.
        Similarly for x < 0, and for z. */
    if( x >= p->Lx )
    {
        x = x - p->Lx;
        if (x >= p->Lx)
        {           
            div = (double)((int)(x/p->Lx));  /* get the floor of x/Lx */
            x = x - div*p->Lx;
        }
    }
    else if (x < 0 )
    {
        x = x + p->Lx;
        if (x < 0)
        {
            div = (double)((int)abs(x/p->Lx));  /* get the floor of abs(x/Lx) */
            x = x + (1.0 + div)*p->Lx;
        }
    }
    if( z >= p->Lz )
    {
        z = z - p->Lz;
        if (z >= p->Lz)
        {           
            div = (double)((int)(z/p->Lz));  /* get the floor of z/Lz */
            z = z - div*p->Lz;
        }
    }
    else if (z < 0 )
    {
        z = z + p->Lz;
        if (z < 0)
        {
            div = (double)((int)abs(z/p->Lz));  /* get the floor of abs(z/Lz) */
            z = z + (1.0 + div)*p->Lz;
        }
    }
    
    gx = x/p->dr[0];
    gy = y/p->dr[1];
    gz = z/p->dr[2];
    i = (int)gx; 
    j = (int)gy;
    k = (int)gz; 
    //mexPrintf("i = %d j = %d k = %d\n",i,j,k);
    gx -= (double)i; /* fractional part */
    gy -= (double)j; /* fractional part */
    gz -= (double)k; /* fractional part */
    cgx = 1.0 - gx;
    cgy = 1.0 - gy;
    cgz = 1.0 - gz;
    
    /* get the 8 coefficents for interpolating.  These are the same for all three components of output.
       Do this with 12 multiplications. */
    a[0] = cgx*cgy; //000
    a[1] = a[0]*gz; //001 /* so this is cgx*cgy*gz */
    a[0] *= cgz; /* so this is cgx*cgy*cgz.  We managed to get two coeffs with three multiplications instead of four. */
    a[2] = cgx*gy; //010
    a[3] = a[2]*gz; //011
    a[2] *= cgz;
    a[4] = gx*cgy;//100
    a[5] = a[4]*gz;//101
    a[4] *= cgz;
    a[6] = gx*gy;//110
    a[7] = a[6]*gz;//111
    a[6] *= cgz;
    
    //mexPrintf("gx = %g gy = %g gz = %g\n",gx,gy,gz);    
    
    /* It is at least as efficient for us to compute the offset into the 4-d array ourselves, than use the 4-d indexing that c would do itself,
       because we have precomputed the multipliers. */
    
    pr_vgrid = p->ptr_vgrid + i*3 + j*p->j_mult + k*p->k_mult; /* bottom-left-lower corner of grid_box */
    
#if 0
    /* NOTE: Here is an experiment with repacing the offsts[] with immediate numbers and recompiling every time they change
        since these don't change once the number of gridpoints has been set.
        The overall speed increase for trajectory_c_runge seems about 1 to 2% (not much). */

    kk[0] = a[0]*pr_vgrid[0] + a[1]*pr_vgrid[45360] + a[2]*pr_vgrid[432] + 
            a[3]*pr_vgrid[432+45360] + a[4]*pr_vgrid[3] + a[5]*pr_vgrid[3+45360] +
            a[6]*pr_vgrid[432+45360] + a[7]*pr_vgrid[3+432+45360];
    pr_vgrid += 1; /* move to next output component. */
    kk[1] = a[0]*pr_vgrid[0] + a[1]*pr_vgrid[45360] + a[2]*pr_vgrid[432] + 
            a[3]*pr_vgrid[432+45360] + a[4]*pr_vgrid[3] + a[5]*pr_vgrid[3+45360] +
            a[6]*pr_vgrid[432+45360] + a[7]*pr_vgrid[3+432+45360];
    pr_vgrid += 1; /* move to next output component. */
    kk[2] = a[0]*pr_vgrid[0] + a[1]*pr_vgrid[45360] + a[2]*pr_vgrid[432] + 
            a[3]*pr_vgrid[432+45360] + a[4]*pr_vgrid[3] + a[5]*pr_vgrid[3+45360] +
            a[6]*pr_vgrid[432+45360] + a[7]*pr_vgrid[3+432+45360];
#else
    
    kk[0] = a[0]*pr_vgrid[0] + a[1]*pr_vgrid[p->offsts[1]] + a[2]*pr_vgrid[p->offsts[2]] + 
            a[3]*pr_vgrid[p->offsts[3]] + a[4]*pr_vgrid[3] + a[5]*pr_vgrid[p->offsts[5]] +
            a[6]*pr_vgrid[p->offsts[6]] + a[7]*pr_vgrid[p->offsts[7]];
    pr_vgrid += 1; /* move to next output component. */
    kk[1] = a[0]*pr_vgrid[0] + a[1]*pr_vgrid[p->offsts[1]] + a[2]*pr_vgrid[p->offsts[2]] + 
            a[3]*pr_vgrid[p->offsts[3]] + a[4]*pr_vgrid[3] + a[5]*pr_vgrid[p->offsts[5]] +
            a[6]*pr_vgrid[p->offsts[6]] + a[7]*pr_vgrid[p->offsts[7]];
    pr_vgrid += 1; /* move to next output component. */
    kk[2] = a[0]*pr_vgrid[0] + a[1]*pr_vgrid[p->offsts[1]] + a[2]*pr_vgrid[p->offsts[2]] + 
            a[3]*pr_vgrid[p->offsts[3]] + a[4]*pr_vgrid[3] + a[5]*pr_vgrid[p->offsts[5]] +
            a[6]*pr_vgrid[p->offsts[6]] + a[7]*pr_vgrid[p->offsts[7]];
#endif                       
    /* Note this has taken a total of 36 double-precision multiplications. */                      
}
    

    
    
    
    