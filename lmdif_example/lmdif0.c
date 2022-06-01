/* lmdif0.c -- driver for lmdif */
#include <math.h>
#if !defined(__APPLE__)
#include <malloc.h>
#endif
#include <stdlib.h>
#include "cminpak.h"

int lmdif0(void fcn(int,int,double [],double [],int *),int m, int n,double x[],int msk[],
     double fvec[],double tol,int *info,int *nfev)
{
    int j,maxfev,mode;
    int *ipvt;
    double ftol,xtol,gtol,epsfcn,factor;
    double *diag,**fjac,*qtf,*wa1,*wa2,*wa3,*wa4;

/* Check input parameters */
   if (n <= 0 || m < n || tol < 0.0) {
        *info = 0;
        return(1);
   }
/* Allocate memory for working arrays. */
    ipvt = (int *)calloc(n,sizeof(int));
    diag = (double *)calloc(n,sizeof(double));
    qtf = (double *)calloc(n,sizeof(double));
    wa1 = (double *)calloc(n,sizeof(double));
    wa2 = (double *)calloc(n,sizeof(double));
    wa3 = (double *)calloc(n,sizeof(double));
    wa4 = (double *)calloc(m,sizeof(double));


/* Create 2d matrix for Jacobian */
    fjac = (double **)calloc(n,sizeof(double *));
    for (j=0;j<n;j++)
        fjac[j] = (double *)calloc(m,sizeof(double));

/* Set convergence tolerances */
    ftol = tol;
    xtol = tol;
    gtol = 0.0;

    maxfev = 800;
    epsfcn = 0.0;
    mode = 1;
    factor = 10;
    *nfev = 0;

    lmdif(fcn,m,n,x,msk,fvec,ftol,xtol,gtol,maxfev,epsfcn,diag,mode,
        factor,info,nfev,fjac,ipvt,qtf,wa1,wa2,wa3,wa4);

    if (*info == 8) *info = 4;
    for (j=0;j<n;j++)
        free(fjac[j]);
    free(fjac);
    free(wa4);
    free(wa3);
    free(wa2);
    free(wa1);
    free(qtf);
    free(diag);
    free(ipvt);

    printf("[x[0]]     %.5e\n",x[0]);
    printf("[x[1]]     %.5e\n",x[1]);
    printf("[x[2]]     %.5e\n",x[2]);
    for (j=0;j<m;j++)
        printf("[fvec]     %.5e\n",fvec[j]);

    return(0);
}
    



    
