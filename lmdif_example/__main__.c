#include <math.h>
#if !defined(__APPLE__)
#include <malloc.h>
#endif
#include <stdlib.h>
#include <stdio.h>
#include "cminpak.h"
#include "enorm.c"
#include "fdjac2.c"
#include "lmdif.c"
#include "lmdif0.c"
#include "lmpar.c"
#include "qrfac.c"
#include "qrsolv.c"

/*  reference : https://www.math.utah.edu/software/minpack/minpack/lmdif1.html
 3. Parameters.

       Parameters designated as input parameters must be specified on
       entry to LMDIF1 and are not changed on exit, while parameters
       designated as output parameters need not be specified on entry
       and are set to appropriate values on exit from LMDIF1.

       FCN is the name of the user-supplied subroutine which calculates
         the functions.  FCN must be declared in an EXTERNAL statement
         in the user calling program, and should be written as follows.

         SUBROUTINE FCN(M,N,X,FVEC,IFLAG)
         INTEGER M,N,IFLAG
         DOUBLE PRECISION X(N),FVEC(M)
         ----------
         CALCULATE THE FUNCTIONS AT X AND
         RETURN THIS VECTOR IN FVEC.
         ----------
         RETURN
         END

         The value of IFLAG should not be changed by FCN unless the
         user wants to terminate execution of LMDIF1.  In this case set


                                                                 Page 2

         IFLAG to a negative integer.

       M is a positive integer input variable set to the number of
         functions.

       N is a positive integer input variable set to the number of
         variables.  N must not exceed M.

       X is an array of length N.  On input X must contain an initial
         estimate of the solution vector.  On output X contains the
         final estimate of the solution vector.

       FVEC is an output array of length M which contains the functions
         evaluated at the output X.

       TOL is a nonnegative input variable.  Termination occurs when
         the algorithm estimates either that the relative error in the
         sum of squares is at most TOL or that the relative error
         between X and the solution is at most TOL.  Section 4 contains
         more details about TOL.

       INFO is an integer output variable.  If the user has terminated
         execution, INFO is set to the (negative) value of IFLAG.  See
         description of FCN.  Otherwise, INFO is set as follows.

         INFO = 0  Improper input parameters.

         INFO = 1  Algorithm estimates that the relative error in the
                   sum of squares is at most TOL.

         INFO = 2  Algorithm estimates that the relative error between
                   X and the solution is at most TOL.

         INFO = 3  Conditions for INFO = 1 and INFO = 2 both hold.

         INFO = 4  FVEC is orthogonal to the columns of the Jacobian to
                   machine precision.

         INFO = 5  Number of calls to FCN has reached or exceeded
                   200*(N+1).

         INFO = 6  TOL is too small.  No further reduction in the sum
                   of squares is possible.

         INFO = 7  TOL is too small.  No further improvement in the
                   approximate solution X is possible.

         Sections 4 and 5 contain more details about INFO.

       IWA is an integer work array of length N.

       WA is a work array of length LWA.

       LWA is a positive integer input variable not less than
*/

//■■■■■■■■■■■■■■■■■■■■■■■■■
double sol = 0;
double fric_pred =0;
double resi = 0;
//■■■■■■■■■■■■■■■■■■■■■■■■■
/* independent variable */
extern double data[];
/*  measured value of dependent quantity */
extern double y[];
extern double dy[];
extern double fvecs[];
extern double x[];
//■■■■■■■■■■■■■■■■■■■■■■■■■
// data reference : https://heasarc.gsfc.nasa.gov/ftools/caldb/help/HDmpfit.html
double ey[10];
double data[]  = {-1.7237128E+00,1.8712276E+00,-9.6608055E-01,
                  -2.8394297E-01,1.3416969E+00,1.3757038E+00,
                  -1.3703436E+00,4.2581975E-02,-1.4970151E-01,
                   8.2065094E-01};
/*  measured value of dependent quantity */
double y[]     = {1.9000429E-01,6.5807428E+00,1.4582725E+00,
                  2.7270851E+00,5.5969253E+00,5.6249280E+00,
                  0.787615,3.2599759E+00,2.9771762E+00,
                  4.5936475E+00};
double fvecs[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
double x[3]    = {1.0, 1.0, 1.0};           /* Initial conditions */
//■■■■■■■■■■■■■■■■■■■■■■■■■
int msks        = -1;
double tols     = 1E-10;
int infoes      = 1;
int nfevs       = 0;
double results  = 0;
double dpmpar[] = { 2.220446049250313e-16,
                    2.225073858507201e-308,
                    1.79769313486231e+308  };
//■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■
void residual(int m, int n, double *x, double *fvecs, int *info ){
    
    int i;
    double temp;

    for (i=0; i<m; i++) {
        temp = x[0] + x[1]* data[i]+ x[2]* data[i]* data[i]  ;     /* Linear fit function; note f = a + b*x */
        fvecs[i] = (y[i] - temp)/ey[i];
        //fvecs[i] = (y[i] - temp);
    }

}
//■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■
int main() {

    int i;
    for (i=0; i<10; i++) ey[i] = 1.0;   /* Data errors factor */ 

    results = lmdif0(residual, 10, 3 , x,  &msks, fvecs, tols, &infoes, &nfevs);

    return 0;

}