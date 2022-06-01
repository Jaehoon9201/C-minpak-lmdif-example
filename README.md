# MINPACK LMDIF example in C
**Example for MINPACK subroutine LMDIF in C language**

-----------------
**[Reference1](https://www.math.utah.edu/software/minpack/minpack/lmdif1.html)**  (https://www.math.utah.edu/software/minpack/minpack/lmdif1.html)

**[Reference2](https://heasarc.gsfc.nasa.gov/ftools/caldb/help/HDmpfit.html)**  (https://heasarc.gsfc.nasa.gov/ftools/caldb/help/HDmpfit.html)


----------------

```

Description of Parameters

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
```

-----------------

## Results in VS code(C)

![image](https://user-images.githubusercontent.com/71545160/171332707-4c8bd1f3-3066-444f-948c-aa6d426d6bf9.png)


## Evaluation with visualization in MATLAB

```matlab
data  = [-1.7237128,1.8712276,-0.96608055,...
      -0.28394297,1.3416969,1.3757038,...
      -1.3703436,0.042581975,-0.14970151,...
       0.82065094];

y     = [0.19000429,6.5807428,1.4582725,...
      2.7270851,5.5969253,5.6249280,...
      0.787615,3.2599759,2.9771762,...
      4.5936475];


error_factor = 1;
a0 = 3.18970;
a1 = 1 ;
a2 = 0.0149020 ; 

opt = a0 + a1 .* data + a2.* data.* data ;

plot(y, 'k'); hold on;
plot(opt, 'b');
plot(error_factor.*(y-opt), 'r--'); grid on;
legend('real vals', 'opt vals', 'residuals');
```

![image](https://user-images.githubusercontent.com/71545160/171334344-57441d74-1d62-4ed9-956d-04986f6f22b8.png)
