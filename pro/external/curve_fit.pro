; $Id: curvefit.pro,v 1.20 1998/01/15 18:41:24 scottm Exp $
;
; Copyright (c) 1982-1998, Research Systems, Inc.  All rights reserved.
;	Unauthorized reproduction prohibited.
;
FUNCTION CURVE_FIT, x, y, weights, a, sigma, FUNCTION_NAME = Function_Name, $
                        ITMAX=itmax, ITER=iter, TOL=tol, CHI2=chi2, $
                        NODERIVATIVE=noderivative, CHISQ=chisq, _EXTRA=E
;+
; NAME:
;       CURVEFIT
;
; PURPOSE:
;       Non-linear least squares fit to a function of an arbitrary 
;       number of parameters.  The function may be any non-linear 
;       function.  If available, partial derivatives can be calculated by 
;       the user function, else this routine will estimate partial derivatives
;       with a forward difference approximation.
;
; CATEGORY:
;       E2 - Curve and Surface Fitting.
;
; CALLING SEQUENCE:
;       Result = CURVE_FIT(X, Y, Weights, A, SIGMA, FUNCTION_NAME = name, $
;                         ITMAX=ITMAX, ITER=ITER, TOL=TOL, /NODERIVATIVE)
;
; INPUTS:
;       X:  A row vector of independent variables.  This routine does
;           not manipulate or use values in X, it simply passes X
;           to the user-written function.
;
;       Y:  A row vector containing the dependent variable.
;
;  Weights:  A row vector of weights, the same length as Y.
;            For no weighting,
;                 Weights(i) = 1.0.
;            For instrumental (Gaussian) weighting,
;                 Weights(i)=1.0/sigma(i)^2
;            For statistical (Poisson)  weighting,
;                 Weights(i) = 1.0/y(i), etc.
;
;       A:  A vector, with as many elements as the number of terms, that 
;           contains the initial estimate for each parameter.  IF A is double-
;           precision, calculations are performed in double precision, 
;           otherwise they are performed in single precision. Fitted parameters
;           are returned in A.
;
; KEYWORDS:
;       FUNCTION_NAME:  The name of the function (actually, a procedure) to 
;       fit.  IF omitted, "FUNCT" is used. The procedure must be written as
;       described under RESTRICTIONS, below.
;
;       ITMAX:  Maximum number of iterations. Default = 20.
;       ITER:   The actual number of iterations which were performed
;       TOL:    The convergence tolerance. The routine returns when the
;               relative decrease in chi-squared is less than TOL in an 
;               interation. Default = 1.e-3.
;       CHI2:   The value of chi-squared on exit (obselete)
;     
;       CHISQ:   The value of reduced chi-squared on exit
;       NODERIVATIVE:   IF this keyword is set THEN the user procedure will not
;               be requested to provide partial derivatives. The partial
;               derivatives will be estimated in CURVEFIT using forward
;               differences. IF analytical derivatives are available they
;               should always be used.
;
; OUTPUTS:
;       Returns a vector of calculated values.
;       A:  A vector of parameters containing fit.
;
; OPTIONAL OUTPUT PARAMETERS:
;       Sigma:  A vector of standard deviations for the parameters in A.
;
; COMMON BLOCKS:
;       NONE.
;
; SIDE EFFECTS:
;       None.
;
; RESTRICTIONS:
;       The function to be fit must be defined and called FUNCT,
;       unless the FUNCTION_NAME keyword is supplied.  This function,
;       (actually written as a procedure) must accept values of
;       X (the independent variable), and A (the fitted function's
;       parameter values), and return F (the function's value at
;       X), and PDER (a 2D array of partial derivatives).
;       For an example, see FUNCT in the IDL User's Libaray.
;       A call to FUNCT is entered as:
;       FUNCT, X, A, F, PDER
; where:
;       X = Variable passed into CURVEFIT.  It is the job of the user-written
;           function to interpret this variable.
;       A = Vector of NTERMS function parameters, input.
;       F = Vector of NPOINT values of function, y(i) = funct(x), output.
;       PDER = Array, (NPOINT, NTERMS), of partial derivatives of funct.
;               PDER(I,J) = DErivative of function at ith point with
;               respect to jth parameter.  Optional output parameter.
;               PDER should not be calculated IF the parameter is not
;               supplied in call. IF the /NODERIVATIVE keyword is set in the
;               call to CURVEFIT THEN the user routine will never need to
;               calculate PDER.
;
; PROCEDURE:
;       Copied from "CURFIT", least squares fit to a non-linear
;       function, pages 237-239, Bevington, Data Reduction and Error
;       Analysis for the Physical Sciences.  This is adapted from:
;       Marquardt, "An Algorithm for Least-Squares Estimation of Nonlinear
;       Parameters", J. Soc. Ind. Appl. Math., Vol 11, no. 2, pp. 431-441,
;       June, 1963.
;
;       "This method is the Gradient-expansion algorithm which
;       combines the best features of the gradient search with
;       the method of linearizing the fitting function."
;
;       Iterations are performed until the chi square changes by
;       only TOL or until ITMAX iterations have been performed.
;
;       The initial guess of the parameter values should be
;       as close to the actual values as possible or the solution
;       may not converge.
;
; EXAMPLE:  Fit a function of the form f(x) = a * exp(b*x) + c to
;           sample pairs contained in x and y.
;           In this example, a=a(0), b=a(1) and c=a(2).
;           The partials are easily computed symbolicaly:
;           df/da = exp(b*x), df/db = a * x * exp(b*x), and df/dc = 1.0
;
;           Here is the user-written procedure to return F(x) and
;           the partials, given x:
;
;       pro gfunct, x, a, f, pder      ; Function + partials
;         bx = exp(a(1) * x)
;         f= a(0) * bx + a(2)         ;Evaluate the function
;         IF N_PARAMS() ge 4 THEN $   ;Return partials?
;         pder= [[bx], [a(0) * x * bx], [replicate(1.0, N_ELEMENTS(f))]]
;       end
;
;         x=findgen(10)                  ;Define indep & dep variables.
;         y=[12.0, 11.0,10.2,9.4,8.7,8.1,7.5,6.9,6.5,6.1]
;         Weights=1.0/y            ;Weights
;         a=[10.0,-0.1,2.0]        ;Initial guess
;         yfit=curvefit(x,y,Weights,a,sigma,function_name='gfunct')
;         print, 'Function parameters: ', a
;         print, yfit
;       end
;
; MODIFICATION HISTORY:
;       Written, DMS, RSI, September, 1982.
;       Does not iterate IF the first guess is good.  DMS, Oct, 1990.
;       Added CALL_PROCEDURE to make the function's name a parameter.
;              (Nov 1990)
;       12/14/92 - modified to reflect the changes in the 1991
;            edition of Bevington (eq. II-27) (jiy-suggested by CreaSo)
;       Mark Rivers, U of Chicago, Feb. 12, 1995
;           - Added following keywords: ITMAX, ITER, TOL, CHI2, NODERIVATIVE
;             These make the routine much more generally useful.
;           - Removed Oct. 1990 modification so the routine does one iteration
;             even IF first guess is good. Required to get meaningful output
;             for errors. 
;           - Added forward difference derivative calculations required for 
;             NODERIVATIVE keyword.
;           - Fixed a bug: PDER was passed to user's procedure on first call, 
;             but was not defined. Thus, user's procedure might not calculate
;             it, but the result was THEN used.
;
;      Steve Penton, RSI, June 1996.
;            - Changed SIGMAA to SIGMA to be consistant with other fitting 
;              routines.
;            - Changed CHI2 to CHISQ to be consistant with other fitting 
;              routines.
;            - Changed W to Weights to be consistant with other fitting 
;              routines.
;            _ Updated docs regarding weighing.
;
;	Vinay Kashyap, CfA, Dec 1998
;	     - Changed name to CURVE_FIT.
;	     - Added keyword _EXTRA=E to calling sequence and to all
;	       calls to CALL_PROCEDURE.
;           
;-
       ON_ERROR,2              ;Return to caller IF error

       ;Name of function to fit

       IF n_elements(function_name) LE 0 THEN function_name = "FUNCT"

       IF n_elements(tol) EQ 0 THEN tol = 1.e-3      ;Convergence tolerance
       IF n_elements(itmax) EQ 0 THEN itmax = 20     ;Maximum # iterations
       type = size(a)
       type = type[type[0]+1]
       double = type EQ 5

       IF (type ne 4) AND (type ne 5) THEN a = float(a)  ;Make params floating

       ; IF we will be estimating partial derivatives THEN compute machine
       ; precision

       IF keyword_set(NODERIVATIVE) THEN BEGIN
          res = machar(DOUBLE=double)
          eps = sqrt(res.eps)
       ENDIF

       nterms = n_elements(a)         ; # of parameters
       nfree = n_elements(y) - nterms ; Degrees of freedom

       IF nfree LE 0 THEN message, 'Curvefit - not enough data points.'

       flambda = 0.001                   ;Initial lambda
       diag = lindgen(nterms)*(nterms+1) ; Subscripts of diagonal elements

;      Define the partial derivative array

       IF double THEN pder = dblarr(n_elements(y), nterms) $
       ELSE pder = fltarr(n_elements(y), nterms)
;
       FOR iter = 1, itmax DO BEGIN      ; Iteration loop

;         Evaluate alpha and beta matricies.

          IF keyword_set(NODERIVATIVE) THEN BEGIN

;            Evaluate function and estimate partial derivatives
             CALL_PROCEDURE, Function_name, x, a, yfit, _EXTRA=E

             FOR term=0, nterms-1 DO BEGIN

                p = a       ; Copy current parameters

                ; Increment size for forward difference derivative
                inc = eps * abs(p[term])    
                IF (inc EQ 0.) THEN inc = eps
                p[term] = p[term] + inc
                CALL_PROCEDURE, function_name, x, p, yfit1, _EXTRA=E
                pder[0,term] = (yfit1-yfit)/inc

             ENDFOR
          ENDIF ELSE BEGIN

             ; The user's procedure will return partial derivatives
             call_procedure, function_name, x, a, yfit, pder , _EXTRA=E

          ENDELSE

          IF nterms EQ 1 THEN pder = reform(pder, n_elements(y), 1)

          beta = (y-yfit)*Weights # pder
          alpha = transpose(pder) # (Weights # (fltarr(nterms)+1)*pder)

          ; save current values of return parameters

          sigma1 = sqrt( 1.0 / alpha[diag] )           ; Current sigma.
          sigma  = sigma1

          chisq1 = total(Weights*(y-yfit)^2)/nfree     ; Current chi squared.
          chisq = chisq1

          yfit1 = yfit                                 

          done_early = chisq1 LT total(abs(y))/1e7/NFREE 
          IF done_early THEN GOTO, done

          c = sqrt(alpha[diag])
          c = c # c

          lambdaCount = 0

          REPEAT BEGIN

             lambdaCount = lambdaCount + 1

             ; Normalize alpha to have unit diagonal.

             array = alpha / c

             ; Augment the diagonal.

             array[diag] = array[diag]*(1.+flambda) 

             ; Invert modified curvature matrix to find new parameters.

             IF n_elements(array) EQ 1 THEN array = (1.0 / array) $
             ELSE array = invert(array)

             b = a + array/c # transpose(beta)          ; New params

             call_procedure, function_name, x, b, yfit, _EXTRA=E  ; Evaluate function
             chisq = total(Weights*(y-yfit)^2)/nfree    ; New chisq
             sigma = sqrt(array[diag]/alpha[diag])      ; New sigma

             IF (finite(chisq) EQ 0) OR $
                  (lambdaCount GT 30 AND chisq GE chisq1) THEN BEGIN

                ; Reject changes made this iteration, use old values.

                yfit  = yfit1
                sigma = sigma1
                chisq = chisq1

                message, 'Failed to converge', /INFORMATIONAL

                GOTO, done 

             ENDIF             

             flambda = flambda*10.               ; Assume fit got worse

          ENDREP UNTIL chisq LE chisq1

          flambda = flambda/100.  

          a=b                                    ; Save new parameter estimate.

          IF ((chisq1-chisq)/chisq1) LE tol THEN GOTO,done   ;Finished?
       ENDFOR                        ;iteration loop
;
       MESSAGE, 'Failed to converge', /INFORMATIONAL
;
done:  chi2 = chisq         ; Return chi-squared (chi2 obsolete-still works)
       IF done_early THEN iter = iter - 1
       return,yfit          ; return result
END
