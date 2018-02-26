pro levmarq,x,y,a,chisq,ifit=ifit,alamda=alamda,jumpup=jumpup,jumpdn=jumpdn,$
	ochisq=ochisq,alpha=alpha,covar=covar,bvec=bvec,svdthr=svdthr,$
	trial=trial,yfunc=yfunc,tiptoe=tiptoe,verbose=verbose, _extra=e
;+
;procedure	levmarq
;	Levenberg-Marquardt method, attempting to reduce the value chi^2 of a
;	fit between a set of data points Y(X) with individual std.dev.s SIG,
;	and a nonlinear function dependent of A parameters.  The input array
;	IFIT indicates by non-zero entries those components of A that should
;	be fitted for, and by 0 entries those components that should be held
;	fixed at their input values.  The program returns current best-fit
;	values for the parameters A and the value of CHISQ.  The arrays
;	COVAR and ALPHA are used as working space during most iterations.
;	Supply a routine FUNCS(x,a,yfit,dYdA) that evaluates the fitting
;	function YFIT, and its derivatives dY/dA wrt fitting parameters A@X.
;	On the first call, provide an initial guess for parameters A, and
;	set ALAMDA < 0 for initialization.  If a step succeeds, CHISQ becomes
;	smaller and ALAMDA decreases by a factor 10.  If a step fails, ALAMDA
;	grows by a factor 3.  Call this routine repeatedly until convergence
;	is achieved.  Then, make one final call with ALAMDA=0, so that COVAR
;	returns the covariance matrix, and ALPHA the curvature matrix.  Those
;	parameters held fixed will return zero covariances.
;
;syntax
;	levmarq,x,y,a,chisq,funcs=funcs,sig=sig,ifit=ifit,alamda=alamda,$
;	jumpup=jumpup,jumpdn=jumpdn,ochisq=ochisq,alpha=alpha,covar=covar,$
;	bvec=bvec,svdthr=svdthr,trial=trial,yfunc=yfunc,/tiptoe,$
;	funcs=funcs,sig=sig,ties=ties, FUNCS_KEYWORDS
;
;parameters
;	x	[INPUT; required]
;	y	[INPUT; required] size must match that of X
;	a	[INPUT; required]
;	chisq	[OUTPUT; required] the chi-sq statistic denoting degree of
;		agreement between model and data
;
;keywords
;	ifit	[INPUT] integer array of same size as A, with 0's indicating
;		frozen parameters and 1's indicating thawed parameters
;		* if IFIT does not match A for any reason, all parameters
;		  are assumed to be thawed
;	jumpup	[INPUT; default=3] factor by which to increase ALAMDA if
;		current trial fails
;	jumpdn	[INPUT; default=0.1] factor by which to decrease ALAMDA if
;		current trial succeeds
;	alamda	[I/O] if not defined on input or is -ve, initializes the
;		procedure.
;	ochisq	[OUTPUT] the old value of the chi^2
;	alpha	[OUTPUT] if ALAMDA=0, the curvature matrix
;	covar	[OUTPUT] if ALAMDA=0, the covariance matrix
;	bvec	[OUTPUT] the BETA fitting-vector
;	svdthr	[INPUT; default=1e-6] threshold for singular values of
;		diagonal SVD matrix
;	trial	[OUTPUT] current trial values of the parameters
;	yfunc	[OUTPUT] best-fit Y(X;A)
;	tiptoe	[INPUT] if set, forces small steps on A
;	verbose	[INPUT] controls chatter
;	_extra	[INPUT] pass defined keywords to subroutines
;		LMCOEFF: FUNCS, SIG
;		ADJUSTIE: TIES
;		FUNCS: whatever is needed
;
;subroutines
;	LMCOEFF
;	SVDC
;	SVSOL
;	ADJUSTIE
;
;history
;	(C) Copr. 1986-92 Numerical Recipes Software =$j!]Y'1,).
;	translated from MRQMIN.F to IDL by Vinay Kashyap (Oct98)
;	changes: (i) matrix solution method from Gauss-Jordan elimination
;	  to Singular-Value Decomposition; (ii) if trial fails, ALAMDA
;	  increased by x{JUMPUP=3} instead of x10 and if trial succeeds,
;	  decreased by x{JUMPDN=0.1} to avoid obvious oscillations;
;	  (iii) call to COVSRT avoided; (iv) call to ADJUSTIE included;
;	  (v) no resetting chisq within procedure (VK; Oct98)
;	added keyword YFUNC (VK; Dec98)
;	allowed "inversion" of 1-element "matrices" bypassing
;	  SVDC and SVSOL (VK; Aug99)
;	now returns correct chisq if all params are frozen (VK; JanMM)
;	added keyword VERBOSE (VK; Nov04)
;-

;	initialize
ok='ok'
nx=n_elements(x) & ny=n_elements(y) & na=n_elements(a) & np=n_params()

;	usage
if nx eq 0 then ok='X not defined' else $
 if ny eq 0 then ok='Y(X) not defined' else $
  if na eq 0 then ok='model parameters not defined' else $
   if nx ne ny then ok='Y and X do not match' else $
    if np lt 4 then ok='insufficient parameters'
if ok ne 'ok' then begin
  print,'Usage: levmarq,x,y,a,chisq,funcs=funcs,sig=sig,ifit=ifit,alamda=alamda,$'
  print,'       jumpup=jumpup,jumpdn=jumpdn,ochisq=ochisq,alpha=alpha,covar=covar,bvec=bvec,$'
  print,'       svdthr=svdthr,trial=trial,yfunc=yfunc,/tiptoe,verbose=verbose'
  print,'  implements one iteration of Levenberg-Marquardt fitting method'
  message,ok,/info
  return
endif

;	keywords
vv=0L & if keyword_set(verbose) then vv=long(verbose[0])>1L
;
if n_elements(alamda) eq 0 then alamda=-1.
;
ns=n_elements(sig)
if ns ne ny then ysig=sqrt(abs(y)+0.75)+1. else ysig=sig
;
ni=n_elements(ifit) & kfit=intarr(na)+1
if ni ne na then message,$
  'thawing all parameters because IFIT has gone bad',/info else kfit=ifit
ofit=where(kfit ne 0,nfit)
if nfit eq 0 then begin
  message,'all parameters are frozen!',/info
  ;	just for the sake of chisq, go through LMCOEFF
  chisq=0.
  lmcoeff,x,y,a,ifit=kfit,chisq=chisq,yfunc=yfunc, _extra=e
  ochisq=chisq
  return
endif

;	singularities of SV Decomposition
if not keyword_set(svdthr) then svdthr=1e-6

;	initialize algorithm
if alamda lt 0 then begin
  alamda=0.001
  lmcoeff,x,y,a,ifit=kfit,alpha=alpha,bvec=bvec,chisq=chisq,yfunc=yfunc,$
	_extra=e
  ochisq=chisq
endif
atry=a
if not keyword_set(ochisq) then ochisq=chisq

if vv gt 1000 then stop,'HALTing point 1'

;	alter linearizing fitting matrix, by augmenting diagonal elements
covar=alpha		;save ALPHA in case new trial goes bad
da=bvec & jj=lindgen(nfit) & covar(jj,jj)=alpha(jj,jj)*(1.+alamda)

;	matrix solution
;the NumRec program uses Gauss-Jordan elimination to obtain the solution
;VK has changed it to SVD			;gaussj,covar,nfit,da,1,1

;	get SVD of matrix
mcvr=max(abs(covar))
if mcvr lt 1e-10 or mcvr gt 1e10 then double=1 else double=0
if n_elements(covar) gt 1 then begin
  oy=where(finite(covar) eq 0,moy)
  if moy gt 0 then covar(oy)=0.
  svdc,covar,w,u,v,double=double
endif else begin
  w=[covar(0)] & u=0*covar+1. & v=u
endelse

if vv gt 1000 then stop,'HALTing point 2'

;	catch singular w
ow=where(abs(w) lt svdthr,mow) & if mow gt 0 then w(ow)=0.

;	once converged, evaluate covariance matrix
;the NumRec progam called subroutine COVSRT here, but since VK has changed
;to using the SVD method, we'll get the covariance matrix directly from
;the SVD.. cf. SVDVAR, NumRec in C, 2nd ed., p679	;covsrt,covar,kfit
if alamda eq 0 then begin
  wti=0.*w & hw=where(abs(w) ge svdthr,mhw)
  if mhw gt 0 then wti(hw)=1./w(hw)^2			;weighting factors
  for i=0,nfit-1 do begin
    for j=0,i do begin
      sum=total(v(i,*)*v(j,*)*wti(*))
      covar(i,j)=sum & covar(j,i)=sum
    endfor
  endfor
  return
endif

if vv gt 1000 then stop,'HALTing point 3'

;	solve for step size
if n_elements(u) gt 1 then begin
  da=svsol(u,w,v,da)		;cf. NumRec in C, 2nd ed., eqn 15.5.14
endif else begin
  da=da/w			;U=[1], V=[1]
endelse

;	force small steps if TIPPYTOE is set
if keyword_set(tiptoe) then begin
  tippy=0.01	;1% change at a time by default
  if tiptoe(0) gt 0 and tiptoe(0) lt 1 then tippy=tiptoe(0) else $
   if tiptoe(0) gt 1 and tiptoe(0) lt 100 then tippy=tiptoe(0)/100. else $
    if tiptoe(0) gt 100 or tiptoe(0) lt 0 then tippy=abs(tiptoe(0))/a(ofit)
  delta_A=tippy*a(ofit)
  if vv gt 100 then print,da,delta_A
  for i=0,nfit-1 do if abs(da(i)) gt abs(delta_A(i)) then da(i)=(da(i)/abs(da(i)))*delta_A(i)
  if vv gt 100 then print,da,delta_A
endif

atry(ofit)=a(ofit)+da	;new parameters
adjustie,atry,_extra=e	;adjust as needed
trial=atry		;save for output

if vv gt 1000 then stop,'HALTing point 4'

;	test new parameters
lmcoeff,x,y,atry,ifit=kfit,alpha=covar,bvec=da,chisq=chisq,yfunc=yfunc,$
	_extra=e

if vv gt 100 then print,'LMCOEFF:',atry,da,w,chisq,ochisq

if vv gt 1000 then stop,'HALTing point 5'

;	did the trial succeed?
if chisq lt ochisq then begin		;(success, accept new solution
  alamda=0.1*alamda
  ;ochisq=chisq		;deliberately commented out -- done in FIT_LEVMAR
  alpha=covar & bvec=da & a=atry
endif else begin			;)(Failure, increase ALAMDA and return
  alamda=3.*alamda
  ;chisq=ochisq		;deliberately commented out -- handled in FIT_LEVMAR
endelse					;CHISQ v/s OCHISQ)

return
end
