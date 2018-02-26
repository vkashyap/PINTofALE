pro lmcoeff,x,y,a,funcs=funcs,sig=sig,ifit=ifit,alpha=alpha,bvec=bvec,$
	chisq=chisq,yfunc=yfunc,poisson=poisson, _extra=e
;+
;procedure	lmcoeff
;	used by LEVMARQ (formerly MRQMIN.F) to evaluate the linearized
;	fitting matrix ALPHA and vector BVEC (BETA in Numerical Recipes)
;	and calculate chi^2
;
;syntax
;	lmcoeff,x,y,a,funcs=funcs,sig=sig,ifit=ifit,alpha=alpha,bvec=bvec,$
;	chisq=chisq,yfunc=yfunc,poisson=poisson
;
;parameters
;	x	[INPUT; required] data points
;	y	[INPUT; required] Y(X)
;		* sizes of X and Y must match
;	a	[INPUT; required] parameters for user-supplied procedure
;
;keywords
;	funcs	[INPUT] name of user-defined procedure that takes as input X
;		and A, and returns YMODEL(X;A), and the partial derivatives
;		of the model function wrt A.  Any procedure that was written
;		to work with CURVEFIT or GHRS' WFIT will do.
;		* default is set to X3MODEL
;		* why is it a procedure and not a function?
;		  well, ask the person who wrote CURVEFIT.
;	sig	[INPUT; default=sqrt(abs(Y)+0.75)+1] standard deviations on Y
;	ifit	[INPUT] integer array of same size as A, with 0's indicating
;		frozen parameters and 1's indicating thawed parameters
;		* if IFIT does not match A for any reason, all parameters
;		are assumed to be thawed
;	alpha	[OUTPUT] fitting matrix
;	bvec	[OUTPUT] fitting vector
;	chisq	[OUTPUT] chi^2, or ln(p(D|M)) if POISSON is set
;	yfunc	[OUTPUT] Y(X;A) -- function values
;	poisson	[INPUT] if set, returns the log of the poisson likelihood,
;		p(D|M) instead of the chi-sq
;		* do _not_ use this keyword for background-subtracted data!!
;	_extra	[JUNK] here only to prevent crashing the program
;
;subroutines
;	LNPOISSON
;
;history
;	(C) Copr. 1986-92 Numerical Recipes Software =$j!]Y'1,).
;	translated from MRQCOF.F to IDL by Vinay Kashyap (Oct98)
;	now returns correct chisq if all params are frozen (VK; JanMM)
;	added keyword POISSON and call to LNPOISSON (VK; Aug01)
;-

;	initialize
ok='ok'
nx=n_elements(x) & ny=n_elements(y) & na=n_elements(a)

;	usage
if nx eq 0 then ok='X not defined' else $
 if ny eq 0 then ok='Y(X) not defined' else $
  if na eq 0 then ok='model parameters not defined' else $
   if nx ne ny then ok='Y and X do not match'
if ok ne 'ok' then begin
  print,'Usage: lmcoeff,x,y,a,funcs=funcs,sig=sig,ifit=ifit,$'
  print,'       alpha=alpha,bvec=bvec,chisq=chisq,yfunc=yfunc'
  print,'  subroutine of LEVMARQ to evaluate linearized fitting matrix'
  return
endif

;	keywords
;
if not keyword_set(funcs) then funcs='x3model'
;
ns=n_elements(sig)
if ns ne ny then ysig=sqrt(abs(y)+0.75)+1. else ysig=sig
sig2i=0.*ysig+1.
oo=where(ysig gt 0,moo) & if moo gt 0 then sig2i(oo)=1./ysig(oo)^2
;
ni=n_elements(ifit) & kfit=intarr(na)+1
if ni ne na then message,$
  'thawing all parameters because IFIT has gone bad',/info else kfit=ifit

;	initialize
ofit=where(kfit ne 0,nfit)
if nfit eq 0 then begin
  message,'all parameters are frozen',/info
  ;	get chisq and return
  chisq=0.
  call_procedure,funcs,x,a,yfunc,dyda, _extra=e
  if keyword_set(poisson) then begin
    tmp=fltarr(ny)
    for i=0L,ny-1L do tmp[i]=lnpoisson(y[i],yfunc[i],/chilike)
    chisq=total(tmp)
  endif else begin
    dy=y-yfunc & chisq=total(dy^2*sig2i)
  endelse
  return
endif
alpha=fltarr(nfit,nfit) & bvec=fltarr(nfit) & chisq=0.

;	call the external function
call_procedure,funcs,x,a,yfunc,dyda, _extra=e
dy=y-yfunc

;	generate the output matrices by summing over each data point
for j=0,nfit-1 do begin
  jj=ofit(j)
  wt=dyda(*,jj)*sig2i
  for k=0,j do alpha(j,k)=total(wt*dyda(*,ofit(k)))	;fill in ALPHA
  bvec(j)=total(dy*wt)					;fill in BVEC
  if finite(bvec(j)) eq 0 then bvec(j)=0.
endfor

;	fill in symmetric side of ALPHA
for j=1,nfit-1 do for k=0,j-1 do alpha(k,j)=alpha(j,k)

;	get CHISQ
if keyword_set(poisson) then begin
  tmp=fltarr(ny)
  for i=0L,ny-1L do tmp[i]=lnpoisson(y[i],yfunc[i],/chilike)
  chisq=total(tmp)
endif else chisq=total(dy^2*sig2i)

return
end
