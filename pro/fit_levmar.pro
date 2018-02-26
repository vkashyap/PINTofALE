pro fit_levmar,x,y,a,yfunc,freeze=freeze,$
	erra=erra,chisq=chisq,itmax=itmax,chithr=chithr,dumb=dumb,$
	tiptoe=tiptoe,verbose=verbose, _extra=e
;+
;procedure	fit_levmar
;	uses the Levenberg-Marquardt method to find best-fit parameters
;
;syntax
;	fit_levmar,x,y,a,yfunc,freeze=freeze,erra=erra,chisq=chisq,$
;	itmax=itmax,chithr=chithr,/dumb,ties=ties,vname=vname,$
;	jumpup=jumpup,jumpdn=jumpdn,svdthr=svdthr,funcs=funcs,sig=sig,$
;	/poisson
;
;parameters
;	x	[INPUT; required] data points
;	y	[INPUT; required] Y(X)
;		* sizes of X and Y must match
;	a	[I/O; required] parameters for user-supplied function
;		* on input, these are the initial guesses
;		* on output, these contain the best-fit values
;	yfunc	[OUTPUT] best-fit Y(X;A)
;
;keywords
;	freeze	[INPUT] freeze numbered parameters (starting from 0!)
;	erra	[OUTPUT] formal errors on the best-fit parameters
;	chisq	[OUTPUT] the chi-sq statistic denoting degree of agreement
;		between model and data
;	itmax	[INPUT; default=100] maximum number of iterations
;	chithr	[INPUT; default=0.1] stopping rule: ignore changes in CHISQ
;		smaller than this amount
;	dumb	[INPUT] if set, skips the part where the user can readjust
;		the parameter values in case the fit has gone bad.
;	tiptoe	[INPUT] if set, forces small steps on A
;	verbose	[INPUT] controls chatter
;	_extra	[INPUT] use this to pass defined variables to subroutines:-
;		ADJUSTIE: TIES, VNAME
;		LEVMARQ: JUMPUP, JUMPDN, SVDTHR
;		LMCOEFF: FUNCS, SIG
;		(note: FUNCS is actually the name of a _procedure_ --
;		CURVEFIT is to blame for the confusion)
;
;subroutines
;	ADJUSTIE
;	LEVMARQ
;	    LMCOEFF
;	    SVDC
;	    SVSOL
;
;history
;	vinay kashyap (Oct98)
;	added parameter YFUNC (VK; Dec98)
;	bug fix: crashing if freeze not set (VK; 99Aug)
;	now returns correct chisq if all params are frozen (VK; JanMM)
;	also adjust ERRA for ties (VK; FebMM)
;	what if adjusted ERRA becomes 0? (VK; MarMM)
;	added keywords TIPTOE and VERBOSE, changed behavior if
;	  dead end to now tiptoe around point (VK; Nov04)
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
  print,'Usage: fit_levmar,x,y,a,yfunc,freeze=freeze,erra=erra,chisq=chisq,$'
  print,'       itmax=itmax,chithr=chithr,/dumb,/tiptoe,ties=ties,vname=vname,$'
  print,'       jumpup=jumpup,jumpdn=jumpdn,svdthr=svdthr,funcs=funcs,sig=sig,$'
  print,'       /poisson'
  print,'  uses Levenberg-Marquardt method to obtain best-fit model parameters'
  message,ok,/info
  return
endif
erra=0.*a & chisq=1e20

;	keywords
tippy=0
if keyword_set(tiptoe) then tippy=tiptoe
;
vv=0L & if keyword_set(verbose) then vv=long(verbose[0])>1L
;
adjustie,a,_extra=e
;
nf=n_elements(freeze) & ifit=intarr(na)+1	;default is to fit all
if nf ge na then begin
  message,'all parameters frozen',/info
  ;	get CHISQ and YFUNC and quit
  levmarq,x,y,a,chisq,ifit=intarr(na),alamda=-1,ochisq=ochisq,yfunc=yfunc,$
	tiptoe=tippy,verbose=vv, _extra=e
  return
endif else begin
  if nf gt 0 then begin
    if freeze(0) ne -1 then ifit(long(freeze))=0
  endif
endelse
;
if not keyword_set(itmax) then itmax=100
;
if not keyword_set(chithr) then chithr=0.1

;	go forth and fit
iter=0 & dchi=2*chithr & alamda=-1.
while dchi gt chithr do begin 		;{iterate fit parameters

  iter=iter+1
  levmarq,x,y,a,chisq,ifit=ifit,alamda=alamda,ochisq=ochisq,$
	alpha=alpha,covar=covar,bvec=bvec,trial=trial,yfunc=yfunc,$
	tiptoe=tippy,verbose=vv, _extra=e

  if vv gt 0 then print,'>>>',iter,a,chisq,ochisq
  if vv gt 1000 then stop,'HALTing, type .CON to continue'
  if ochisq gt chisq then begin		;(if trial succeeded..
    dchi=ochisq-chisq
    ochisq=chisq
  endif					;CHISQ < OCHISQ)

  ;(what if fit keeps getting worse..?
  mxalp=max(abs(alpha),min=mnalp)
  oalp=where(alpha ne 0,moalp)
  if moalp gt 0 then mxalp=max(abs(alpha(oalp)),min=mnalp)
  ;mxalp=max(abs(alpha(where(alpha ne 0))),min=mnalp)
  if alamda gt 1e10*mxalp or alamda lt 1e-15*mnalp then begin
    message,'dead end.. bailing out',/info
    if vv gt 10 then print,alamda
    ;message,$
    ;  'WARNING: guess was hopelessly bad, and we are going nowhere fast.',/info
    if not keyword_set(dumb) then begin
      print,'current parameter values are:'
      for i=0,na-1 do begin
        cc='A('+strtrim(i,2)+') = '+strtrim(a(i),2)
        if ifit(i) eq 0 then cc=cc+' frozen ' else cc=cc+'   free '
        print,cc
      endfor
      print,'type a to adjust A, z to halt momentarily, any key to quit'
      c='' & c=strlowcase(get_kbrd(1))
    endif else c=''
    ;help,1e10*mxalp,alamda,1e-15*mnalp,dchi,chisq,iter
    ;stop
    case c of
      'a': begin
	k=0
	print,'type parameter indices and corresponding value, end with <CR> or q'
	while k ge 0 and k le na-1 do begin
	  ck='' & read,ck & ck=strlowcase(strtrim(ck,2))
	  k=-1
	  if ck ne '' and ck ne 'q' then begin
	    isp=strpos(cc,' ',1) & icm=strpos(cc,',',1)
	    if isp ge 0 or icm ge 0 then reads,ck,k,val
	  endif
	  if k ge 0 and k le na-1 then a(k)=val
	endwhile
	adjustie,a,_extra=e
        iter=1 & old_alamda=alamda & alamda=-1
      end
      'z': begin
        iter=1 & old_alamda=alamda & alamda=-1
        stop,'HALTING.. type .CON to continue'
      end
      else: begin
	;iter=itmax
	;if n_elements(tiptoe) eq 0 then tippy=1
	if not keyword_set(tiptoe) then iter=itmax
	if finite(alamda) eq 0 then alamda=1e10*mxalp
      end
    endcase
  endif					;ALAMDA has unreasonable value)

  ;	stopping rules
  if dchi lt chithr and chisq gt 2*(nx-nf) then dchi=2*chithr	;keep going
  if iter ge itmax then dchi=0.
  if chisq lt (1.-chithr)*(nx-nf) then dchi=0.
  ;adjustie,a,_extra=e
endwhile					;d(CHI)>CHITHR}

;	but..
if iter ge itmax and ochisq lt chisq then begin
  message,'WARNING: too many iterations -- this may not be a good fit',/info
endif

;	and one final call
levmarq,x,y,a,chisq,ifit=ifit,alamda=0,ochisq=ochisq,$
	alpha=alpha,covar=covar,bvec=bvec,yfunc=yfunc,$
	verbose=vv, _extra=e

;	get the formal errors
;the diagonal elements of the matrix COVAR are the standard errors
;(see NumRec in C, 2nd ed., p677).  COVAR is a square matrix of size
;NFITxNFIT, so be careful when putting it into ERRA
ofit=where(ifit gt 0,nfit) & diag=lindgen(nfit) * (nfit+1)
erra(ofit)=sqrt(covar(diag))

erra = (erra > 1e-6) < 1e12
o=where(finite(erra) eq 0,mo)
if mo gt 0 then erra(o)=0.1

;	adjust the errors also to match the ties
olderra=erra & tmpu=a+erra & tmpl=a-erra
adjustie,tmpu,_extra=e & adjustie,tmpl,_extra=e
erra=0.5*(tmpu-tmpl)
for i=0,na-1 do if ifit(i) ne 0 and erra(i) eq 0 then erra(i)=olderra(i)

return
end
