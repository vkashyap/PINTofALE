pro erors,x,y,a,erru,errl,ysig=ysig,freeze=freeze,dchi=dchi,algo=algo,$
	maxstep=maxstep,verbose=verbose,erra=erra,yfunc=yfunc, _extra=e
;+
;procedure	erors
;	compute asymmetric error bars on fit parameters by projecting the
;	chi-sq surface onto each parameter axis.  what this means is that
;	for each non-frozen parameter, step through a range of values of
;	the parameter, finding the best-fit chi-square calculated by fitting
;	the rest of the parameters, and report that range where the chi-square
;	increases by a set amount from the minimum.
;
;	this program tries to be a wee bit clever by starting from the
;	best-fit value and then stumbling about on either side until
;	the chi-sq goes above the mark; the length of the strides change
;	according to the projected location of said mark.
;
;syntax
;	erors,x,y,a,erru,errl,ysig=ysig,freeze=freeze,dchi=dchi,$
;	  dchi=dchi,algo=algo,maxstep=maxstep,verbose=verbose,$
;	  erra=erra,yfunc=yfunc,$
;	  itmax=itmax,chithr=chithr,/dumb, jumpup=jumpup,jumpdn=jumpdn,$
;	  svdthr=svdthr,funcs=funcs, function_name=function_name,$
;	  type=type, missing=missing,/fwhm,/norm,betap=betap, /poisson
;
;parameters
;	x	[INPUT; required] data points
;	y	[INPUT; required] Y(X)
;		* sizes of X and Y must match
;	a	[I/O; required] parameters for user-supplied function
;		* on input, these are assumed to be initial guesses
;		* on output, these contain the best-fit values
;	erru	[OUTPUT; required] upper limits of confidence range interval
;	errl	[OUTPUT] lower limits of confidence range interval
;		* if ERRL is not specified on input, ERRU will contain the
;		  average value of the upper and lower >>deviations<<.
;
;keywords
;	ysig	[INPUT] standard deviations on Y
;		* default=sqrt(abs(Y)+0.75)+1
;		* if single element, then sig(Y[*])=YSIG(0)
;		* if -ve, taken to be the fractional error
;	freeze	[INPUT] freeze numbered parameters (index starts from 0!)
;	dchi	[INPUT] how big a change in chi-square to look for?
;		* default is 2.7 (corresponding to 90% CL)
;	algo	[INPUT] fitting algorithm
;		* only the following are implemented:
;		-- LevMarq+SVD (default; calls FIT_LEVMAR)
;		-- IDL-Curvefit (calls CURVE_FIT)
;	maxstep	[INPUT] maximum number of steps to take before giving up
;		on actually finding the bounds
;		* default is 100
;	verbose	[INPUT] verbosity level
;	erra	[OUTPUT] formal "curvature" errors on the best-fit parameters
;	yfunc	[OUTPUT] best-fit Y(X;A)
;	_extra	[INPUT] pass defined variables to subroutines:-
;		FIT_LEVMAR: ITMAX, CHITHR, DUMB
;		ADJUSTIE: TIES, VNAME
;		LEVMARQ: JUMPUP, JUMPDN, SVDTHR
;		LMCOEFF: FUNCS, POISSON
;		CURVE_FIT: FUNCTION_NAME
;		note:-	FUNCS and FUNCTION_NAME refer to name of
;			user-defined function that takes as input X
;			and A, and returns YMODEL(X;A), and the
;			partial derivatives of the model function
;			wrt A.  Any function that was written to work
;			with CURVEFIT or GHRS' WFIT will do.
;			The default for FIT_LEVMAR is X3MODEL.
;		MK_3MODEL: TYPE
;		MK_GAUSS: MISSING, FWHM, NORM
;		MK_LORENTZ: BETAP, MISSING, NORM
;
;history
;	vinay kashyap (MM.I) (yes, I _do_ know how to spell "error")
;	added call to ADJUSTIE to handle constraints on ERRU,ERRL (VK; FebMM)
;	what if adjusted ERRA becomes 0? (VK; MarMM)
;	allowed halting with either "q" or "x" (VK; SepMM)
;	force extra iteration after "r"; increased wait time (VK; JanMMI)
;	added confirmation check for too many thawed params (VK; FebMMI)
;	interpol was crashing because it was getting only 1 element
;	  arrays; now correctly updates all parameters if better fit is found
;	  (VK; Aug01)
;	made various bug fixes that was causing program to go bonkers for
;	  some special cases, such as small fluxes (VK; Apr02)
;	added hooks into MPFIT (Liwei Lin/VK; Oct02)
;	bug correction re VERBOSE (LL; Apr03)
;	bug correction, first step was failing sometimes when input A were
;	  integers (VK; Mar08)
;-

forward_function mpfitfun

;	usage
ok='ok' & np=n_params() & nx=n_elements(x) & ny=n_elements(y) & na=n_elements(a)
if np lt 4 then ok='insufficient parameters' else $
 if nx eq 0 then ok='Input X undefined' else $
  if ny eq 0 then ok='Input Y undefined' else $
   if nx eq 1 then ok='Input not an array?' else $
    if nx ne ny then ok='X and Y(X) not compatible' else $
     if na eq 0 then ok='model parameters undefined'
if ok ne 'ok' then begin
  print,'Usage: erors,x,y,a,erru,errl,ysig=ysig,freeze=freeze,dchi=dchi,$'
  print,'       dchi=dchi,algo=algo,verbose=verbose,erra=erra,yfunc=yfunc,$'
  print,'       itmax=itmax,chithr=chithr,/dumb,$'
  print,'       jumpup=jumpup,jumpdn=jumpdn,svdthr=svdthr,$'
  print,'       funcs=funcs, function_name=function_name,$'
  print,'       type=type, missing=missing,/fwhm,/norm,betap=betap,/poisson'
  print,'  compute error bars on fit parameters'
  if np ne 0 then message,ok,/info
  return
endif

;	inputs
;	VERBOSE
vv=0 & if keyword_set(verbose) then vv=long(verbose(0)) > 1
;	YSIG
sigy=sqrt(abs(y)+0.75)+1. & nys=n_elements(ysig)
if nys eq ny then sigy=ysig else begin		;default is Poisson
  if nys eq 1 then sigy(*)=ysig(0)		;constant
endelse
om=where(sigy lt 0,mom)
if mom gt 0 then sigy(om)=y(om)*abs(sigy(om))	;fractional
ywt=sigy & ow=where(ywt gt 0,mow) & if mow gt 0 then ywt(ow)=1./ywt(ow)
;	FREEZE
ifit=intarr(na)+1 & nf=n_elements(freeze)	;default is fit all
if nf gt 0 then begin
  if freeze(0) eq -1 then nf=0 else begin	;ignore keyword
    if nf gt na then ifit=intarr(na) else ifit(long(freeze))=0	;as aprop
  endelse
endif
;	DCHI
delchi=2.7 & if keyword_set(dchi) then delchi=float(dchi(0))
if delchi lt 0 then begin
  delchi=-delchi	;placeholder for grander designs
endif
;	ALGO
algo_type='LevMarq+SVD'			;default
if keyword_set(algo) then fitalgo=strtrim(algo(0),2) else fitalgo=algo_type
i0=-1L
i00=strpos(strlowcase(fitalgo),'lev',0) & i0 = i0 > i00
i01=strpos(strlowcase(fitalgo),'mar',0) & i0 = i0 > i01
i02=strpos(strlowcase(fitalgo),'lm',0) & i0 = i0 > i02
i03=strpos(strlowcase(fitalgo),'svd',0) & i0 = i0 > i03
i10=strpos(strlowcase(fitalgo),'idl',0) & i0 = i0 > i10
i11=strpos(strlowcase(fitalgo),'curve',0) & i0 = i0 > i11
i20=strpos(strlowcase(fitalgo),'mp',0) & i0 = i0 > i20
i21=strpos(strlowcase(fitalgo),'craig',0) & i0 = i0 > i21
if i00 ge 0 or i01 ge 0 or i02 ge 0 or i03 ge 0 then algo_type='LevMarq+SVD'
if i10 ge 0 or i11 ge 0 then algo_type='IDL-Curvefit'
if i20 ge 0 or i21 ge 0 then algo_type='MPFIT'
if algo_type eq 'IDL-Curvefit' and nf ne 0 then begin
  message,algo_type+' does not understand parameter freezing',/info
endif
if i0 lt 0 then message,$
	fitalgo+': not understood.  reverting to '+algo_type,/info
if vv ge 1 then message,'using fitting algorithm: '+algo_type,/info
;	MAXSTEP
stepmax=100L & if keyword_set(maxstep) then stepmax=long(maxstep(0)) > 0
stepmaxmax=stepmax*100L > 100000L

;	get the best-fit
if algo_type eq 'LevMarq+SVD' then fit_levmar,x,y,$
	a,yfunc,freeze=freeze,erra=erra,chisq=chisq,sig=sigy, _extra=e
if algo_type eq 'IDL-Curvefit' then yfunc=curve_fit(x,y,$
	ywt,a,erra,chisq=chisq, _extra=e)
if algo_type eq 'MPFIT' then begin
  ;freeze is where thaw=0
  PARINFO=replicate({fixed:0,limited:[0,0],limits:[0,0]},n_elements(aa))
  ;freeze it
  for j=0,nf-1 do parinfo(freeze(j)).fixed(0)=1
  userfuncs=funcs
  functarg=create_struct(e,'normflx',znorm,'type',type)
  if strpos(userfuncs,'_f',0) lt 0 then userfuncs=userfuncs+'_f'
  results=MPFITFUN(userfuncs,x,y,sigy,aa,PERROR=erra,yfit=yfit,$
	parinfo=parinfo,functarg=functarg, _extra=e)
  ;chisq=results(n_elements(results)-1L)
  ;aa=results[0:n_elements(results)-2L]
  aa=results
  chisq=total( ((yfit-zz)^2)/((sigy)^2) )
endif

;	outputs
erru=0.*erra & errl=erru

;	how many thawed parameters?
ot=where(ifit gt 0,mot)		;which parameters are thawed
if vv gt 1 then message,'There are '+strtrim(mot,2)+' thawed parameters',/info
if mot gt 15 then begin
  message,'This will eons to complete and the answer will still be 42',$
	/info
  message,'type "x" NOW to return without further damage',/info
  wait,3 & c=get_kbrd(0)
  if strlowcase(c) eq 'x' then begin
    if vv gt 5 then message,'quitting',/info
    return
  endif else message,'.. proceeding to find errors for all '+$
	strtrim(mot,2)+' thawed parameters',/info
endif

;	project the chi-sq surface
for i=0L,mot-1L do begin		;{find error for each parameter

  k=ot(i)
  aval=[1.0*a(k)] & chival=[chisq]
  kfit=ifit & kfit(k)=0 & kfreeze=where(kfit eq 0)

  starthere:		;((return here only if new minimum
  stretch=0.5		;start off by checking every point at half-sigma
  astep=stretch*erra(k)

  go_on=1 & kstep=0L & zstep=0L
  while go_on do begin			;{hunt for the bounds
    kstep=kstep+1L & zstep=zstep+1L
    while (aval(0)+astep)-aval(0) eq 0 do begin
      stretch=stretch*sqrt(2.) & astep=stretch*erra(k)
      if vv ge 2 then print,'stretch,astep:',stretch,astep
    endwhile
    if vv le 1 then kilroy else print,$
	'par,step,a-,a,a+ :',$
	strtrim(k,2)+' '+strtrim(kstep,2),$
	a(k)-kstep*astep,a(k),a(k)+kstep*astep

    ;	the forward step
    aa=1.0*a & aa(k)=a(k)+kstep*astep
    ok=where(aval ne aa(k),mok) & if mok eq 0 then message,'BUG!'
    aval=aval(ok) & chival=chival(ok)
    if algo_type eq 'LevMarq+SVD' then fit_levmar,x,y,$
	aa,tmp,freeze=kfreeze,erra=ea,chisq=chiv,sig=sigy, _extra=e
    if algo_type eq 'IDL-Curvefit' then tmp=curve_fit(x,y,$
	ywt,aa,ea,chisq=chiv, _extra=e)
    if algo_type eq 'MPFIT' then begin
      ;freeze is where thaw=0
      ;PARINFO=replicate({fixed:0,limited:[0,0],limits:[0,0]},n_elements(aa))
      ;freeze it
      ;for j=0,nf-1 do parinfo(freeze(j)).fixed(0)=1
      userfuncs=funcs
      functarg=create_struct(e,'normflx',znorm,'type',type)
      if strpos(userfuncs,'_f',0) lt 0 then userfuncs=userfuncs+'_f'
      results=MPFITFUN(userfuncs,x,y,sigy,aa,PERROR=erra,yfit=yfit,$
	parinfo=parinfo,functarg=functarg, _extra=e)
      ;chisq=results(n_elements(results)-1L)
      ;aa=results[0:n_elements(results)-2L]
      aa=results
      chisq=total( ((yfit-zz)^2)/((sigy)^2) )
    endif
	;check to eliminate NaNs and Infs from parabola (VK; Nov04)
    ;if finite(chiv) ne 0 then begin
      aval=[aval,aa(k)] & chival=[chival,chiv]
    ;endif else begin
      ;if vv ge 5 then message,'Problem with NaNs',/informational
      ;chiv=max(chival,/nan)
    ;endelse
    if chiv lt chisq then begin
      if vv ge 1 then message,'found a smaller minimum @ +'+$
	strtrim(kstep*astep,2)+' w. chisq '+$
	strtrim(chiv,2)+', cf. '+strtrim(chisq,2),/info
      ;a(k)=aa(k)
      if vv gt 5 then begin
	print,'Old params: ',a
	print,'New params: ',aa
	if vv ge 10 then begin
	  plot,x,y,/xs,psym=10 & oplot,x,tmp,thick=2,color=150
	endif
      endif
      a=aa & yfunc=tmp & chisq=chiv & kstep=0L
      adjustie,a,_extra=e	;constraint corrections, if needed
      aval=[aa(k)] & chival=[chiv]
      ok=where(ea gt 0,mok) & if mok gt 0 then erra(ok)=ea(ok)
      goto,starthere			;back up)
    endif

    ;	the backward step
    aa=1.0*a & aa(k)=a(k)-kstep*astep
    ok=where(aval ne aa(k),mok) & if mok eq 0 then message,'BUG!'
    aval=aval(ok) & chival=chival(ok)
    if algo_type eq 'LevMarq+SVD' then fit_levmar,x,y,$
	aa,tmp,freeze=kfreeze,erra=ea,chisq=chiv,sig=sigy, _extra=e
    if algo_type eq 'IDL-Curvefit' then tmp=curve_fit(x,y,$
	ywt,aa,ea,chisq=chiv, _extra=e)
    if algo_type eq 'MPFIT' then begin
      ;freeze is where thaw=0
      PARINFO=replicate({fixed:0,limited:[0,0],limits:[0,0]},n_elements(aa))
      ;freeze it
      for j=0,nf-1 do parinfo(freeze(j)).fixed(0)=1
      userfuncs=funcs
      functarg=create_struct(e,'normflx',znorm,'type',type)
      if strpos(userfuncs,'_f',0) lt 0 then userfuncs=userfuncs+'_f'
      results=MPFITFUN(userfuncs,x,y,sigy,aa,PERROR=erra,yfit=yfit,$
	parinfo=parinfo,functarg=functarg, _extra=e)
      ;chisq=results(n_elements(results)-1L)
      ;aa=results[0:n_elements(results)-2L]
      aa=results
      chisq=total( ((yfit-zz)^2)/((sigy)^2) )
    endif
	;check to eliminate NaNs and Infs from parabola (VK; Nov04)
    ;if finite(chiv) ne 0 then begin
      aval=[aval,aa(k)] & chival=[chival,chiv]
    ;endif else begin
      ;if vv ge 5 then message,'Problem with NaNs',/informational
      ;chiv=max(chival,/nan)
    ;endelse
    if chiv lt chisq then begin
      if vv ge 1 then message,'found a smaller minimum @ -'+$
	strtrim(kstep*astep,2)+' w. chisq '+$
	strtrim(chiv,2)+', cf. '+strtrim(chisq,2),/info
      ;a(k)=aa(k)
      if vv gt 5 then begin
	print,'Old params: ',a
	print,'New params: ',aa
	if vv ge 10 then begin
	  plot,x,y,/xs,psym=10 & oplot,x,tmp,thick=2,color=150
	endif
      endif
      a=aa & yfunc=tmp & chisq=chiv & kstep=0L
      adjustie,a,_extra=e	;constraint corrections, if needed
      aval=[aa(k)] & chival=[chiv]
      ok=where(ea gt 0,mok) & if mok gt 0 then erra(ok)=ea(ok)
      goto,starthere			;back up)
    endif

    ;	find the range
    oe=sort(aval) & xe=aval(oe) & ye=chival(oe)
    chimin=min(ye,i0) & i1=n_elements(xe)
    if xe(i0) ne a(k) then begin
      ;this appears when a parameter is tied to another but
      ;is not frozen as it ought to be
      message,'Parameter #'+strtrim(k,2)+' ought to have been frozen!',/info
    endif

    ;	(this bit to make sure we're not bumping up onto a lee shore
    leelo=0 & leehi=0
    if xe(0) eq xe(1) then leelo=1 & if xe(i1-2L) eq xe(i1-1L) then leehi=1
    if ye(0) gt chimin+delchi then leelo=0
    if ye(i1-1L) gt chimin+delchi then leehi=0
    ;xep=a(k)+(1-leehi)*astep & xem=a(k)-(1-leelo)*astep
    xep=xe(i1-1L) & xem=xe(0)
    oi=lindgen(i1)
    if leehi eq 1 then oi(i1-2L)=-1
    ;if leelo eq 1 then oi(i1-1L)=-1	;eh?  this doesn't make sense.. (VK/Apr02)
    if leelo eq 1 then oi(1L)=-1
    oi=where(oi ge 0)	;these are basically the "legit" points so far
    ;	the above bit to ensure we're within proper bounds)
    xxe=xe(i0:*) & yye=ye(i0:*) & oe=uniq(xxe) & xxe=xxe(oe) & yye=yye(oe)
    ;if leehi eq 0 then xep=(interpol(xe(i0:*),ye(i0:*),chisq+delchi))(0)
    ;if leehi eq 0 and n_elements(xep) gt 1 then xep=(interpol(xxe,yye,chisq+delchi))(0)
    if leehi eq 0 and n_elements(xxe) gt 1 then xep=(interpol(xxe,yye,chisq+delchi))(0)
    xxe=xe(0:i0) & yye=ye(0:i0) & oe=uniq(xxe) & xxe=xxe(oe) & yye=yye(oe)
    ;if leelo eq 0 then xem=(interpol(xe(0:i0),ye(0:i0),chisq+delchi))(0)
    ;if leelo eq 0 and n_elements(xem) gt 1 then xem=(interpol(xxe,yye,chisq+delchi))(0)
    if leelo eq 0 and n_elements(xxe) gt 1 then xem=(interpol(xxe,yye,chisq+delchi))(0)

    ;	reset the step size if necessary
    xstep=0.5*((xep-aa(k))+(aa(k)-xem))
    if xstep le 0 then begin
      ;	this happens if (a) the chi-sq surface is non-parabolic and/or
      ;	ratty, or (b) we have only found a local minimum
      astep=astep*1.1
      ;	don't end at this step
      if stepmax-kstep lt 3 then kstep=kstep-1L
      message,'non-parabolic CHISQ',/informational
      if vv ge 2 then print,'New step size:',astep
    endif
    if xstep/astep gt 10 then begin
      astep=astep*2.
      ;	throw away useless points to avoid hitting stopping rule
      aval=aval(oi) & chival=chival(oi)
      kstep=kstep-1L & leehi=0 & leelo=0	;don't end at this step
      if vv ge 2 then print,'New step size:',astep
    endif else begin
      if abs(xstep/astep) lt 1. then begin
        astep=astep/2.5	;not 2.0 because we don't want to repeat points
	if xstep gt 0 then begin
          ;	throw away useless points to avoid hitting stopping rule
          aval=aval(oi) & chival=chival(oi)
          kstep=kstep-1L & leehi=0 & leelo=0	;don't end at this step
	endif
        if vv ge 2 then print,'New step size:',astep
      endif
    endelse

    ;	stopping rules
    yep=max(ye(i0:*)) & yem=max(ye(0:i0)) & chimax=chimin+1.5*delchi
    ok='ok'
    if yep gt chimax and yem gt chimax then ok='done'
    if leelo eq 1 and yep gt chimax then ok=$
	'found upper bound; lower bound beyond reach'
    if leehi eq 1 and yem gt chimax then ok=$
	'found lower bound; upper bound beyond reach'
    if leelo eq 1 and leehi eq 1 then ok=$
	'hit the maximum in allowed range'
    if kstep gt stepmax then ok='too many iterations; giving up'
    if zstep gt stepmaxmax then ok='BEYOND HELP!  Stuck in infinte loop?'
    if ok ne 'ok' and ok ne 'done' then message,ok,/info
    if ok ne 'ok' then go_on=0

    ;	populate output
    if go_on eq 0 then begin
      erru(k)=xep & errl(k)=xem
      if vv ge 4 then begin
	tt=strtrim(k,2)+': '+strtrim(xem,2)+' < '+strtrim(a(k),2)+$
		' < '+strtrim(xep,2)
	plot,xe,ye,yr=chimin+[-1,1.5*delchi],title=tt,ytitle='!4v!3!u2!n'
	oplot,a(k)*[1,1],chimin+[-2,3*delchi],line=1
	oplot,xe,0*ye+chimin+delchi,line=1
	if vv ge 10 then print,'type q or x to halt, r to repeat'
	if vv gt 5 then wait,vv/2.
	if vv ge 10 then c=get_kbrd(0) else c=''
	if strlowcase(c) eq 'x' or strlowcase(c) eq 'q' then $
		stop,'Halting.. type .CON to continue'
	if strlowcase(c) eq 'r' then go_on=1
      endif
    endif

  endwhile				;GO_ON}

endfor					;I=0,MOT-1}

;	adjust the errors to suit the constraints
olderru=erru & olderrl=errl
adjustie,erru,_extra=e & adjustie,errl,_extra=e
for i=0,n_elements(ifit)-1 do begin
  if ifit(i) ne 0 and erru(i) eq 0 then erru(i)=olderru(i)
  if ifit(i) ne 0 and errl(i) eq 0 then errl(i)=olderrl(i)
endfor

;	convert to symmetric error bars if necessary
if np eq 4 then begin
  tmpu=erru-a & tmpl=a-errl & erru=0.5*(abs(tmpu)+abs(tmpl))
endif

if vv ge 1 then print,''
if vv ge 1 then print,'Done finding projected errors'

return
end
