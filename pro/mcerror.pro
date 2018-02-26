pro mcerror,x,y,a,erru,errl,ysig=ysig,freeze=freeze,algo=algo,$
	verbose=verbose,erra=erra,yfunc=yfunc,mpfunc=mpfunc,$
        mciter=mciter,mcconf=mcconf,mcsims=mcsims, _extra=e
;+
;procedure	mcerror
;	procedure that estimates confidence bounds for model fits by
;	fitting to different monte-carlo realizations of the data and
;	estimating the bounds from the distribution of the resulting
;	parameter values.
;
;	generates NN number of synthetic data sets by randomly varying
;	each data point within its estimated error. fits these NN
;	synthetic data  sets to generate parameter sets a(1),a(2),..a(NN),
;	so that we can estimate distribution of a(x)-a(best).  Uses this
;	as a basis for confidence bound estimation.  Finds smallest possible
;	bounds giving required confidence (as opposed to 'counting' up
;	and down from best fit by required fraction)
;            
;syntax
;	mcerror,x,y,a,erru,errl,ysig=ysig,freeze=freeze,mcfuncs=mpfuncs,$
;	algo=algo,verbose=verbose, erra=erra,yfunc=yfunc,mciter=mciter,$
;	itmax,chitr=chitr,mcconf=mcconf,/dumb, jumpup=jumpup,jumpdn=jumpdn,$
;	=svdthr,funcs=funcs, function_name=function_name,$
;	ype, missing=missing,/fwhm,/norm,betap=betap, /poisson
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
;       mciter  [INPUT] number of additional iterations for each parameter:
;               total iterations = (mciter x number of parameters) > 2xmciter+1                               
;               default: mciter=30
;       mcconf  [INPUT] (two-sided) confidence for bound estimation,
;                       default: mpfunc =0.68
;       mcsims  [OUTPUT] 2 dimensional array MCSIMS==MCSIMS(free,sims)
;               where free are the number of thawed paramters and sims 
;               is the number of simulations. contains all simulations 
;               on output 
;       mpfunc  [INPUT] contains name of input function to fit
;                       REQUIRED WHEN USING MPFIT
;	freeze	[INPUT] freeze numbered parameters (index starts from 0!)
;	algo	[INPUT] fitting algorithm
;		* only the following are implemented:
;		-- LevMarq+SVD 
;		-- IDL-Curvefit 
;               -- Mardkwardt's MPFIT
;	verbose	[INPUT] verbosity level
;	erra	[OUTPUT] formal "curvature" errors on the best-fit parameters
;	yfunc	[OUTPUT] best-fit Y(X;A)
;	_extra	[INPUT] pass defined variables to subroutines:-
;		FIT_LEVMAR: ITMAX, CHITHR, DUMB
;		ADJUSTIE: TIES, VNAME
;		LEVMARQ: JUMPUP, JUMPDN, SVDTHR
;		LMCOEFF: FUNCS, POISSON
;		CURVE_FIT: FUNCTION_NAME
;		note:-	FUNCS, FUNCTION_NAME, and MPFUNC refer to name
;			of user-defined function that takes as input X
;			and A, and returns YMODEL(X;A), and the partial
;			derivatives of the model function wrt A.  Any
;			function that was written to work with CURVEFIT
;			or GHRS' WFIT will do.
;			The default for FIT_LEVMAR is X3MODEL.
;		MK_3MODEL: TYPE
;		MK_GAUSS: MISSING, FWHM, NORM
;		MK_LORENTZ: BETAP, MISSING, NORM
;
;history
;	Liwei Lin (Apr03) idenentical to ERORS.pro in input/output/parameter checks
;	  added keyword MPFUNC for compatability with MPFIT            
;	  added keyword MCCONF so user can toggle confidence required
;                 BUG FIX, mcerror not mcerors LL/SC JUL03
;                 BUG FIX, yfunc was being recursivly overwritten in
;                          synthetic loop LL/SC JUL03
;	  added forward_function (VK; Jul03)
;	some cosmetic changes (VK; Feb04)
;         added MCSIMS (LL; Mar04) 
;       BUG FIX algo keyword incorrectly processed (LL ; Sep05)
;-

forward_function mpfit,mpfitfun

;	usage
ok='ok' & np=n_params() & nx=n_elements(x) & ny=n_elements(y) & na=n_elements(a)
if np lt 4 then ok='insufficient parameters' else $
 if nx eq 0 then ok='Input X undefined' else $
  if ny eq 0 then ok='Input Y undefined' else $
   if nx eq 1 then ok='Input not an array?' else $
    if nx ne ny then ok='X and Y(X) not compatible' else $
     if na eq 0 then ok='model parameters undefined'
if ok ne 'ok' then begin
  print,'Usage: mcerror,x,y,a,erru,errl,ysig=ysig,freeze=freeze,mciter=mciter,$'
  print,'       algo=algo,verbose=verbose,erra=erra,yfunc=yfunc,$'
  print,'       ,chithr=chithr, itmax = itmax,mcconf=mcconf,/dumb,mpfunc=mpfunc,$'
  print,'       jumpup=jumpup,jumpdn=jumpdn,svdthr=svdthr,$'
  print,'       funcs=funcs, function_name=function_name,$'
  print,'       type=type, missing=missing,/fwhm,/norm,betap=betap,/poisson'
  print,'  compute error bars on fit parameters'
  if np ne 0 then message,ok,/info
  return
endif
aN1=a

;	inputs
;	VERBOSE
v=0 & if keyword_set(verbose) then v=long(verbose(0)) > 1
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

;	ALGO
algo_type='LevMarq+SVD'			;default
algval=['LevMarq','Curvefit','MPFIT','tifsrhg'] ; 
if keyword_set(algo) then fitalgo=strtrim(algo(0),2) else fitalgo=algo_type                  
i0=-1L
i00=strpos(strlowcase(fitalgo),'lev',0) & i0 = i0 > i00
i01=strpos(strlowcase(fitalgo),'mar',0) & i0 = i0 > i01
i02=strpos(strlowcase(fitalgo),'lm',0) & i0 = i0 > i02
i03=strpos(strlowcase(fitalgo),'svd',0) & i0 = i0 > i03
i10=strpos(strlowcase(fitalgo),'idl',0) & i0 = i0 > i10
i11=strpos(strlowcase(fitalgo),'curve',0) & i0 = i0 > i11
i12=strpos(strlowcase(fitalgo),'mpfit',0) & i0 = i0 > i12 

if i00 ge 0 or i01 ge 0 or i02 ge 0 or i03 ge 0 then algo_type='LevMarq+SVD'
if i10 ge 0 or i11 ge 0 then algo_type='IDL-Curvefit'
if i12 ge 0 then algo_type ='MPFIT'
if algo_type eq 'IDL-Curvefit' and nf ne 0 then begin
  message,algo_type+' does not understand parameter freezing',/info
endif
if i0 lt 0 then message,$
	fitalgo+': not understood.  reverting to '+algo_type,/info
if v ge 1 then message,'using fitting algorithm: '+algo_type,/info

;	get the best-fit
if algo_type eq 'LevMarq+SVD' then fit_levmar,x,y,$
	a,yfunc,freeze=freeze,erra=erra,chisq=chisq,sig=sigy, _extra=e
if algo_type eq 'IDL-Curvefit' then yfunc=curve_fit(x,y,ywt,a,erra,chisq=chisq, _extra=e)

if algo_type eq 'MPFIT' then begin
         if keyword_set(mpfunc) then userfuncs=mpfunc else begin 
         message, 'mpfunc keyword required for mcerror calculation with MPFIT',/info
         return	;(VK: WAS quit)
         endelse 
         if v ge 5 then quiet = 0 else quiet = 1
         a=MPFITFUN(userfuncs,x,y,sigy,a,PERROR=erra,yfit=yfit,$
	    	parinfo=parinfo,functarg=functarg,quiet=quiet, _extra=e)
endif

;	outputs
erru=0.*erra & errl=erru

;	how many thawed parameters?
ot=where(ifit gt 0,mot)		;which parameters are thawed
if v gt 1 then message,'There are '+strtrim(mot,2)+' thawed parameters',/info
if mot gt 15 then begin
  message,'This will eons to complete and the answer will still be 42',$
	/info
  message,'type "x" NOW to return without further damage',/info
  wait,3 & c=get_kbrd(0)
  if strlowcase(c) eq 'x' then begin
    if v gt 5 then message,'quitting',/info
    return
  endif else message,'.. proceeding to find errors for all '+$
	strtrim(mot,2)+' thawed parameters',/info
endif


IF keyword_set(mciter) then nmc = mciter else nmc = 30  ; set mc iteration increment
IF keyword_set(mcconf) then conf=(mcconf<1d)/2d else conf= 0.68269/2d ;set confidence
IF conf eq 0 then conf= 0.68269/2d
ice = where(ifit eq 0,mice)                                  ; identify frozen parms for init config
fir = where(ifit eq 1,mfir)                                  ; identify thawed parms for init config
nmc = nmc*mfir >(2*(mfir+1)+1)  			; scale nmc by number of thawed parms
arg = DBLARR(n_elements(fir),nmc)		        ; initialize array to hold monte carlo results

IF algo_type eq 'LevMarq+SVD' THEN BEGIN
    for j = 0, nmc - 1, 1 do begin       ; initiate monte carlo sim  
        if v gt 2 then print, 'mc iter#',j+1
        randomdis=DOUBLE(randomn(seed, /normal, n_elements(x), 1))  ; create synthetic data set  
        yys = y + sigy*randomdis & aN = aN1
        fit_levmar, x, yys, aN, yfunc0, freeze=freeze, erra=erra, chisq=chisq, sig=sigy,verbose = 0, _extra=e
        arg[*,j] = aN(fir)
     endfor 
ENDIF

IF algo_type eq 'IDL-Curvefit' THEN BEGIN
     for j = 0, nmc - 1, 1 do begin      ; initiate monte carlo sim 
        if v gt 2 then print, 'mc iter#',j+1
        randomdis=DOUBLE(randomn(seed, /normal, n_elements(x), 1))  ; create synthetic data set  
        yys = y + sigy*randomdis & aN = aN1
        yfunc0=curve_fit(x,yys,ywt,aN,erra,chisq=chisq, _extra=e)  
        arg[*,j] = aN(fir)
     endfor
 ENDIF

IF algo_type eq 'MPFIT' THEN BEGIN	  
   for j = 0, nmc-1,1 do begin           ; initiate montecarlo sim 
       if v gt 2 then print, 'mc iter#',j+1
       randomdis=DOUBLE(randomn(seed, /normal, n_elements(x), 1))   ; create synthetic data set  
       yys = y + sigy*randomdis & aN = aN1
            results=MPFITFUN(userfuncs,x,yys,sigy,aN,PERROR=erra,yfit=yfit,$
	    	parinfo=parinfo,functarg=functarg,quiet=quiet,_extra=e)
	    arg[*,j] = results(fir) 
   endfor  
ENDIF 

;       Now get the upper/lower one sigma confidence bounds err
;       sort each row in arg (this contains mc results), after placing original fit in it 

For q =0, mfir-1 do begin 
  g = a(fir(q)) 
  auxar = [g, transpose(arg[q,*])] & argg = auxar(sort(auxar)) & a0w = where(argg eq g) 
  a0w = median(a0w) & intt = round(nmc*conf) & nmcc = nmc ;not nmc-1 becuase you added g 
  btbd = a0w - 2*intt > 0     ;lowest possible lower bound in IDL index space
  tpbd = a0w + 2*intt < nmcc  ;highest possible high bound in IDL index space
  testint = fltarr((tpbd-btbd-2*intt)+1);tricky? how many test intervals necessary  
  for j = btbd, tpbd-2*intt do testint(j-btbd) = argg(j+2*intt)-argg(j);calculate width of each interval
    jj = median(where(testint eq min(testint)))            ; identify smallest width
    eu = argg(jj+2*intt+btbd) & el = argg(jj+btbd) ; and corresponding bounds
errl(fir(q))=el & erru(fir(q))=eu 
endfor

;       ouput sims if necessary 
if keyword_set(mcsims) then  mcsims = arg 

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

if v ge 1 then print,''
if v ge 1 then print,'Done finding errors'

return
end
