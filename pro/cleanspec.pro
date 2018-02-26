function cleanspec,spec,bspec,bounds,bkgscl,srcscl=srcscl,srcest=srcest,$
	bkgest=bkgest,getave=getave,clev=clev,verbose=verbose, _extra=e
;+
;function	cleanspec
;	return a "clean" spectrum, one with the best estimates of the
;	source count without explicitly subtracting the background but
;	rather by marginalizing over the background model intensities,
;	as in van Dyk, D.A., Connors, A., Kashyap, V.L., & Siemiginowska,
;	A., 2001, ApJ 548, 224.
;
;syntax
;	srcct=cleanspec(spec,bspec,bounds,bkgscl,srcscl=srcscl,srcest=srcest,$
;	bkgest=bkgest,/getave,clev=clev,verbose=verbose, nsgrid=nsgrid,$
;	srcmax=srcmax,srcmin=srcmin,jmax=jmax,/double,k=k,eps=eps)
;
;parameters
;	spec	[INPUT; required] array of counts in the source region
;	bspec	[INPUT; required] array of counts in the background region
;		* if size is smaller than SPEC, gets filled out with the
;		  value of the last element of the array
;	bounds	[OUTPUT] array of size [2,N(SPEC)] containing the lower [0,*]
;		and upper [1,*] bounds on the cleaned source count, at a
;		level defined by CLEV
;	bkgscl	[INPUT] the area (or exptime) in which BSPEC were collected
;		* if not given, assumed to be same as SRCSCL
;		* if size is smaller than SPEC, fills out with last element
;
;keywords
;	srcscl	[INPUT] area (or exptime) in which SPEC were collected
;		* default is 1.0
;		* if size is smaller than SPEC, fills out with last element
;	srcest	[INPUT] a priori information about the source strength --
;		this value is used to define ALPHA_S (=SRCEST+1) in the
;		gamma prior and BETA_S (1 if SRCEST is non-zero, 0 otherwise)
;		* if illegal (i.e., < 0), then defaults to the
;		  non-informative prior of ALPHA_S=1,BETA_S=0
;	bkgest	[INPUT] a priori information about the background strength --
;		this is used to define ALPHA_B (=BKGEST+1)
;		* used iff BSPEC is 0 or -ve
;		* default value is <BSPEC> or <SPEC>*BKGSCL/SRCSCL,
;		  if the former is zero
;	getave	[INPUT] if set, returns the mean value computed from
;		the posterior probability distribution rather than
;		the mode.  beware that this will be off by 1 from
;		what one may normally expect from gaussian intuition.
;	clev	[INPUT] level at which to determine the bounds on the
;		source intensities
;		* default is 0.68
;		* if < 0, abs value is used
;		* if > 1 and < 100, then assumed to be given as a percentage
;		* if > 100, then 1-1/CLEV is used as the true value
;	verbose	[INPUT] controls chatter
;	_extra	[INPUT ONLY] pass defined keywords to subroutines
;		-- POST_SRCGAMMA: NSGRID,SRCMAX,SRCMIN
;		-- QROMB: JMAX,DOUBLE,K,EPS
;
;subroutines
;	POST_SRCGAMMA
;	PROB_GAMMADIST
;	KILROY
;
;history
;	vinay kashyap (Aug01)
;-

;	usage
ok='ok' & np=n_params() & ns=n_elements(spec) & nb=n_elements(bspec)
if np lt 2 then ok='Insufficient parameters' else $
 if ns eq 0 then ok='SPEC is undefined' else $
  if nb eq 0 then ok='BSPEC is undefined'
if ok ne 'ok' then begin
  print,'Usage: srcct=cleanspec(spec,bspec,bounds,bkgscl,srcscl=srcscl,$'
  print,'       srcest=srcest,bkgest=bkgest,/getave,clev=clev,verbose=verbose,$'
  print,'       nsgrid=nsgrid,srcmax=srcmax,srcmin=srcmin,jmax=jmax,$'
  print,'       /double,k=k,eps=eps)'
  print,'  return a "clean" source spectrum determined by marginalizing'
  print,'  the background intensities'
  if np ne 0 then message,ok,/info
  return,-1L
endif

;	inputs and keywords
vv=0 & if keyword_set(verbose) then vv=long(verbose[0]) > 1
;
sct=[spec[*]] & bct=fltarr(ns) & sscl=fltarr(ns)+1.0 & bscl=sscl
nss=n_elements(srcscl) & nbs=n_elements(bkgscl)
nse=n_elements(srcest) & nbe=n_elements(bkgest)
;
if nb eq ns then bct=[bspec[*]] else begin
  bct[*]=bspec[nb-1L]
  if nb le ns then bct[0L:nb-1L]=bspec[*] else bct[*]=bspec[0L:ns-1L]
endelse
;
if nss gt 0 then begin
  sscl[*]=srcscl[nss-1L]
  if nss le ns then sscl[0L:nss-1L]=srcscl[*] else sscl[*]=srcscl[0L:ns-1L]
endif
o0=where(sscl eq 0,mo0) & if mo0 ne 0 then sscl[o0]=1.0
;
if nbs gt 0 then begin
  bscl[*]=bkgscl[nbs-1L]
  if nbs le ns then bscl[0L:nbs-1L]=bkgscl[*] else bscl[*]=bkgscl[0L:ns-1L]
endif else bscl=sscl
o0=where(bscl eq 0,mo0) & if mo0 ne 0 then bscl[o0]=1.0
;
sprior=fltarr(ns)
if nse gt 0 then begin
  sprior[*]=srcest[nse-1L]
  if nse le ns then sprior[0L:nse-1L]=srcest[*] else sprior[*]=srcest[0L:ns-1L]
endif
;
bprior=bct & bp=bct & dbp=bprior
if nbe gt 0 then begin
  o0=where(bprior eq 0,mo0)
  if mo0 gt 0 then dbp[*]=total(bct)/nb
  if dbp[0] eq 0 then dbp=total(sct)*bscl/sscl
  if nbe le ns then bp[0L:nbe-1L]=bkgest[*] else bp[*]=bkgest[0L:ns-1L]
  oo=where(bprior eq 0 and bp gt 0,moo)
  if moo gt 0 then bprior[oo]=bp[oo]
  oo=where(bprior eq 0 and dbp gt 0,moo)
  if moo gt 0 then bprior[oo]=dbp[oo]
endif
;
credlev=0.68 & if keyword_set(clev) then credlev=0.0+clev[0]
if credlev lt 0 then credlev=abs(credlev)
if credlev gt 1 and credlev lt 100 then credlev=credlev/100.
if credlev gt 100 then credlev=1.0D - 1.0D/credlev

;	and now determine the output
sspec=fltarr(ns) & bounds=fltarr(2,ns)
for i=0L,ns-1L do begin			;{for each bin
  
  yct=sct[i] & ybkg=bprior[i] & asrc=sscl[i] & abkg=bscl[i]
  agsrc=sprior[i]+1. & if nse gt 0 then bgsrc=1. else bgsrc=0.
  if sprior[i] lt 0 or finite(sprior[i]) eq 0 then begin
    agsrc=1. & bgsrc=0.			 ;non-informative prior
  endif
  if vv gt 0 then kilroy

  ;	detour for a time saver
  oo=where(yct eq sct and ybkg eq bprior and $
	agsrc eq sprior+1 and asrc eq sscl and abkg eq bscl,moo)

  if oo[0] ge 0 and oo[0] lt i then begin	;(aha!

    ; this means there was at least one previous calculation performed
    ; for the same set of parameters

    if vv gt 3 then kilroy,dot='('+strtrim(yct,2)+','+strtrim(ybkg,2)+$
      ';'+strtrim(asrc,2)+','+strtrim(abkg,2)+')'

    sspec[i]=sspec[oo[0]] & bounds[*,i]=bounds[*,oo[0]]

  endif else begin		;)(need to compute..

    if vv gt 2 then kilroy,dot='('+strtrim(yct,2)+','+strtrim(ybkg,2)+$
      ';'+strtrim(asrc,2)+','+strtrim(abkg,2)+')'

    psrc=post_srcgamma(yct,ybkg,pstr,asrc=asrc,abkg=abkg,agsrc=agsrc,$
	bgsrc=bgsrc,clev=credlev,verbose=vv, _extra=e)

    if keyword_set(getave) then begin
      sspec[i]=pstr.MEAN & bounds[*,i]=pstr.CREDLEV
    endif else begin
      sspec[i]=pstr.MODE & bounds[0,i]=pstr.LLEV & bounds[1,i]=pstr.ULEV
    endelse

  endelse			;)

  if vv gt 1 then kilroy,dot=strtrim(ns-i,2)
  if vv gt 10 then kilroy,dot='['+strtrim(sspec[i],2)+']'

endfor					;I=0,NS-1}

return,sspec
end
