function lnvandyk,s,b,Y,alfaS=alfaS,betaS=betaS,alfaB=alfaB,betaB=betaB,$
	Ybkg=Ybkg,asrc=asrc,abkg=abkg,chilike=chilike,verbose=verbose, _extra=e
;+
;function	lnvandyk
;	return the natural log of the joint posterior probability
;	distribution of the source and background built up using
;	the Poisson likelihood and Gamma-distribution priors for
;	the source and background,
;	  p(s,b|Y,I) = exp(-b*(betaB+1)+s*(betaS+1)) * 
;		       (b+s)^(Y) * b^(alfaB-1) * s^(alfaS-1) /
;		       (Y! * Gamma(alfaB) * Gamma(alfaS))
;	where Y are the observed counts, s,b are the model source and
;	background intensities, and alfaB,alfaS,betaB,betaS are the
;	hyperparameters of the gamma-distribution priors p(s|I) and
;	p(b|I) (see eqn 5 of van Dyk et al., 2001, ApJ 548, 224).
;	alfaB is generally set to the number of observed background
;	counts + 1, betaB is the ratio of the background area (or
;	exposure time) to the source area (or exposure time), and
;	alfaA and betaA are set based on prior knowledge (if any)
;	of the source intensity.
;
;syntax
;	lnp=lnvandyk(s,b,Y,alfaS=alfaS,betaS=betaS,alfaB=alfaB,betaB=betaB,$
;	Ybkg=Ybkg,asrc=asrc,abkg=abkg,/chilike,verbose=verbose)
;
;parameters
;	s	[INPUT; required] model source intensity
;	b	[INPUT; required] model background intensity
;	Y	[INPUT; required] observed counts
;		* output will be array of size [N(Y),N(S),N(B)]
;		  (beware that IDL automatically compresses trailing
;		  dimensions if any of them are equal to 1)
;
;keywords
;	alfaS	[INPUT; default=1] alpha parameter for source prior
;	betaS	[INPUT; default=0] beta parameter for source prior
;	alfaB	[INPUT; default=1] alpha parameter for background prior
;		* if YBKG is set, this is ignored completely
;	betaB	[INPUT; default=1] beta parameter for background prior
;		* if ASRC and ABKG are set, this is ignored completely
;	Ybkg	[INPUT; default=0] observed background counts
;		* if given, ALFAB is taken to be YBKG+1
;		* if even one element is specified, overrides ALFAB entirely
;	asrc	[INPUT; default=1] area (or exposure time) over which
;		source (and some background) counts are collected
;	abkg	[INPUT; default=1] area (or exposure time) over which
;		background counts are collected
;		* if ASRC _or_ ABKG is given, overrides BETAB entirely
;		  (unless, of course, they are 0 or less, in which case
;		  they are ignored)
;
;		* size of all of the above keywords are expanded or pruned
;		  to match that of Y.  i.e., there is freedom to specify
;		  a separate prior for each data point, but not for each
;		  model grid point.
;
;	chilike	[INPUT] if set, returns -2*SUM_{Y} (ln(p(s,b|Y)))
;	verbose	[INPUT] controls chatter
;
;	_extra	[JUNK] here only to prevent crashing the program
;
;common blocks
;
;history
;	vinay kashyap (Aug01)
;-

;	usage
ok='ok' & np=n_params() & ns=n_elements(s) & nb=n_elements(b) & nY=n_elements(Y)
if np lt 3 then ok='Insufficient parameters' else $
 if ns eq 0 then ok='S is undefined' else $
  if nb eq 0 then ok='B is undefined' else $
   if nY eq 0 then ok='Y is undefined'
if ok ne 'ok' then begin
  print,'Usage: lnp=lnvandyk(s,b,Y,alfaS=alfaS,betaS=betaS,alfaB=alfaB,$'
  print,'       betaB=betaB,Ybkg=Ybkg,asrc=asrc,abkg=abkg,/chilike,$'
  print,'       verbose=verbose)'
  print,'  return ln(p(s,b|Y)) assuming gamma-priors on s and b, as in
  print,'  van Dyk et al. (2001; ApJ 548, 224)'
  if np ne 0 then message,ok,/info
  return,-1L
endif

;	keywords
chi=1. & if keyword_set(chilike) then chi=-2.
vv=0 & if keyword_set(verbose) then vv=long(verbose[0])>1
;
nas=n_elements(alfaS) & nbs=n_elements(betaS)
nab=n_elements(alfaB) & nbb=n_elements(betaB)
nyb=n_elements(Ybkg) & nasa=n_elements(asrc) & naba=n_elements(abkg)
aS=fltarr(ny)+1. & bS=fltarr(ny) & aB=aS & bB=bS
yB=fltarr(ny) & sar=fltarr(ny)+1. & bar=sar
;
if nas gt 0 then begin
  aS[*]=alfaS[nas-1L]
  if nas lt ny then aS[0L:nas-1L]=alfaS[*] else aS[*]=alfaS[0L:ny-1L]
endif
;
if nbs gt 0 then begin
  bS[*]=betaS[nbs-1L]
  if nbs lt ny then bS[0L:nbs-1L]=betaS[*] else bS[*]=betaS[0L:ny-1L]
endif
;
if nab gt 0 then begin
  aB[*]=alfaB[nab-1L]
  if nab lt ny then aB[0L:nab-1L]=alfaB[*] else aB[*]=alfaB[0L:ny-1L]
endif
;
if nbb gt 0 then begin
  bB[*]=betaS[nbb-1L]
  if nbb lt ny then bB[0L:nbb-1L]=betaS[*] else bB[*]=betaS[0L:ny-1L]
endif
;
if nyb gt 0 then begin
  if vv gt 10 and nab gt 0 then message,'overriding ALFAB with YBKG+1',/info
  aB[*]=Ybkg[nyb-1L]+1.
  if nyb lt ny then aB[0L:nyb-1L]=Ybkg[*] else aB[*]=Ybkg[0L:ny-1L]
endif
;
if nasa gt 0 then begin
  sar[*]=asrc[nasa-1L]
  if nasa lt ny then sar[0L:nasa-1L]=asrc[*] else sar[*]=asrc[0L:ny-1L]
endif
if naba gt 0 then begin
  bar[*]=abkg[naba-1L]
  if naba lt ny then bar[0L:naba-1L]=asrc[*] else bar[*]=asrc[0L:ny-1L]
endif
if nasa gt 0 or naba gt 0 then begin
  if vv gt 10 and nbb gt 0 then message,'overriding BETAB with ABKG/ASRC',/info
  oo=where(bar gt 0 and sar gt 0,moo)
  if moo gt 0 then bB[oo]=bar[oo]/sar[oo]
endif

;	define the output
lnp=fltarr(ny,ns,nb)

;	populate the output
oy0=where(y eq 0,moy0) & oy1=where(y gt 0,moy1) & oyx=where(y lt 0,moyx)
ob0=where(aB eq 0,mob0) & ob1=where(aB gt 0,mob1) & obx=where(aB lt 0,mobx)
os0=where(aS eq 0,mos0) & os1=where(aS gt 0,mos1) & osx=where(aS lt 0,mosx)
den1=!values.F_INFINITY & if moy1 gt 0 then den1[oy1]=lngamma(y[oy1]+1)
den2=!values.F_INFINITY & if mob1 gt 0 then den2[ob1]=lngamma(aB[ob1])
den3=!values.F_INFINITY & if mos1 gt 0 then den3[os1]=lngamma(aS[os1])
;
for j=0L,ns-1L do begin			;{for each source model point

  num4=0.*y
  if s[j] gt 0 then num4=(aS-1.)*alog(s[j]) else begin
    num4[*]=-!values.F_INFINITY & if mos1 gt 0 then num4[os1]=0.
  endelse

  for k=0L,nb-1L do begin		;{for each background model point

    num1=-b[k]*(bB+1.)-s[j]*(bS+1.)

    num2=0.*y
    if s[j]+b[k] gt 0 then num2=y*alog(b[k]+s[j]) else begin
      if moy0 gt 0 then num2[oy0]=0.
      if moy1 gt 0 then num2[oy1]=-!values.F_INFINITY
    endelse

    num3=0.*y
    if b[k] gt 0 then num3=(aB-1.)*alog(b[k]) else begin
      num3[*]=-!values.F_INFINITY & if mob1 gt 0 then num3[ob1]=0.
    endelse

    lnp[*,j,k]=chi*(num1+num2+num3+num4-den1-den2-den3)

  endfor				;K=0,NB-1}
endfor					;J=0,NS-1}

if keyword_set(chilike) then begin
  if ny gt 1 then begin
    tmp=fltarr(ns,nb)
    for i=0L,ns-1L do for j=0L,nb-1L do tmp[i,j]=total(lnp[*,i,j])
    lnp=tmp
  endif
endif

return,lnp
end
