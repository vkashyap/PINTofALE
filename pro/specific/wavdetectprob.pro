function wavdetectprob,src,bkg,signi,wscale=wscale,$
	psfwdt=psfwdt,verbose=verbose,$
	nbkgsim=nbkgsim,prbarr=prbarr,bkgarr=bkgarr,$
	_extra=e
;+
;function wavdetectprob
;	Computes the probability that a source of given strength
;	will be detected at a specified significance for a given
;	background.  Uses the analytical and numerically derived
;	formulae in Appendix B of the Wavdetect paper,
;	Freeman, P.E., Kashyap, V., Ronser, R., \& Lamb, D.Q.,
;	2002, ApJS, 138, 185
;
;syntax
;	prob=wavdetectprob(src,bkg,signi,wscale=wscale,psfwdt=psfwdt,$
;	nbkgsim=nbkgsim,prbarr=prbarr,bkgarr=bkgarr,verbose=verbose)
;
;parameters
;	src	[INPUT; required] intrinsic source counts for which to
;		compute probability of detection
;		* may be an array
;		* ignored if 0 or -ve
;	bkg	[INPUT; required] background counts per pixel
;		* must be scalar
;		* quits if 0 or -ve
;	signi	[INPUT] significance at which to compute detection prob
;		* must be scalar
;		* default is 1e-6
;		* if <0, abs(SIGNI) is assumed
;		* if =1, taken to be 1e-6
;		* if >1, reciprocal is used
;
;keywords
;	wscale	[INPUT] wavelet scale
;		* default is 4
;	psfwdt	[INPUT] sigma of Gaussian that describes the PSF
;		* default is 1
;	nbkgsim	[INPUT] number of realizations of the background
;		* default is 0
;		* if set, generates Poisson deviates of q (the
;		  background under the wavelet) and computes
;		  the detection probability for each case
;	prbarr	[OUTPUT] if NBKGSIM>0, then the output arrays
;		for each simulation are written into an array
;		of size (NSRC,NBKGSIM+1), with the first row
;		corresponding to the input BKG, i.e., identical
;		to the primary output for the case NBKGSIM=0
;		* if NBKSIM>0, the primary output will be set to the
;		  average of PRBARR[*,1:*] across the simulations
;	bkgarr	[OUTPUT] background values for which PRBARR rows are computed
;	verbose	[INPUT] controls chatter
;	_extra	[JUNK] here only to control chatter
;
;example
;	.run wavdetectprob
;
;history
;	Vinay Kashyap (2014apr22)
;	bug fix (Nick W; 2014apr24)
;	added keywords NBKGSIM, PRBARR, BKGARR; changed IGAMMA call
;	  to explicit summation (VK; 2014apr25)
;-

;	usage
ok='ok' & np=n_params() & ns=n_elements(src) & nb=n_elements(bkg) & nt=n_elements(signi)
if np lt 2 then ok='Insufficient parameters' else $
 if ns eq 0 then ok='SRC is not defined' else $
  if nb eq 0 then ok='BKG is not defined' else $
   if nb gt 1 then ok='BKG must be a scalar' else $
    if bkg[0] le 0 then ok='BKG must be +ve definite' else $
     if nt gt 1 then ok='SIGNI must be a scalar'
if ok ne 'ok' then begin
  print,'Usage: prob=wavdetectprob(src,bkg,signi,wscale=wscale,psfwdt=psfwdt,$'
  print,'  nbkgsim=nbkgsim,prbarr=prbarr,bkgarr=bkgarr,verbose=verbose)'
  print,'       compute probability of detection for wavdetect'
  if np ne 0 then message,ok,/informational
  return,-1L
endif

;	inputs
vv=0L & if keyword_set(verbose) then vv=long(verbose[0])>1L
;
ss=src
;
bb=bkg[0]
;
thr=signi[0]
if thr lt 0 then thr=abs(thr)
if thr eq 1 then thr=1e-6
if thr gt 1 then thr=1.D / thr
;
nw=n_elements(wscale)
wx=4. & wy=4.
if nw ne 0 then begin
  wx=0.+abs(wscale[0]) & wy=wx
  if nw gt 1 then wy=0.+abs(wscale[1])	;allow for asymmetric wavelets sometime in future
endif
wy=wx	;but ignore asymmetric wavelets for now
;
npsf=n_elements(psfwdt)
wp=1.
if npsf ne 0 then begin
  wp=0.+abs(psfwdt[0])
endif
;
nsim=0L & nbkg=n_elements(nbkgsim)
if nbkg eq 0 then nsim=0L else $
 if nbkg ge 1 then nsim=long(nbkgsim[0])>0

;	output
pp=0.*ss	;is in same form as input
prbarr=fltarr(ns,nsim+1L)

;	compute q
qq=bb*2.*!pi*wx*wy/exp(1.)
lqq = alog10(qq)

;	compute expected threshold correlation value from Eqns B2 and B3
A_lo = 0.00462 & B_lo = 0.0661
C_lo = -0.0154*(alog10(thr))^2 - 0.252*alog10(thr) - 0.031
A_mid = 0.00182*(alog10(thr))^3+0.0279*(alog10(thr))^2+0.158*alog10(thr)+0.607
if alog10(thr) lt -7 then B_mid = -0.064*alog10(thr)+0.612 else $
 if alog10(thr) le -2 then B_mid = -0.064*alog10(thr)+0.612-0.00085*(alog10(thr)+7.)^3.5 else $
  if alog10(thr) le -1 then B_mid = -0.064*alog10(thr)+0.612-0.00015*(alog10(thr)+7.)^4.75 else $
   message,'SIGNI is too high, this regime is uncalibrated'
if alog10(thr) ge -7 then A_hi = -0.509*alog10(thr) + 1.897 - 0.00172*(alog10(thr)+7)^3.606 else $
	A_hi = -0.509*alog10(thr) + 1.897
B_hi = -1.115*alog10(thr) - 1.038
;
if lqq lt -0.5 then Corthr = 10.^(A_lo * lqq^2 + B_lo * lqq + C_lo) else $
  if lqq lt 1 then Corthr = 10.^(A_mid * lqq + B_mid) else $
    Corthr = (A_hi * sqrt(qq) + B_hi)

;	this corresponds to an amplitude of (from Eqn 6)
corcor = (wx * wy) / sqrt((wp^2+wx^2)*(wp^2+wy^2))
corcor = corcor * (2. - ((1.0*wp)^2/(wp^2+wx^2)) - ((1.0*wp)^2/(wp^2+wy^2)) )
if corcor le 0 then message,'BUG!'
ctampl = Corthr / corcor
;	which corresponds to a count of
cttot = ctampl ;* 2. * !pi * wp^2 (because the Gaussian is already normalized)

;	what is the probability of seeing >CTTOT counts if the true intensity is SS?
;for i=0L,ns-1L do if ss[i] gt 0 then pp[i] = 1.D - igamma(ss[i],floor(cttot+1))
;for i=0L,ns-1L do if ss[i] gt 0 then pp[i] = igamma(floor(cttot+1),ss[i])
kk=lindgen(floor(cttot)+1L) & lnk=lngamma(kk+1)
for i=0L,ns-1L do begin
  if ss[i] gt 0 then begin
    tmp=kk*alog(ss[i])-ss[i]-lnk & tmpmax=max(tmp,/nan) & tmp2=exp(tmp-tmpmax)
    pp[i]=1.D - exp(tmpmax)*total(tmp2,/nan)
  endif
endfor

if vv gt 10 then print,'wscale,psfwdt,bkg,q,logq,Corthr,corcor,cttot',wx,wp,bkg,qq,lqq,Corthr,corcor,cttot

;	compute error
prbarr[*,0]=pp & bkgarr=[bb]
if nsim gt 0 then begin
  bkgarr=randomu(seed,nsim,poisson=qq)*(float(bb)/qq)
  for i=1L,nsim do prbarr[*,i]=wavdetectprob(src,bkgarr[i-1L],signi,wscale=wscale,psfwdt=psfwdt)
  if nsim gt 1 then pp=total(prbarr[*,1:*],2)/nsim
endif

if vv gt 1000 then stop,'HALTing; type .CON to continue'

return,pp
end

;===============================================================================
;example calling sequence
;===============================================================================

if n_elements(src) eq 0 then src=findgen(200)+1
if n_elements(src) eq 0 then src=10.^(findgen(31)*0.1-1.5)
if not keyword_set(bkg) then bkg=0.27
if not keyword_set(bkg) then bkg=0.1
if not keyword_set(signi) then signi=1e-6
if not keyword_set(wscale) then wscale=8
if not keyword_set(wscale) then wscale=4
if not keyword_set(psfwdt) then psfwdt=3.66
if not keyword_set(psfwdt) then psfwdt=1
if n_elements(nbkgsim) eq 0 then nbkgsim=1000
if not keyword_set(nbkgsim) then nbkgsim=0
if not keyword_set(verbose) then verbose=11
if not keyword_set(verbose) then verbose=1

prob=wavdetectprob(src,bkg,signi,wscale=wscale,psfwdt=psfwdt,nbkgsim=nbkgsim,prbarr=prbarr,bkgarr=bkgarr,verbose=verbose)

plot,src,prob,/xl,/xs,/ys,xtitle='source counts',ytitle='detection probability',$
	subtitle='background='+strtrim(bkg,2)+' [ct/pix] ; scale='+strtrim(wscale,2)

end
