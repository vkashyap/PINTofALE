pro mcmc_burn,x,y,ain,rngpar,aout,siga,nburn=nburn,rhat=rhat,qhat=qhat,$
	qfrac=qfrac,rthresh=rthresh,qthresh=qthresh,nbatch=nbatch,$
	aburn=aburn, _extra=e
;+
;procedure	mcmc_burn
;	burns in the parameters in a Markov-Chain Monte Carlo sequence
;
;syntax
;	mcmc_burn,x,y,ain,rngpar,aout,siga,nburn=nburn,rhat=rhat,qhat=qhat,$
;	qfrac=qfrac,rthresh=rthresh,qthresh=qthresh,$
;	sigy=sigy,funcs=funcs,fnprob=fnprob,thaw=thaw,/singly,$
;	sampar=sampar,sclpar=sclpar,seed=seed,testtyp=testtyp,ties=ties,$
;	FUNCS;	type=type,/fwhm,/norm,betap=betap,vrot=vrot,angle=angle,$
;		phase=phase,group=group,delp=delp,missing=missing,$
;	FNPROB;	ulim=ulim,/chi2,/binom,/cash,/castor
;
;parameters
;	x	[INPUT; required] data points
;	y	[INPUT; required] Y(X)
;		* sizes of X and Y _must_ match
;	ain	[INPUT; required] initial values of the parameters of
;		the model to fit to Y(X)
;		* parameters should be 1-D array
;		* if 2-D, the size of the 2nd dimension describes the
;		  number of chains that are concurrently running
;	rngpar	[INPUT; required] the allowed range for each parameter
;		as a 2-D array of size (N(PAR),2)
;		* if 1-D array, range assumed to be symmetrical
;		  -- additive if +ve
;		  -- abs value multiplicative if -ve
;		* if size does not match, assumed to be uniformly (-100)
;	aout	[OUTPUT; required] values of the parameters at conclusion
;		of burn-in phase -- will be same size as AIN
;	siga	[OUTPUT; required] the width of the allowed deviation
;		to generate new parameter values from old values.  this
;		is usually the standard deviation of the samples, unless
;		abs(kurtosis)>0.2, in which case SIGA is set to -ve of
;		the maximum possible deviation.
;
;keywords
;	nburn	[INPUT] maximum number of steps to go through
;		* default is 100
;	rhat	[OUTPUT] ratio of rms(inter-chain) to rms(intra-chain)
;		for each parameter -- the chains have converged if
;		RHAT drops to ~ 1
;		* will be set to 0 if there is only one chain
;	qhat	[OUTPUT] a measure of the stability of the persistent
;		distribution of the parameters, which indicates the
;		convergence of a given chain.  this is simply the
;		ratio of the product of an older frequency distribution
;		of a parameter with a newer one, normalized by number
;		of samples, and averaged over all the parameters.
;		* will return a value for each chain
;		* the chain will have converged if QHAT goes up to ~ 1
;	qfrac	[INPUT] the fractions of the parameter chain to be used
;		to define the "older" and "newer" frequency distributions
;		in the calculation of QHAT.
;		* must be composed of numbers >0 and <1, and must
;		  sum to <1
;		 -- if .le.0 or .ge.1, or is not set, is set to 0.25
;		* should be 2-element vector describing the fractional
;		  lengths of the chains to consider: QFRAC[0] is for the
;		  "new" distribution and runs backwards from newest, and
;		  QFRAC[1] is for the "old" distribution, and runs backwards
;		  from the end of "old"
;		  -- if QFRAC[0]+QFRAC[1] > 1, then QFRAC[1] is set to
;		     1-QFRAC[0]
;	rthresh	[INPUT] the threshold value of RHAT below which it can be
;		assumed that the chains have converged
;		* default is 1.5
;	qthresh	[INPUT] the threshold value of QHAT above which it can be
;		assumed that the chain has converged
;		* default is 0.95
;		* notes:
;		  -- generally, both RTHRESH and QTHRESH must be satisfied,
;		     and the latter for all chains, for burning to stop
;		     before NBURN is reached
;		  -- if the number of chains is 1, RTHRESH is ignored
;	nbatch	[INPUT] caught and discarded, because it will be set to 1
;		in call to MCMC_STEP
;	aburn	[OUTPUT] all the parameters looked at in each step
;		* is of size {NPAR,NCHAIN,NBURN}
;	_extra	[INPUT ONLY] pass defined keywords to subroutines:
;		sigy	[MCMC_STEP]
;		funcs	[MCMC_STEP]
;		fnprob	[MCMC_STEP]
;		thaw	[MCMC_STEP]
;		singly	[MCMC_STEP]
;		sampar	[MCMC_STEP]
;		sclpar	[MCMC_STEP]
;		seed	[MCMC_STEP]
;		testtyp	[MCMCM]
;		ties	[ADJUSTIE]
;		*	[FUNCS]
;		*	[FNPROB]
;
;subroutines
;	ADJUSTIE
;	ARITHTOGRAM()
;	LOGIT()
;	MCMC_STEP()
;	MCMCM()
;	UNLOGIT()
;	-FUNCS-
;	-FNPROB()-
;
;history
;	vinay kashyap (Jun05)
;-

;	usage
ok='ok' & np=n_params() & nx=n_elements(x) & ny=n_elements(y)
sza=size(ain) & nsza=n_elements(sza) & szr=size(rngpar) & nszr=n_elements(szr)
if np lt 6 then ok='Insufficient parameters' else $
 if nx eq 0 then ok='X is undefined' else $
  if ny eq 0 then ok='Y(X) is undefined' else $
   if nx ne ny and nx+1L ne ny then ok='X and Y(X) are incompatible' else $
    if sza[0] eq 0 then ok='input parameters are undefined' else $
     if sza[0] gt 2 then ok='cannot understand parameter arrays >2-D'
if ok ne 'ok' then begin
  print,'Usage: mcmc_burn,x,y,ain,rngpar,aout,siga,nburn=nburn,$'
  print,'       rhat=rhat,qhat=qhat,qfrac=qfrac,rthresh=rthresh,qthresh=qthresh,$'
  print,'       sigy=sigy,funcs=funcs,fnprob=fnprob,nbatch,thaw=thaw,/singly,$'
  print,'       sampar=sampar,sclpar=sclpar,seed=seed,testtyp=testtyp,ties=ties,$'
  print,'       type=type,/fwhm,/norm,betap=betap,vrot=vrot,angle=angle,$'
  print,'       phase=phase,group=group,delp=delp,missing=missing,$'
  print,'       ulim=ulim,/chi2,/binom,/cash,/castor'
  print,'  burns in the parameters in a Markov-Chain Monte Carlo sequence'
  if np ne 0 then message,ok,/informational
  return
endif

;	check keywords
  ;
;  number of parameters
npar=sza[1]
  ;
;  number of chains
if sza[0] eq 2 then nchain=sza[2] else nchain=1
  ;
;  nburn
mburn=100L & if keyword_set(nburn) then mburn=long(nburn[0])>1
  ;
;  rthresh
rthr=1.5 & if keyword_set(rthresh) then rthr=abs(rthresh[0])
if rthr le 0 then rthr=1.5
  ;
;  qthresh
qthr=0.95 & if keyword_set(qthresh) then qthr=abs(qthresh[0])
if qthr le 0 then qthr=0.95
  ;
;  qfrac
fracq=[0.25,0.25] & nq=n_elements(qfrac)
if nq eq 1 then if qfrac[0] gt 0 and qfrac[0] lt 1 then fracq=qfrac[0]*[1,1]
if nq eq 2 then begin
  fracq=qfrac
  if fracq[0] le 0 or fracq[0] ge 1 then fracq[0]=0.25
  if fracq[1] le 0 or fracq[1] ge 1 then fracq[1]=0.25
endif
if total(fracq) gt 1 then fracq[1]=1.-fracq[0]
  ;
;  range of PAR
szr=size(rngpar) & userng=0 & parrng=reform([par/100.,par*100.],npar,2)
case szr[0] of
  1: begin
    if szr[1] eq npar then begin
      for i=0L,npar-1L do begin
	if rngpar[i] gt 0 then parrng[i,*]=par[i]+rngpar[i]*[-1,1]
	if rngpar[i] lt 0 then parrng[i,*]=par[i]*$
		[1./abs(rngpar[i]),abs(rngpar[i])]
      endfor
    endif
  end
  2: if szr[1] eq npar and szr[2] eq 2 then parrng=rngpar
  else: ;nothing
endcase

;	run the chains
iburn=0L & sigpar=-abs(parrng[*,1]-parrng[*,0])
allpar=fltarr(npar,nchain,mburn) & qhat=fltarr(nchain)
rflag=0 & qflag=0
while iburn lt mburn do begin		;{burn in
  
  jburn=iburn
  aout=mcmc_step(x,y,ain,sigpar,rngpar=parrng,nbatch=1, _extra=e)
  for j=0L,nchain-1L do allpar[*,j,iburn]=aout[*,j]

  ; compute RHAT
  if nchain gt 1 and jburn gt 2 then begin
    sigin=fltarr(npar) & meain=sigin
    for ip=0L,npar-1L do begin
      tmps=fltarr(nchain) & tmpm=tmps
      for ic=0L,nchain-1L do begin
	tmps[ic]=stddev(allpar[ip,ic,0L:jburn])
	tmpm[ic]=mean(allpar[ip,ic,0L:jburn])
      endfor
      sigin[ip]=sqrt(total(tmps^2))
      sigou[ip]=stddev(tmpm)
    endfor
    sdevin=sqrt(total(sigin^2)) & sdevou=sqrt(total(sigou^2))
    if sdevin gt 0 then rhat=sdevou/sdevin
    if rhat le rthr then rflag=1
  endif
  if nchain eq 1 then rflag=1

  ; compute QHAT
  nn=jburn+1L
  i1=long(nn*(1-fracq[0]))
  i0=long(nn*(1-fracq[0]-fracq[1]))
  if nn-i1 gt 2 and i1-i0 gt 2 then begin
    qhat=fltarr(nchain)
    for ic=0L,nchain-1L do begin
      tmpqhat=0. & tmpnorm=0.
      for ip=0L,npar-1L do begin
        tmppre=allpar[ip,ic,i0:i1]
        tmppost=allpar[ip,ic,i1:jburn]
	nitem=(i1-i0+1L) < (jburn-i1+1L)
	nbin=((nitem/10L) > 5L) < 101L	;on avg, about 10 items/bin
        h1=arithtogram(tmppre,tmppre,junk,'*',nbin=nbin,$
	  xmin=parrng[ip,0],xmax=parrng[ip,1],/nonorm,$
	  _extra=e)
        h2=arithtogram(tmppost,tmppost,junk,'*',nbin=nbin,$
	  xmin=parrng[ip,0],xmax=parrng[ip,1],/nonorm,$
	  _extra=e)
        h12=arithtogram(tmppre,tmppost,junk,'*',nbin=nbin,$
	  xmin=parrng[ip,0],xmax=parrng[ip,1],/nonorm,$
	  _extra=e)
	tmphat=total(abs(h1-h2)/h12^2,/nan)
	if tmphat gt 0 then begin
	  tmpqhat=tmpqhat+tmphat
	  tmpnorm=tmpnorm+1.
	endif
      endfor
      if tmpnorm gt 0 then qhat[ic]=tmpqhat/tmpnorm
    endfor
    ok=where(qhat le qthr,mok)
    if mok eq 0 then qflag=1
  endif

  ; stopping rules
  iburn=iburn+1L
  ok='ok'
  if iburn eq mburn then ok='burnt up all the firewood'
  if rflag eq 1 and qflag eq 1 then $
    ok='RHAT drops below threshold and QHAT exceeds threshold'
  if ok ne 'ok' then begin
    message,ok,/informational
    iburn=mburn
  endif
endwhile				;IBURN<MBURN}

;	recompute SIGA
nn=jburn & aburn=allpar[*,*,0L:nn] & mm=long(0.5*nn)
siga=fltarr(npar)
for ip=0L,npar-1L do begin
  tmp=aburn[ip,*,mm:nn]
  tmps=stddev(tmp,/nan) & tmpk=abs(kurtosis(tmp,/nan))
  tmpd=max(tmp)-min(tmp)
  if tmpk le 0.2 then siga[ip]=tmps else siga[ip]=-tmpd
endfor

return
end
