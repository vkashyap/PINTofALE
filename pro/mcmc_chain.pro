function mcmc_chain,x,y,parval,parsig,parinit=parinit,parrng=parrng,$
	ulim=ulim,useMH=useMH,ysig=ysig,freeze=freeze,funcs=funcs,$
	type=type,ydc=ydc,$
	arfil=arfil,areff=areff,rmfil=rmfil,ties=ties,xfilter=xfilter,perbin=perbin,$
	nsim=nsim,nburn=nburn,onepar=onepar,verbose=verbose,$
	bgx=bgx,bgy=bgy,bparval=bparval,bparsig=bparsig,bparini=bparini,$
	bparrng=bparrng,bysig=bysig,bkgscal=bkgscal,bfreeze=bfreeze,$
	btype=btype,bties=bties,$
	mcmc_trace=mcmc_trace,mcmc_prb=mcmc_prb,$
	mcmc_flux=mcmc_flux,mcmc_counts=mcmc_counts,$
	mcmc_btrace=mcmc_btrace,mcmc_bprb=mcmc_bprb,$
	mcmc_bflux=mcmc_bflux,mcmc_bcounts=mcmc_bcounts,$
	solnup=solnup,solndn=solndn,solnpr=solnpr,solnsimul=solnsimul,$
	denspr=denspr,densxpr=densxpr,$
	bsolnup=bsolnup,bsolndn=bsolndn,bsolnpr=bsolnpr,bsolnsimul=bsolnsimul,$
	bdenspr=bdenspr,bdensxpr=bdensxpr,$
	_extra=e
;+
;function	mcmc_chain
;	a simple MCMC chain that generates new parameter values
;	and adds them on to a chain, and returns the chains in
;	a structure of the form
;	{PROB: likelihood, PAR1: parameter 1, ..., PARn: parameter N}
;	where each parameter is stored in a structure of the form
;	{TRACE,MEAN,MEDIAN,MODE,MIN,MAX,STDDEV,HIP68,CLEV68,CLEV90}
;	This is based on the XSPEC "chain" command, which carries
;	out a Metropolis-Hastings chain using the best-fit values
;	and the covariance matrix sigmas as Gaussian proposal
;	distributions.
;
;syntax
;	parstr=mcmc_chain(x,y,parval,parsig,parinit=parinit,parrng=parrng,$
;	ulim=ulim,/useMH,ysig=ysig,freeze=freeze,funcs=funcs,$
;	type=type,ydc=ydc,arfil=arfil,areff=areff,rmfil=rmfil,ties=ties,xfilter=xfilter,$
;	/perbin,nsim=nsim,nburn=nburn,onepar=onepar,verbose=verbose,$
;	bgx=bgx,bgy=bgy,bparval=bparval,bparsig=bparsig,bparini=bparini,$
;	bparrng=bparrng,bysig=bysig,bkgscal=bkgscal,bfreeze=bfreeze,$
;	btype=btype,$
;	mcmc_trace=mcmc_trace,mcmc_prb=mcmc_prb,$
;	mcmc_flux=mcmc_flux,mcmc_counts=mcmc_counts,$
;	mcmc_btrace=mcmc_btrace,mcmc_bprb=mcmc_bprb,$
;	mcmc_bflux=mcmc_bflux,mcmc_bcounts=mcmc_bcounts,$
;	solnup=solnup,solndn=solndn,solnpr=solnpr,solnsimul=solnsimul,$
;	denspr=denspr,densxpr=densxpr,$
;	bsolnup=bsolnup,bsolndn=bsolndn,bsolnpr=bsolnpr,bsolnsimul=bsolnsimul,$
;	bdenspr=bdenspr,bdensxpr=bdensxpr,$
;	...)
;
;parameters
;	x	[INPUT; required] the abscissa array defining the data
;		* if an RMF is given, then this _must_ be in channels
;	y	[INPUT; required] the data array
;		* size of X and Y must match
;		* counts preferred
;	parval	[INPUT; required] array of starting parameter values
;		* if USEMH is set, these are assumed to be the
;		  best-fit values
;	parsig	[INPUT; required] stddev of the Gaussian distribution to
;		use as the proposal distribution.  typically, use values
;		determined with the covariance matrix, inflated by 2-10x
;		as needed.
;		* use of functions other than Gaussians is possible,
;		  but not yet implemented
;		* size must match that of PARVAL, and all elements
;		  must be +ve.  if not,
;		  - if scalar, assumed to be the same error on all PARVAL:
;		    	if 0<PARSIG<1, then set to PARSIG*PARVAL
;		    	if 1<PARSIG<100, then set to (PARSIG/100)*PARVAL
;		    	if PARSIG>100, then set to PARSIG
;		    	if PARSIG<0, then set to abs(PARSIG)
;		  - if vector, then each element controls the width of
;		    the distribution for that parameter.  if any of the
;		    elements are -ve, then the absolute value is treated
;		    as the fractional error on that parameter, same as
;		    the scalar case above.
;
;keywords
;	parinit	[INPUT] use these to set the starting values of the
;		chain.  if given and has the right size, will always
;		override PARVAL.  useful when USEMH is set.
;	parrng	[INPUT] 2-D array of size [2,N(PARVAL)] to set the
;		limits on how far the chain can range.  points that
;		fall outside the range will be forced on to the boundary.
;		* if not given, then no bounds are set
;		* set individual elements to !VALUES.F_NAN to selectively
;		  toggle the bounds
;		* if scalar, then the same fractional range is applied to
;		  all the parameters, e.g., setting "parrng=10" causes the
;		  range for all parameters to go from PARVAL/10 to PARVAL*10
;		  - PARRNG=1 is the same as PARRNG=1e5
;		  - if -ve, then PARVAL+-abs(PARRNG) is used as the range
;		* if 2-element vector, it is the same as the scalar, except
;		  that the range can be asymmetrical, with the 1st element
;		  defining the lower range and the 2nd element defining the
;		  upper range.
;	ulim	[INPUT] array of indices specifying which of the input Y
;		are upper limits (1: UL, 0: not).
;		* An upper limit is treated differently in the calculation of the
;		  likelihood: any parameter value that predicts a model value less
;		  than Y will have no effect on the computed likelihood, but if it
;		  predicts a model value greater than Y, the likelihood becomes 0
;		  (or the chi^2 blows up)
;	useMH	[INPUT] if set, uses a Metropolis-Hastings scheme, with
;		a proposal distribution centered at PARVAL and with
;		width PARSIG.
;		* if not set, uses a Metropolis scheme, with a proposal
;		  distribution that is centered on the current parameter
;		  value and with width PARSIG
;	ysig	[INPUT] the error on Y
;		* default is to use sqrt(abs(Y)+0.75)+1, corresponding
;		  Gehrels' approximation
;	freeze	[INPUT] array of indices of parameters that must be frozen
;		(note: indexing begins from zero!)
;	funcs	[INPUT] name of user-defined procedure that takes as
;		input X and A, and returns YMODEL(X;A).  Any procedure
;		that was written to work with CURVEFIT or GHRS' WFIT or
;		Craig Marwardt's MPFIT will do
;		* default is set to LIBMODEL
;		* why is it a procedure and not a function?
;		  well, ask the person who wrote CURVEFIT.
;	type	[INPUT] string array of models to pass down to LIBMODEL
;	ydc	[INPUT] an array given as an offset to Y(X)
;		* if scalar, assumed to be a constant offset at all X
;		* if vector, size must match that of Y, else is ignored
;	arfil	[INPUT] ancillary response file
;	areff	[INPUT] an array of effective areas may be input in lieu
;		of ARFIL, on the same grid as X
;		* ignored if size does not match X or if ARFIL is given
;	rmfil	[INPUT] response matrix file
;	ties	[INPUT] constraints on parameters, passed w/o comment to ADJUSTIE
;	xfilter	[INPUT] array of X indices to use as a filter
;		* if not set, the entire array is used
;	perbin	[INPUT] if set, computes the model in units of /bin instead
;		of /unit (e.g., /Ang, /keV, etc.)
;	nsim	[INPUT] number of simulations to keep after burn-in
;		* default is 1000L
;	nburn	[INPUT] minimum number of simulations to use for burn-in
;		* default is 0.1*NSIM
;		* some rudimentary checks are made to see whether burn-in
;		  has succeeded, and if it hasn't, then burn-in will
;		  continue.  these simulations are not saved.
;		  *** NOT YET IMPLEMENTED ***
;		* to override the fuzziness of the burn-in, set NBURN to
;		  a negative number, e.g., -100, and then the burn-in phase
;		  will end regardless after abs(NBURN) iterations.
;	onepar	[INPUT] if set, kicks one parameter at a time
;		* the default is to get new deviates for all the
;		  thawed parameters at once
;		* If input as an integer vector, this can also be used to
;		  vary bunches of parameters at once.  the first element
;		  contains the number of bunches defined, NBUNCH.  The
;		  succeeding elements contain, in sequence, the number of
;		  parameters in a given bunch, and their indices.  e.g.,
;		  in the case of a DEM with 81 temperature bins and 30
;		  abundance bins where each set is varied separate bunches,
;		  set this to
;		  	[2, 81, lindgen(81), 30, lindgen(30)+81]
;		  Parameters not included in this list are varied one by
;		  one after all those listed here are dealt with.
;	verbose	[INPUT] controls chatter
;
;	bgx	[INPUT] the abscissa array to define the background
;		* if not given, but BGY is given, assumed to be same as X
;	bgy	[INPUT] the background array
;		* if not given, it is assumed that there is no background
;		* size must match BGX
;		* counts preferred
;		* WARNING: separate ARFs/RMFs for background are not
;		  yet implemented
;	bparval	[INPUT] like PARVAL, for the background model
;		* if not given, but BGY is given, then the background
;		  is not fit, but rather the source model is fit to
;		  background-subtracted data
;		* SOMEDAY: convert to using background marginalization
;	bparsig	[INPUT] like PARSIG, for the background model
;		* ignored if BPARVAL is not defined
;		* size must match that of BPARVAL
;		* if not given, but BPARVAL is given, then the background
;		  model is frozen at BPARVAL (or BPARINI, if given)
;	bparini	[INPUT] like PARINIT, for the background model
;		* ignored if BPARVAL is not defined
;		* size must match that of BPARVAL
;		* if given, overrides BPARVAL
;	bparrng	[INPUT] like PARRNG, for the background model
;		* ignored if BPARVAL is not defined
;	bysig	[INPUT] the error on BGY
;		* default is to use sqrt(abs(BGY)+0.75)+1, corresponding
;		  Gehrels' approximation
;	bkgscal	[INPUT] the background scaling factor for normalization
;		corrections, is the ratio of the geometric area of the
;		background to source, but could also include exposure
;		times, effective area, etc.
;		* correction applied post RMF convolution
;		* if vector, must match the size of BGX
;	bfreeze	[INPUT] array of indices of background parameters that
;		must be frozen
;		* note: indexing starts with 0
;	btype	[INPUT] string array of background models to pass down
;		to LIBMODEL
;	bties	[INPUT] constraints on bakground model parameters, passed
;		w/o comment to ADJUSTIE
;	barfil	[INPUT] ARFIL for background data (NOT IMPLEMENTED)
;	brmfil	[INPUT] RMFIL for background data (NOT IMPLEMENTED)
;	mcmc_trace	[OUTPUT] trace of source parameters
;	mcmc_btrace	[OUTPUT] trace of background parameters
;	mcmc_flux	[OUTPUT] trace of source model fluxes
;	mcmc_bflux	[OUTPUT] trace of background model fluxes
;	mcmc_counts	[OUTPUT] trace of source model counts
;	mcmc_bcounts	[OUTPUT] trace of background model counts
;	mcmc_prb	[OUTPUT] trace of samples from source posterior pdf
;	mcmc_bprb	[OUTPUT] trace of samples from background posterior pdf
;	solnup	[INPUT] a previously obtained set of MCMC solutions that
;		serve as an upper limit to the current calculation
;		* should be an array of size [N(PARVAL),N(SIMULATIONS)]
;	solndn	[INPUT] as SOLNUP, but as a lower limit to the current calc
;		* should be an array of size [N(PARVAL),N(SIMULATIONS)]
;	solnpr	[INPUT] a previously obtained set of MCMC solutions that
;		can be used as a direct prior for the current calculation
;		* should be an array of size [N(PARVAL),N(SIMULATIONS)]
;	solnsimul	[INPUT] a previously obtained set of MCMC solutions
;			that are incorporated into the current MCMC simultaneously
;			(not implemented)
;			* should be an array of size [N(PARVAL),N(SIMULATIONS)]
;	denspr	[INPUT] a formal probability density to act as a prior on the
;			parameters (not implemented)
;	densxpr	[INPUT] the absissa grid for DENSPR (not implemented)
;	bsolnup	[INPUT] as SOLNUP, for the background parameters (not implemented)
;		* should be an array of size [N(BPARVAL),N(SIMULATIONS)]
;	bsolndn	[INPUT] as SOLNDN, for the background parameters (not implemented)
;		* should be an array of size [N(BPARVAL),N(SIMULATIONS)]
;	bsolnpr	[INPUT] as SOLNPR, for the background parameters (not implemented)
;		* should be an array of size [N(BPARVAL),N(SIMULATIONS)]
;	bsolnsimul	[INPUT] as SOLNSIMUL, for the background parameters
;			(not impelmented)
;			* should be an array of size [N(BPARVAL),N(SIMULATIONS)]
;	bdenspr	[INPUT] as DENSPR, for the background parameters (not implemented)
;	bdensxpr	[INPUT] as DENSXPR, for the background parameters
;			(not implemented)
;	_extra	[INPUT ONLY] pass defined keywords to subroutines
;		LIBMODEL: NORM,FWHM,etc.
;		ADJUSTIE: VNAME
;		LIKELI: SOFTLIM,CASH,CASTOR,BINOM
;
;subroutines
;	RDARF()
;	RD_OGIP_RMF()
;	BINERSP()
;	ADJUSTIE
;	KILROY
;
;history
;	vinay kashyap (Apr2008)
;	various updates, bug fixes, and enhancements (VK; Sep2008)
;	changed behavior of vector PARSIG and BPARSIG to match scalar, so now
;	  -ve inputs are treated as absolute sigma; bugfix when chain was
;	getting stuck, the proposal was turning into a delta fn (VK; Oct2009)
;-

;	usage
ok='ok' & np=n_params() & nx=n_elements(x) & ny=n_elements(y)
npar=n_elements(parval) & mpar=n_elements(parsig)
if np eq 0 then ok='Insufficient parameters' else $
 if nx eq 0 then ok='X is not defined' else $
  if ny eq 0 then ok='Y is not defined' else $
   if npar eq 0 then ok='PARVAL is not defined' else $
    if mpar eq 0 then ok='PARSIG is not defined' else $
     if ny ne nx then ok='Y(X) is incompatible with X'
if ok ne 'ok' then begin
  print,'Usage: parstr=mcmc_chain(x,y,parval,parsig,parinit=parinit,parrng=parrng,$'
  print,'       ulim=ulim,/useMH,ysig=ysig,freeze=freeze,funcs=funcs,$'
  print,'       type=type,ydc=ydc,arfil=arfil,areff=areff,rmfil=rmfil,ties=ties,xfilter=xfilter,$'
  print,'       /perbin,nsim=nsim,nburn=nburn,onepar=onepar,verbose=verbose,$'
  print,'       bgx=bgx,bgy=bgy,bparval=bparval,bparsig=bparsig,bparini=bparini,$'
  print,'       bparrng=bparrng,bysig=bysig,bkgscal=bkgscal,bfreeze=bfreeze,$'
  print,'       btype=btype,$'
  print,'       mcmc_trace=mcmc_trace,mcmc_btrace=mcmc_btrace,$'
  print,'       mcmc_flux=mcmc_flux,mcmc_bflux=mcmc_bflux,$'
  print,'       mcmc_counts=mcmc_counts,mcmc_bcounts=mcmc_bcounts,$'
  print,'       mcmc_prb=mcmc_prb,mcmc_bprb=mcmc_bprb,$'
  print,'       solnup=solnup,solndn=solndn,solnpr=solnpr,solnsimul=solnsimul,$'
  print,'       denspr=denspr,densxpr=densxpr,$'
  print,'       bsolnup=bsolnup,bsolndn=bsolndn,bsolnpr=bsolnpr,bsolnsimul=bsolnsimul,$'
  print,'       bdenspr=bdenspr,bdensxpr=bdensxpr,$'
  print,'       ...)'
  if np ne 0 then message,ok,/informational
  return,-1L
endif

tstart=systime()

;	the background
nbx=n_elements(bgx) & nby=n_elements(bgy) & usebg='none'
if nby gt 0 then begin
  if nbx eq nx or nbx eq 0 then begin
    if nby eq ny then usebg='yes'
  endif
endif

;	check keywords
vv=0L & if keyword_set(verbose) then vv=long(verbose[0])>1L
;
npar0=n_elements(parinit) & apar0=parval & if npar0 eq npar then apar0=parinit
;
sigy=sqrt(abs(Y)+0.75)+1. & nye=n_elements(ysig)
if nye gt 0 then begin
  if nye eq ny then sig=ysig else begin
    if nye eq 1 then begin
      ye=ysig[0]
      if ye lt 0 or ye gt 100 then sigy[*]=abs(ye)
      if ye gt 1 and ye le 100 then sigy[*]=y*(ye/100.)
      if ye gt 0 and ye le 1 then sigy[*]=y*ye
    endif
  endelse
endif
;
ipar=lonarr(npar) & nf=n_elements(freeze)
if nf gt 0 then if freeze[0] ne -1 then ipar[freeze]=-1
okpar=where(ipar ge 0,mokpar)
if mokpar eq 0 then begin
  message,'all parameters frozen',/informational
  tmp=create_struct('PAR0',-1L) & for i=1L,npar-1L do tmp=create_struct(tmp,'PAR'+strtrim(i,2),-1L) & return,tmp
endif
;
numsim=1000L & if keyword_set(NSIM) then numsim=long(abs(nsim[0]))>1
numburn=long(0.1*numsim)>1
burnlimit=0	;burn limit is fuzzy
if keyword_set(NBURN) then begin
  numburn=long(abs(nburn[0]))>1
  if NBURN[0] lt 0 then burnlimit=1	;force a hard burn-in limit
endif
;
apar=parval & asig=0.*apar+parsig[0]
if mpar ne npar and mpar ne 1 then begin
  message,'PARSIG not kosher; quitting',/informational
  tmp=create_struct('PAR0',-1L) & for i=1L,npar-1L do tmp=create_struct(tmp,'PAR'+strtrim(i,2),-1L) & return,tmp
endif
if mpar eq 1 then begin
  if parsig[0] lt 0 or parsig[0] gt 100 then asig[*]=abs(parsig[0])
  if parsig[0] gt 1 and parsig[0] le 100 then asig[*]=apar*(parsig[0]/100.)
  if parsig[0] gt 0 and parsig[0] le 1 then asig[*]=apar*parsig[0]
endif else begin
  ;o0=where(parsig lt 0,mo0,complement=op,ncomplement=mop)
  ;if mop gt 0 then asig[op]=parsig[op]
  ;for j=0L,mo0-1L do begin
  for j=0L,mpar-1L do begin
    ;k=o0[j] & xsig=parsig[k]
    k=j & xsig=parsig[k]
    if xsig lt 0 or xsig gt 100 then asig[k]=abs(xsig)
    if xsig gt 1 and xsig le 100 then asig[k]=apar[k]*(xsig/100.)
    if xsig gt 0 and xsig le 1 then asig[k]=apar[k]*xsig
    print,'k,parsig[k],asig[k]',k,xsig,asig[k]
  endfor
endelse
;
nbunch=1L & bunchstr=create_struct('BUNCH0',okpar)	;default is to change all at once
if keyword_set(onepar) then begin	;(ONEPAR
  jpar=ipar+1 & n1par=n_elements(onepar) & i1par=0L & n1bunch=onepar[i1par] & k1par=1L
  while i1par+1L lt n1par do begin	;{parse ONEPAR
    k1par=k1par+1L
    j1par=onepar[i1par+1] & jjpar=onepar[i1par+2:i1par+1+j1par] & jpar[jjpar]=jpar[jjpar]*k1par
    i1par=i1par+1+j1par
  endwhile				;parse ONEPAR}
  o1=where(jpar eq 1,mo1) & nbunch=n1bunch+mo1 & tmp=0
  for i=0L,n1bunch-1L do begin		;{make structure to contain all bunched parameters
    oo=where(jpar eq i+2,moo)
    if i eq 0 then tmp=create_struct('BUNCH0',oo) else tmp=create_struct(tmp,'BUNCH'+strtrim(i,2),oo)
  endfor				;I=0,N1BUNCH-1}
  if n_elements(onepar) eq 1 then begin
    ;	if N1BUNCH=0 then TMP will have been undefined, i.e.,
    ;	every parameter is in a bunch by itself
    tmp=create_struct('BUNCH0',okpar[0])
    for i=1,mo1-1L do tmp=create_struct(tmp,'BUNCH'+strtrim(i,2),o1[i])
  endif else begin
    for i=0L,mo1-1L do tmp=create_struct(tmp,'BUNCH'+strtrim(n1bunch+1+i,2),o1[i])
  endelse
  if n_tags(tmp) gt 0 then bunchstr=tmp
endif					;ONEPAR)
numbunch=n_tags(bunchstr)
;
nrng=n_elements(parrng)
parmin=apar & parmax=parmin	;first freeze everything at the best-fit value
parmin[okpar]=-!values.F_INFINITY & parmax[okpar]=!values.F_INFINITY	;by default, every unfrozen parameter is unbounded
if nrng gt 0 then begin
  if nrng eq 1 then begin
    ;	scalar, means same relative range for all parameters
    xfac=abs(parrng[0]) & if xfac eq 1 then xfac=1e5 & if xfac lt 1 then xfac=1./xfac
    if parrng[0] lt 0 then parmin[okpar]=apar[okpar]-abs(parrng[0]) else parmin[okpar]=apar[okpar]/xfac
    if parrng[0] lt 0 then parmax[okpar]=apar[okpar]+abs(parrng[0]) else parmin[okpar]=apar[okpar]*xfac
  endif
  if nrng eq 2 then begin
    ;	2-element vector, means two-sided range common to all parameters
    xfac1=abs(parrng[0]) & xfac2=abs(parrng[1])
    if xfac1 eq 1 then xfac1=1e5 & if xfac1 lt 1 then xfac1=1./xfac1
    if xfac2 eq 1 then xfac2=1e5 & if xfac2 lt 1 then xfac2=1./xfac2
    if parrng[0] lt 0 then parmin[okpar]=apar[okpar]-abs(parrng[0]) else parmin[okpar]=apar[okpar]/xfac1
    if parrng[0] lt 0 then parmax[okpar]=apar[okpar]+abs(parrng[1]) else parmin[okpar]=apar[okpar]*xfac2
  endif
  if nrng eq 2*mokpar then begin
    ;	a specific range for each thawed parameter
    szr=size(parrng)
    if szr[1] eq 2 then begin
      for i=0L,mokpar-1L do parmin[okpar[i]]=min(parrng[*,i])
      for i=0L,mokpar-1L do parmax[okpar[i]]=max(parrng[*,i])
    endif else begin
      for i=0L,mokpar-1L do parmin[okpar[i]]=min(parrng[i,*])
      for i=0L,mokpar-1L do parmax[okpar[i]]=max(parrng[i,*])
    endelse
  endif
  if mokpar ne npar and nrng eq 2*npar then begin
    ;	a specific range for each thawed parameter
    szr=size(parrng)
    if szr[1] eq 2 then begin
      for i=0L,npar-1L do parmin[i]=min(parrng[*,i])
      for i=0L,npar-1L do parmax[i]=max(parrng[*,i])
    endif else begin
      for i=0L,npar-1L do parmin[i]=min(parrng[i,*])
      for i=0L,npar-1L do parmax[i]=max(parrng[i,*])
    endelse
  endif
endif
;
if not arg_present(funcs) then funcs='libmodel'
;
if n_elements(xfilter) eq 0 then xpick=lindgen(nx) else $
	xpick=[long(xfilter[*])]
;
usesolnpr='no' & mpr=n_elements(solnpr)
if mpr gt 0 then begin
  szpr=size(solnpr)
  if szpr[0] eq 2 then begin
    if szpr[1] eq npar then begin
      msim=szpr[2] & usesolnpr='yes' & prsoln=solnpr
    endif else begin
      if szpr[2] eq npar then begin
	msim=szpr[1] & usesolnpr='yes'
	prsoln=0.*solnpr[0]+fltarr(npar,msim)
	for i=0L,npar-1L do prsoln[i,*]=solnpr[*,i]
      endif
    endelse
  endif
endif
;
usesolnup='no' & mup=n_elements(solnup)
if mup gt 0 then begin
  szup=size(solnup)
  if szup[0] eq 2 then begin
    if szup[1] eq npar then begin
      msim=szup[2] & usesolnup='yes' & upsoln=solnup
    endif else begin
      if szup[2] eq npar then begin
	msim=szup[1] & usesolnup='yes'
	upsoln=0.*solnup[0]+fltarr(npar,msim)
	for i=0L,npar-1L do upsoln[i,*]=solnup[*,i]
      endif
    endelse
  endif
endif
;
usesolndn='no' & mdn=n_elements(solndn)
if mdn gt 0 then begin
  szdn=size(solndn)
  if szdn[0] eq 2 then begin
    if szdn[1] eq npar then begin
      msim=szdn[2] & usesolndn='yes' & dnsoln=solndn
    endif else begin
      if szdn[2] eq npar then begin
	msim=szdn[1] & usesolndn='yes'
	dnsoln=0.*solndn[0]+fltarr(npar,msim)
	for i=0L,npar-1L do dnsoln[i,*]=solndn[*,i]
      endif
    endelse
  endif
endif
;
usesolnsimul='no' & msimul=n_elements(solnsimul)
if msimul gt 0 then begin
  szsimul=size(solnsimul)
  if szsimul[0] eq 2 then begin
    if szsimul[1] eq npar then begin
      msim=szsimul[2] & usesolnsimul='yes' & simulsoln=solnsimul
    endif else begin
      if szsimul[2] eq npar then begin
	msim=szsimul[1] & usesolnsimul='yes'
	simulsoln=0.*solnsimul[0]+fltarr(npar,msim)
	for i=0L,npar-1L do simulsoln[i,*]=solnsimul[*,i]
      endif
    endelse
  endif
endif
if usesolnsimul eq 'yes' then message,'SOLNSIMUL not implemented',/informational
;
usedenspr='no' & ok='ok'
mdpr=n_elements(denspr) & mdxpr=n_elements(densxpr)
szdpr=size(denspr) & szdxpr=size(densxpr)
if mdpr eq 0 then ok='DENSPR not defined' else $
 if mdxpr eq 0 then ok='DENSXPR not defined' else $
  if mdpr lt npar then ok='DENSPR incomplete' else $
   if mdpr ne mdxpr then ok='DENSPR and DENSXPR are incompatible' else $
    if szdpr[0] ne 2 then ok='DENSPR must be 2D' else $
     if szdxpr[0] ne 2 then ok='DENSXPR must be 2D' else $
      if szdpr[1] ne szdxpr[1] then ok='DENSPR does not match DENSXPR' else $
       if szdpr[2] ne szdxpr[2] then ok="DENSPR doesn't match DENSXPR" else $
	if szdpr[1] ne npar and szdpr[2] ne npar then ok='DENSPR not the right size'
if ok eq 'ok' then begin
  message,'p(theta) not implemented',/informational
endif

;	solnup=solnup,solndn=solndn,solnpr=solnpr,solnsimul=solnsimul,$
;	denspr=denspr,densxpr=densxpr,$
;	bsolnup=bsolnup,bsolndn=bsolndn,bsolnpr=bsolnpr,bsolnsimul=bsolnsimul,$
;	bdenspr=bdenspr,bdensxpr=bdensxpr,$

;	the data
xx=x & yy=y
xmin=min(xx,max=xmax)

yoffset=0.*yy & ndc=n_elements(ydc)
if ndc eq 1 then yoffset[*]=ydc[0]
if ndc eq ny then yoffset=ydc

;	background keywords
mokbpar=0L
if usebg ne 'none' then begin	;(background data supplied

  ybg=bgy & usebg='fit'
  if nbx ne nx then xbg=xx else xbg=bgx

  bsigy=sqrt(abs(ybg)+0.75)+1. & nbye=n_elements(bysig)
  if nbye gt 0 then begin
    if nbye eq nby then bsigy=bysig else begin
      if nbye eq 1 then begin
        bgye=bysig[0]
        if bgye lt 0 or bgye gt 100 then bsigy[*]=abs(bgye)
        if bgye gt 1 and bgye le 100 then bsigy[*]=ybg*(bgye/100.)
        if bgye gt 0 and bgye le 1 then bsigy[*]=ybg*bgye
      endif
    endelse
  endif

  nbscal=n_elements(bkgscal)	;BKGSCAL
  if nbscal eq 0 then bgscal=1. else begin	;(BKGSCAL is set
    if nbscal eq nbx then bgscal=bkgscal else bgscal=bkgscal[0]
  endelse					;BKGSCAL is set)

  nbpar=n_elements(bparval) & mokbpar=0L
  if nbpar gt 0 then begin	;(background model parameters are given

    abpar=bparval
    ibpar=lonarr(nbpar) & nfb=n_elements(bfreeze)
    if nfb gt 0 then if bfreeze[0] ne -1 then ibpar[bfreeze]=-1
    okbpar=where(ibpar ge 0,mokbpar)

    absig=fltarr(nbpar) & mbpar=n_elements(bparsig)
    if mbpar eq 0 or mbpar ne nbpar then begin	;(all bkg model parameters are frozen
      ibpar[*]=-1 & okbpar=where(ibpar ge 0,mokbpar)
    endif					;MBPAR=0||MBPAR.ne.NBPAR)

    if mokbpar eq 0 then begin	;(all bkg model parameters are frozen
      message,'all background parameters frozen',/informational
    endif else begin		;MOKBPAR=0)(MOKBPAR.ne.0
      absig=0.*abpar+bparsig[0]
      if mbpar eq 1 then begin
	if bparsig[0] lt 0 or bparsig[0] gt 100 then absig[*]=abs(bparsig[0])
	if bparsig[0] gt 1 and bparsig[0] le 100 then absig[*]=abpar*(bparsig[0]/100.)
	if bparsig[0] gt 0 and bparsig[0] le 1 then absig[*]=abpar*bparsig[0]
      endif else begin
	;o0b=where(bparsig lt 0,mo0b,complement=opb,ncomplement=mopb)
	;if mopb gt 0 then absig[opb]=bparsig[opb]
	;for j=0L,mo0b-1L do begin
	for j=0L,mbpar-1L do begin
	  ;k=o0b[j] & xsig=bparsig[k]
	  k=j & xsig=bparsig[k]
	  if xsig lt 0 or xsig gt 100 then absig[k]=abs(xsig)
	  if xsig gt 1 and xsig le 100 then absig[k]=abpar[k]*(xsig/100.)
	  if xsig gt 1 and xsig le 1 then absig[k]=abpar[k]*xsig
	endfor
      endelse

      nrngb=n_elements(bparrng)
      bparmin=abpar & bparmax=bparmin	;first freeze everything at the best-fit value
      bparmin[okbpar]=-!values.F_INFINITY & bparmax[okbpar]=!values.F_INFINITY	;by default, every unfrozen parameter is unbounded
      if nrngb gt 0 then begin
        if nrngb eq 1 then begin
          ;	scalar, means same relative range for all parameters
          xfac=abs(bparrng[0]) & if xfac eq 1 then xfac=1e5 & if xfac lt 1 then xfac=1./xfac
          if bparrng[0] lt 0 then bparmin[okbpar]=abpar[okbpar]-abs(bparrng[0]) else bparmin[okbpar]=abpar[okbpar]/xfac
          if bparrng[0] lt 0 then bparmax[okbpar]=abpar[okbpar]+abs(bparrng[0]) else bparmin[okbpar]=abpar[okbpar]*xfac
        endif
        if nrngb eq 2 then begin
          ;	2-element vector, means two-sided range common to all parameters
          xfac1=abs(bparrng[0]) & xfac2=abs(bparrng[1])
          if xfac1 eq 1 then xfac1=1e5 & if xfac1 lt 1 then xfac1=1./xfac1
          if xfac2 eq 1 then xfac2=1e5 & if xfac2 lt 1 then xfac2=1./xfac2
          if bparrng[0] lt 0 then bparmin[okbpar]=abpar[okbpar]-abs(bparrng[0]) else bparmin[okbpar]=abpar[okbpar]/xfac1
          if bparrng[0] lt 0 then bparmax[okbpar]=abpar[okbpar]+abs(bparrng[1]) else bparmin[okbpar]=abpar[okbpar]*xfac2
        endif
        if nrngb eq 2*mokbpar then begin
          ;	a specific range for each thawed parameter
          szrb=size(bparrng)
          if szrb[1] eq 2 then begin
            for i=0L,mokbpar-1L do bparmin[okbpar[i]]=min(bparrng[*,i])
            for i=0L,mokbpar-1L do bparmax[okbpar[i]]=max(bparrng[*,i])
          endif else begin
            for i=0L,mokbpar-1L do bparmin[okbpar[i]]=min(bparrng[i,*])
            for i=0L,mokbpar-1L do bparmax[okbpar[i]]=max(bparrng[i,*])
          endelse
        endif
        if mokbpar ne nbpar and nrngb eq 2*nbpar then begin
          ;	a specific range for each thawed parameter
          szrb=size(bparrng)
          if szrb[1] eq 2 then begin
            for i=0L,nbpar-1L do bparmin[i]=min(bparrng[*,i])
            for i=0L,nbpar-1L do bparmax[i]=max(bparrng[*,i])
          endif else begin
            for i=0L,nbpar-1L do bparmin[i]=min(bparrng[i,*])
            for i=0L,nbpar-1L do bparmax[i]=max(bparrng[i,*])
          endelse
        endif
      endif

    endelse			;MOKBPAR.ne.0)

    nbpar0=n_elements(bparinit) & abpar0=bparval & if nbpar0 eq nbpar then abpar0=bparinit
    abpar=bparval

  endif else begin		;NBPAR>0)(subtract background instead
    usebg='subtract'
    yy=yy-ybg/bgscal
    sigy=sqrt(sigy^2+bsigy^2/bgscal^2)
  endelse			;NBPAR=0)

  ybfunc0=0.*yy
endif				;USEBG.ne.'NONE')

;	read in ARF if present
xarf=xx & effar=fltarr(nx)+1.0 & dxarf=mid2bound(xarf, _extra=e)
if keyword_set(areff) then begin
  nea=n_elements(areff)
  if nea eq nx then effar=areff
endif
ARFinterp=0	;no need to interpolate to fit the input grid (applies iff RMFMULT=0, see below)
if arg_present(arfil) then begin
  ok='ok'
  if n_elements(arfil) gt 1 then ok='ARFIL must be a scalar' else $
   if size(arfil,/type) ne 7 then ok='ARFIL must be a filename' else $
    arf=findfile(arfil,count=narf) & if narf eq 0 then ok='ARFIL: file not found'
  if ok ne 'ok' then begin
    message,ok,/informational
    message,'using flat ARF',/informational
  endif else begin
    effar=rdarf(arfil[0],arfstr)
    xarf=0.5*(arfstr.ELO+arfstr.EHI)
    dxarf=arfstr.EHI-arfstr.ELO
  endelse
  ;	the model is always calculated on the ARF grid.  is it necessary
  ;	to interpolate it on to the input x-grid?
  if n_elements(xarf) ne nx then ARFinterp=1 else $
   if median(arfstr.EHI-arfstr.ELO) ne median(xx[1:*]-xx) then ARFinterp=1 else $
    if total(xarf-xx) gt min(arfstr.EHI-arfstr.ELO) then ARFinterp=1
endif
narf=n_elements(xarf)

;	read in RMF if present
RMFmult=0
if arg_present(rmfil) then begin
  ok='ok'
  if n_elements(rmfil) gt 1 then ok='RMFIL must be a scalar' else $
   if size(rmfil,/type) ne 7 then ok='RMFIL must be a filename' else $
    rmf=findfile(rmfil,count=nrmf) & if nrmf eq 0 then ok='RMFIL: file not found'
  if ok ne 'ok' then begin
    message,ok,/informational
    message,'using delta-function RMF',/informational
    RMFmult=0
  endif else begin
    RMFmult=1
    rmstr=rd_ogip_rmf(rmfil[0])
    ie=binersp(xarf,rmstr.ELO,rmstr.EHI)
    ;	recast the RMF to speed up the model calculation
    nchan=rmstr.NCHAN & n_grp=rmstr.N_GRP & f_chan=rmstr.F_CHAN & n_chan=rmstr.N_CHAN & matrix=rmstr.MATRIX & firstchan=rmstr.FIRSTCHAN
    szg=size(f_chan)
    rmatrix=fltarr(nx,narf)
    for i=0L,nx-1L do begin
      j=ie[i] & tmp=fltarr(nchan)
      for ig=0L,n_grp[j]-1 do begin
	if szg[0] eq 1 then begin
	  i0=f_chan[j]-firstchan & i1=n_chan[j] & tmp[i0:i0+i1-1]=matrix[0:i1-1,j]
	endif else begin
	  i0=f_chan[ig,j]-firstchan & i1=n_chan[ig,j] & tmp[i0:i0+i1-1]=tmp[i0:i0+i1-1]+matrix[0:i1-1,j]
	endelse
      endfor
      ;rmatrix[*,j]=tmp[long(xmin):long(xmax)]
      rmatrix[*,j]=tmp
    endfor
  endelse
endif

;	now run the MCMC chain
mcmc_trace=fltarr(mokpar,numsim+1L)
mcmc_prb=fltarr(numsim+1L)
mcmc_burn=fltarr(mokpar,numburn)
mcmc_flux=fltarr(numsim+1L)
mcmc_counts=fltarr(numsim+1L)
ybest=0*yy & bestpar=apar
if usebg eq 'fit' and mokbpar gt 0 then begin
  mcmc_btrace=fltarr(mokbpar,numsim+1L)
  mcmc_bprb=fltarr(numsim+1L)
  mcmc_bburn=fltarr(mokbpar,numburn)
  mcmc_bflux=fltarr(numsim+1L)
  mcmc_bcounts=fltarr(numsim+1L)
  ybestbg=0*yy & bestbpar=abpar
endif

if vv gt 10 then begin
  pmulti=!p.multi
  window,0 & window,1 & if mokbpar gt 0 then window,2
endif

go_on=1 & iburn=-1L & isim=-1L & done_burnin=0 & irenormb=-1 & irenorm=-1
kburn=0L
oldprb=0 & oldbprb=0 & yspec0=0 & ybspec0=0
curpar=apar & if usebg eq 'fit' then curbpar=abpar
;
while go_on do begin	;{see chain.  see chain run.  run, chain, run
  if iburn ge numburn then isim=isim+1L else iburn=iburn+1L

  ;{
  ;	ARF and RMF variations to be included here
  ;}

  ;{
  ;	first, one iteration of the background model
  if usebg eq 'fit' then begin		;(use background model
    rdev=randomn(seed,nbpar) & newbpar=curbpar
    if not keyword_set(useMH) then abpar0=curbpar
    jjb=0
    for ii=0L,mokbpar do begin	;{step through each bkg parameter
      if ii eq 0 then begin
	newbpar=abpar0
	;this way the loop will be run at least once and YBSPEC0 will be populated
      endif else begin
        op=okbpar[ii-1L] & newbpar=curbpar & newbpar[op]=abpar0[op]+rdev[op]*absig[op]
        ;	enforce ranges
        newbpar[op]=newbpar[op]>reform(bparmin[op])
        newbpar[op]=newbpar[op]<reform(bparmax[op])
        if arg_present(bties) then adjustie,newbpar,ties=bties, _extra=e
      endelse

      ;	compute background model
      if (mokbpar gt 0) or (not keyword_set(ybspec0)) then begin
	call_procedure,funcs,xarf,newbpar,ybfunc,type=btype, _extra=e

	;	model flux
	bflux=total(ybfunc[xpick])
	if not keyword_set(perbin) then bflux=total(ybfunc[xpick]*dxarf[xpick])

        ;	apply ARF
        ybmod=ybfunc*effar
	if not keyword_set(perbin) then ybmod=ybmod*dxarf

        ;	apply RMF
        if keyword_set(RMFmult) then begin
	  ybspec=fltarr(nx) & for kk=0L,narf-1L do ybspec=ybspec+rmatrix[*,kk]*ybmod[kk]
        endif else begin
	  ybspec=ybmod
	  if keyword_set(ARFinterp) then ybspec=interpol(ybmod,xarf,xx)
        endelse

        ;	compute likelihood
        bprb=likeli(ybg[xpick],ybspec[xpick],dsigma=bsigy[xpick],ulim=ulim,/chi2, _extra=e)
	bprb=bprb/2.	;because LIKELI() returns -2ln(Likelihood)

	;	compute prior and thence posterior
	bprbprior=0.
	bprb=bprb+bprbprior

	;	model counts
	bcounts=total(ybspec[xpick])

	;	force renormalize the first few times
	if irenormb lt mokbpar then begin
	  if bcounts gt 0 then renormb=total(ybg[xpick])/float(bcounts) else renormb=1.0
	  if vv gt 0 then message,'renormalizing background by a factor '+strtrim(renormb,2),/informational
	  call_procedure,funcs,xarf,newbpar,ybfunc,jnk,updatbpar,type=btype,renorm=renormb, _extra=e
	  newbpar=updatbpar & curbpar=updatbpar & abpar0=updatbpar & ybspec=ybspec/renormb
	  irenormb=irenormb+1
	endif

	;	do the MCMC test on whether to keep or discard new parameters
	if keyword_set(oldbprb) then begin		;(OLDBPRB is set
	  accept_fnb = (oldbprb-bprb)<0
	  metro_checkb = alog(randomu(seed))
	  if metro_checkb le accept_fnb then begin	;(accept new value
	    jjb=jjb+1L
	    oldbprb=bprb & curbprb=bprb & curbpar=newbpar
	    ybspec0=ybspec
	    if bprb lt bestbprb then begin
	      ybestbg=ybspec & bestbprb=bprb & bestbpar=curbpar & bestbflux=bflux & bestbct=bcounts
	      if mokbpar gt 0 then begin
	        mcmc_btrace[*,numsim]=bestbpar[okbpar]
	        mcmc_bprb[isim+1:numsim]=bestbprb
	        mcmc_bflux[numsim]=bestbflux
	        mcmc_bcounts[numsim]=bestbct
	      endif
	    endif
	    if ii gt 0 then print,'B:	',op,jjb,curbprb,curbpar[op]
	  endif						;METRO_CHECKB<ACCEPT_FNB)
	endif else begin				;OLDBPRB)(set OLDBPRB
	  oldbprb=bprb & curbprb=bprb & curbpar=newbpar
	  ybspec0=ybspec
          bestbprb=bprb & ybestbg=ybspec & bestbpar=curbpar & bestbflux=bflux & bestbct=bcounts
	endelse						;set OLDPRBB)
      endif

    endfor				;II=0,mokbpar-1}
  endif else begin			;USEBG='fit')(no bkg
    ybspec0=0.*yy & bprb=0.
  endelse				;USEBG='none')
  ;}

  ;{
  ;	next, one iteration of the source model, given the background model

  ;	get new parameter deviates
  rdev=randomn(seed,npar) & newpar=curpar
  if not keyword_set(useMH) then apar0=curpar
  jj=0
  for ii=0L,numbunch-1L do begin

    op=bunchstr.(ii) & newpar=curpar & newpar[op]=apar0[op]+rdev[op]*asig[op]
    mop=n_elements(op)
    ;	enforce ranges
    newpar[op]=newpar[op]>reform(parmin[op])
    newpar[op]=newpar[op]<reform(parmax[op])
    ;	enforce ties
    if arg_present(ties) then adjustie,newpar,ties=ties, _extra=e

    ;	compute model
    call_procedure,funcs,xarf,newpar,yfunc,type=type, _extra=e

    ;	model flux
    flux=total(yfunc[xpick])
    if not keyword_set(perbin) then flux=total(yfunc[xpick]*dxarf[xpick])

    ;	apply ARF
    ymod=yfunc*effar
    if not keyword_set(perbin) then ymod=ymod*dxarf

    ;	apply RMF
    if keyword_set(RMFmult) then begin
      ;yspec=rmatrix#ymod
      yspec=fltarr(nx) & for kk=0L,narf-1L do yspec=yspec+rmatrix[*,kk]*ymod[kk]
    endif else begin
      yspec=ymod
      if keyword_set(ARFinterp) then yspec=interpol(ymod,xarf,xx)
    endelse

    ;	model counts
    counts=total(yspec[xpick])

    ;forcibly renormalize the first few tims
    if irenorm lt mokpar then begin
      if counts gt 0 then renorm=total(yy[xpick])/float(counts) else renorm=1.0
      if vv gt 5 then print,'*',newpar
      if vv gt 0 then message,'renormalizing by a factor '+strtrim(renorm,2),/informational
      call_procedure,funcs,xarf,newpar,yfunc,jnk,updatpar,type=type,renorm=renorm, _extra=e
      if vv gt 5 then print,'**',newpar
      newpar=updatpar & curpar=updatpar & apar0=updatpar & yspec=yspec/renorm
      irenorm=irenorm+1
    endif

    ;	compute likelihood
    yspec=yspec+yoffset
    if usebg ne 'none' then yspec=yspec+ybspec0/bgscal
    prb=likeli(yy[xpick],yspec[xpick],dsigma=sigy[xpick],ulim=ulim,/chi2, _extra=e)
    					; returns -2*alog(p(D|M))
    prb=prb/2.
    print,'PRB:	',prb

    ;	compute prior
    prbprior=0.

    if usesolnup eq 'yes' then begin
      numprob=0. & denomprob=0.
      for j=0L,mop-1L do begin
	tmp=reform(upsoln[op[j],*]) & jnk=where(tmp le curpar[op[j]],mjnk)
	denomprob=denomprob+n_elements(tmp) & numprob=numprob+float(mjnk)
      endfor
      if numprob gt 0 then prbprior=prbprior-alog(numprob)+alog(denomprob) else prbprior=prbprior+30.+alog(denomprob)
    endif
    if usesolndn eq 'yes' then begin
      numprob=0. & denomprob=0.
      for j=0L,mop-1L do begin
	tmp=reform(dnsoln[op[j],*]) & jnk=where(tmp ge curpar[op[j]],mjnk)
	denomprob=denomprob+n_elements(tmp) & numprob=numprob+float(mjnk)
      endfor
      if numprob gt 0 then prbprior=prbprior-alog(numprob)+alog(denomprob) else prbprior=prbprior+30.+alog(denomprob)
    endif
    if usesolnpr eq 'yes' then begin
      ;	this assumes that the probability that a parameter value "belongs" to the
      ;	input distribution is determined as the total number less the difference
      ;	between how many are present above and below, normalized to the total.  i.e.,
      ;	n=10000 & r=randomn(seed,n) & os=sort(r) & r=r[os] & plot,r,lindgen(n),psym=3
      ;	y=(findgen(101)-50.)/10. & p=0.*y
      ;	for i=0,100 do begin & om=where(r le y[i],mom) & op=where(r ge y[i],mop) & $
      ;	p[i]=float((mom+mop)-abs(mom-mop))/float(mom+mop) & endfor
      ;	plot,y,p & yh=histogram(r,min=-5,max=5,binsize=0.1) & oplot,y,yh/float(max(yh)),psym=10
      ;	tyh=total(yh,/cumul) & oplot,y+0.1,tyh/(n/2.),col=2	;similarly other side
      numprob=0. & denomprob=0.
      for j=0L,mop-1L do begin
	tmp=reform(prsoln[op[j],*])
	jnk1=where(tmp ge curpar[op[j]],mjnk1) & jnk2=where(tmp le curpar[op[j]],mjnk2)
	mjnk=float(mjnk1+mjnk2) & dmjnk=abs(mjnk1-mjnk2)
	denomprob=denomprob+(mjnk1+mjnk2) & numprob=numprob+(mjnk-dmjnk)
      endfor
      if numprob gt 0 then prbprior=prbprior-alog(numprob)+alog(denomprob) else prbprior=prbprior+30.+alog(denomprob)
    endif

    ;	and thence posterior
    prb=prb+prbprior

    ;	do the MCMC test on whether to keep or discard new parameters
    if keyword_set(oldprb) then begin
      accept_fn = (oldprb-prb)<0
      metro_check = alog(randomu(seed))
      if metro_check lt accept_fn then begin	;(accept new value
	oldprb=prb & curprb=prb & curpar=newpar
	yspec0=yspec
	if prb lt bestprb then begin
	  ybest=yspec & bestprb=prb & bestpar=curpar & bestflux=flux & bestct=counts
          mcmc_trace[*,numsim]=bestpar[okpar]
          mcmc_prb[isim+1:numsim]=bestprb
          mcmc_flux[numsim]=bestflux
          mcmc_counts[numsim]=bestct
	endif
	jj=jj+1L
	print,'S:	',op,jj,curprb,curpar[op]
      endif else begin				;METRO_CHECK<ACCEPT_FN)(not
        print,'xS:	',op,jj,oldprb,prb,newpar[op]
      endelse					;do not accept)
    endif else begin
      oldprb=prb & curprb=prb & curpar=newpar
      yspec0=yspec
      bestprb=prb & ybest=yspec & bestpar=curpar & bestflux=flux & bestct=counts
    endelse 

  endfor
  ;}

  if jj eq 0 then begin		;(redo iteration if no new parameter value was accepted
    print,'rejected!'
    ;if iburn ge numburn then isim=isim-1L else iburn=iburn-1L
    if iburn le numburn then iburn=iburn-1L
    kburn=kburn+1L
    if kburn gt (((numburn-1L)/5) > 10) then begin
      message,'TOO MANY REJECTS!  Resetting the comparison point',/informational
      print,iburn,isim,oldprb,prb
      oldprb=2*prb	;just a hack to avoid the infinite loop, kick it out to restart again
      kburn=0L
    endif
  endif else begin		;JJ=0)(store if past burn-in

    ;	store if past burn-in
    kburn=0L
    if isim ge 0 then begin
      if mokbpar gt 0 then begin
        mcmc_btrace[*,isim]=curbpar[okbpar]
        mcmc_bprb[isim]=curbprb
        mcmc_bflux[isim]=bflux
        mcmc_bcounts[isim]=bcounts
      endif
      mcmc_trace[*,isim]=curpar[okpar]
      mcmc_prb[isim]=curprb
      mcmc_flux[isim]=flux
      mcmc_counts[isim]=counts
    endif else begin
      if mokbpar gt 0 then mcmc_bburn[*,iburn]=curbpar[okbpar]
      mcmc_burn[*,iburn]=curpar[okpar]
    endelse

    if vv gt 10 then begin

      wset,0 & !p.multi=[0,2,2,0,1]
      if mokbpar gt 0 then print,iburn,isim,bprb,oldbprb,curbpar[okbpar]
      print,iburn,isim,prb,oldprb,curpar[okpar]
      plot,xx[xpick],yy[xpick],psym=10,/xlog,/xs,/ynoz
      oplot,xx[xpick],yspec[xpick],color=2
      oplot,xx[xpick],yspec0[xpick],color=3
      if usebg eq 'fit' then oplot,xx[xpick],ybspec[xpick]/bgscal,color=4
      ;
      plot,xx[xpick],yy[xpick]-yspec[xpick],psym=10,/xlog,/xs,/ynoz
      oplot,xx[xpick],yy[xpick]-yspec0[xpick],color=2,psym=3
      oplot,xx[xpick],0*xpick
      ;
      if usebg eq 'fit' then begin
        plot,xx[xpick],ybg[xpick],psym=10,/xlog,/xs,/ynoz
        oplot,xx[xpick],ybspec[xpick],col=2
        oplot,xx[xpick],ybspec0[xpick],col=3
      ;
        plot,xx[xpick],ybg[xpick]-ybspec[xpick],psym=10,/xlog,/xs,/ynoz
        oplot,xx[xpick],ybg[xpick]-ybspec0[xpick],color=2,psym=3
      endif
 
      wset,1 & !p.multi=[0,2,mokpar/2+1]
      if isim ge 0 then $
        for ip=0L,mokpar-1 do plot,mcmc_trace[okpar[ip],0L:isim],title='PAR'+strtrim(okpar[ip],2) else $
        for ip=0L,mokpar-1 do plot,mcmc_burn[okpar[ip],0L:iburn],title='PAR'+strtrim(okpar[ip],2)
      plot,(mcmc_prb-bestprb)>0,title='PROB';,/ynoz
  
      if mokbpar gt 0 then begin
        wset,2 & !p.multi=[0,2,mokbpar/2+1]
        if isim ge 0 then $
          for ip=0L,mokbpar-1 do plot,mcmc_btrace[okbpar[ip],0L:isim],title='BPAR'+strtrim(okbpar[ip],2) else $
          for ip=0L,mokbpar-1 do plot,mcmc_bburn[okbpar[ip],0L:iburn],title='BPAR'+strtrim(okbpar[ip],2)
        plot,mcmc_bprb-bestbprb,title='BPROB';,/ynoz
      endif
    endif

  endelse			;JJ>0)

  ;	check whether burn-in has converged
  ;	currenly not implemented
  if iburn eq numburn-1L then chain_converged=1 else chain_converged=0
  
  ;	stop burn-in?
  if keyword_set(burnlimit) and iburn eq numburn-1L then done_burnin=1
  if not keyword_set(burnlimit) and iburn eq numburn-1L then begin
    ;	check that the chain has converged
    if keyword_set(chain_converged) then done_burnin=1 else begin
      mcmc_burn=shift(mcmc_burn,-mokpar)
      iburn=iburn-1L
    endelse
  endif
  if keyword_set(done_burnin) then begin
    iburn=numburn+1
    ;	compute sdev(PAR) from the burn and update PARSIG and BPARSIG
    tmpmean=asig & tmpsig=asig
    for jj=0L,mokpar-1L do tmpmean[jj]=mean(mcmc_burn[jj,numburn/2:*])
    for jj=0L,mokpar-1L do tmpsig[jj]=abs(stddev(mcmc_burn[jj,numburn/2:*]/tmpmean[jj])*tmpmean[jj])
    for jj=0L,mokpar-1L do if tmpsig[jj] eq 0. then tmpsig[jj]=asig[jj] ;ignore next step if for some reason the chain gets stuck
    asig=asig < tmpsig
    if mokbpar gt 0 then begin
      tmpmeanb=absig & tmpsigb=absig
      for jj=0L,mokbpar-1L do tmpmeanb[jj]=mean(mcmc_bburn[jj,numburn/2:*])
      for jj=0L,mokbpar-1L do tmpsigb[jj]=abs(stddev(mcmc_bburn[jj,numburn/2:*]/tmpmeanb[jj])*tmpmeanb[jj])
      absig=absig < tmpsigb
    endif
  endif

  ;	stopping rule
  if isim eq numsim-1L then go_on=0	;stop if all iterations are done
endwhile		;GO_ON}

if vv gt 10 then begin
  pmulti=!p.multi
endif

;	format the output
parstr=create_struct('PROB',mcmc_prb)
for i=0L,npar-1L do begin
  tmp=min(abs(i-okpar),ii)
  pmean=-1 & pmed=-1 & pmode=-1 & pmin=-1 & pmax=-1 & psig=-1 & pup=-1 & pdn=-1 & p16=-1 & p84=-1 & p05=-1 & p95=-1
  if tmp eq 0 then begin
    arr=reform(mcmc_trace[ii,*])
    pmean=mean(arr) & pmed=median(arr)
    pmin=min(arr,max=pmax)
    if pmin gt -1e5 and pmax lt 1e5 then psig=stddev(arr) else psig=pmean*stddev(arr/pmean)
    ph=hipd_interval(arr,/fsample,fmode=pmode, _extra=e) & pup=ph[1] & pdn=ph[0]
    pc=hipd_interval(arr,/fsample,clev=0.68) & p16=pc[0] & p84=pc[1]
    pc=hipd_interval(arr,/fsample,clev=0.90) & p05=pc[0] & p95=pc[1]
  endif else arr=-1L
  ststr=create_struct('TRACE',arr,'MEAN',pmean,'MEDIAN',pmed,'MODE',pmode,$
	'MIN',pmin,'MAX',pmax,'STDDEV',psig,$
	'HIP68',[pdn,pup],'CLEV68',[p16,p84],'CLEV90',[p05,p95])
  parstr=create_struct(parstr,'PAR'+strtrim(i,2),ststr)
endfor
parstr=create_struct(parstr,'XFILTER',xpick,'BESTPAR',bestpar,'YBEST',ybest,'FLUX',mcmc_flux,'COUNTS',mcmc_counts)
if usebg eq 'fit' then begin
  bparstr=create_struct('PROB',mcmc_bprb)
  pmean=-1 & pmed=-1 & pmode=-1 & pmin=-1 & pmax=-1 & psig=-1 & pup=-1 & pdn=-1 & p16=-1 & p84=-1 & p05=-1 & p95=-1
  for i=0L,nbpar-1L do begin
    tmp=min(abs(i-okbpar),ii)
    if tmp eq 0 then begin
      arr=reform(mcmc_btrace[ii,*])
      pmean=mean(arr) & pmed=median(arr)
      pmin=min(arr,max=pmax)
      if pmin gt -1e5 and pmax lt 1e5 then psig=stddev(arr) else psig=pmean*stddev(arr/pmean)
      ph=hipd_interval(arr,/fsample,fmode=pmode, _extra=e) & pup=ph[1] & pdn=ph[0]
      pc=hipd_interval(arr,/fsample,clev=0.68) & p16=pc[0] & p84=pc[1]
      pc=hipd_interval(arr,/fsample,clev=0.90) & p05=pc[0] & p95=pc[1]
    endif else arr=-1L
    bststr=create_struct('TRACE',arr,'MEAN',pmean,'MEDIAN',pmed,'MODE',pmode,$
	'MIN',pmin,'MAX',pmax,'STDDEV',psig,$
	'HIP68',[pdn,pup],'CLEV68',[p16,p84],'CLEV90',[p05,p95])
    bparstr=create_struct(bparstr,'PAR'+strtrim(i,2),bststr)
  endfor
  bparstr=create_struct(bparstr,'BESTPAR',bestbpar,'YBEST',ybestbg,'FLUX',mcmc_bflux,'COUNTS',mcmc_bcounts)
endif
if usebg eq 'fit' then parstr=create_struct(parstr,'BKG',bparstr)
if n_tags(e) gt 0 then parstr=create_struct(parstr,'EXTRAPAR',e)	;this is so we can reconstruct the model independently if necessary
tstop=systime() & parstr=create_struct(parstr,'RUNTIME',[tstart,tstop])
parstr=create_struct(parstr,'X',xx,'Y',yy)
if usebg ne 'none' then parstr=$
	create_struct(parstr,'XBG',xbg,'YBG',ybg,'BKGSCAL',bgscal)

return,parstr
end
