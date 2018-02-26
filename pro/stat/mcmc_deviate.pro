function mcmc_deviate,par,sigpar,rngpar=rngpar,thaw=thaw,sampar=sampar,$
	seed=seed,sclpar=sclpar, _extra=e
;+
;function	mcmc_deviate
;	returns a set of parameters as a deviation from the current set
;
;	ideally, this function should be called from MCMC_STEP().
;	unfortunately, while that improves readability, it helps
;	with neither speed nor efficiency of computation.  e.g.,
;	when the parameters are sampled from SAMPAR, using this
;	function will result in a lot of wasted samples of the
;	parameters, and when each parameter is updated singly, a
;	lot of wasted tests.  so, while I will try to keep this
;	as up to date as possible, MCMC_STEP() will retain all the
;	functionality present here, and will likely be the program
;	of choice for the foreseeable future.
;
;syntax
;	testpar=mcmc_deviate(par,sigpar,rngpar=rngpar,thaw=thaw,sampar=sampar,$
;	sclpar=sclpar,seed=seed)
;
;parameters
;	par	[INPUT; required] parameters of the model to fit to Y(X)
;		* parameters should be 1-D array
;		* if 2-D, the size of the 2nd dimension describes the
;		  number of chains that are concurrently running
;	sigpar	[INPUT; required] this is used to determine how to get
;		a new set of parameters with PAR as the starting point
;		* if size is identical to the first dimension of PAR, then
;		  +ve ==> assumed to be the stddev of Gaussian centered on PAR
;		  -ve ==> abs value assumed to be maximum deviation from PAR
;
;keywords
;	rngpar	[INPUT] the allowed range for each parameter as a 2-D
;		array of size (N(PAR),2)
;		* overrides output of SIGPAR
;		* if 1-D array, range assumed to be symmetrical
;		  -- additive if +ve
;		  -- abs value multiplicative if -ve
;		* ignored if size doesn't match
;		* WARNING: RNGPAR[.,0] must be less than RNGPAR[.,1]
;		  in the interests of speed, this is _not_ checked for
;	thaw	[INPUT] array indicating which parameters are frozen (0) and
;		which are thawed (1)
;		* size _must_ match number of parameters, else all
;		  parameters are assumed thawed
;		* overrides SIGPAR and RNGPAR
;	sampar	[INPUT] used to determine how to pick the parameters
;		to vary within the program:
;		-- if not set, then parameters are picked in sequence
;		-- if set to scalar, then parameters are picked randomly
;		-- if set to array of same size as PAR, assumed to be
;		   a sampling distribution for the parameters
;	sclpar	[INPUT] a flag that describes how PAR should be scaled
;		prior to taking a deviate
;		1 == alog(PAR) scaling ;	-1 == exp(PAR) scaling
;		10 == alog10(PAR) scaling ;	-10 == 10^PAR scaling
;		100 == LOGIT() scaling ;	-100 == UNLOGIT() scaling
;		2 == sqrt(PAR) scaling ;	-2 == PAR^2 scaling
;		0, or other == no scaling
;		* NOTE: no checks are made to ensure that values are
;		  defined post-transformation -- it is the responsibility
;		  of the user to ensure that the scaling makes sense
;		* if scalar, the same scaling is applied to all the parameters
;		* if vector, size must match the number of parameters, or
;		  else this is ignored
;	seed	[I/O] seed for random number generator
;	_extra	[JUNK] here only to prevent crashing the program
;
;subroutines
;	LOGIT()
;	UNLOGIT()
;
;history
;	vinay kashyap (Jun2005)
;-

;	usage
ok='ok' & np=n_params()
szp=size(par) & nszp=n_elements(szp) & npar=szp[1] & npe=n_elements(sigpar)
if np lt 2 then ok='Insufficient parameters' else $
 if szp[0] eq 0 then ok='PAR are undefined' else $
  if szp[0] gt 2 then ok='cannot understand PAR > 2-D' else $
   if npe eq 0 then ok='SIGPAR are undefined' else $
    if npe ne npar then ok='PAR and SIGPAR are incompatible'
if ok ne 'ok' then begin
  print,'Usage: testpar=mcmc_deviate(par,sigpar,rngpar=rngpar,thaw=thaw,$'
  print,'       sampar=sampar,seed=seed,sclpar=sclpar)'
  print,'  returns a set of parameters as a deviation from the current set'
  if np ne 0 then message,ok,/informational
  return,-1L
endif

;	check input
  ;
;  number of chains
if szp[0] eq 2 then nchain=szp[2] else nchain=1
  ;
;  range of PAR
szr=size(rngpar) & userng=0 & parrng=reform([par,par],npar,2)
case szr[0] of
  0: ;nothing
  1: begin
    if szr[1] eq npar then begin
      for i=0L,npar-1L do begin
	if rngpar[i] gt 0 then parrng[i,*]=par[i]+rngpar[i]*[-1,1]
	if rngpar[i] lt 0 then parrng[i,*]=par[i]*$
		[1./abs(rngpar[i]),abs(rngpar[i])]
      endfor
    endif else $
      message,'RNGPAR incompatible with PAR; ignoring',/informational
  end
  2: begin
    if szr[1] eq npar and szr[2] eq 2 then parrng=rngpar else $
      message,'RNGPAR incompatible with PAR; ignoring',/informational
  end
  else: message,'RNGPAR incompatible with PAR; ignoring',/informational
endcase
  ;
;  frozen parameters
nthaw=n_elements(thaw) & thawed=intarr(npar)+1
if nthaw eq npar then thawed[*]=thaw[*] else begin
  message,'THAW is incompatible with PAR; thawing all parameters',$
	/informational
  nthaw=npar
endif
opar=where(thaw ne 0,nnpar)
  ;
;  parameter sampling distribution
dsamp=[0.,fltarr(npar)+1./npar] & nsamp=n_elements(sampar)
ksamp = nsamp < npar
if ksamp gt 0 then dsamp[1L:ksamp]=sampar[0L:ksamp-1L]
sampfn=total(abs(dsamp),/cumulative) & sampfn=sampfn/sampfn[npar]
rsamp=reform(randomu(seed,npar*nchain),npar,nchain)
rpar=long(interpol(lindgen(npar+1L),sampfn,rsamp))
  ;
;  parameter scaling transformations
scltyp=intarr(npar) & nscl=n_elements(sclpar)
if nscl eq 1 then scltyp[*]=sclpar[0]
if nscl eq npar then scltyp=sclpar

;	get deviates
for ii=0L,nnpar-1L do begin		;{thawed parameters
  ip=opar[ii]		;look only at parameters that aren't frozen
  kp=ip+lonarr(nchain) & if nsamp gt 0 then kp=reform(rpar[ib,ip,*])
  for ic=0L,nchain-1L do begin		;{chain
    if sigpar[kp[ic]] gt 0 then rr=randomn(seed,nchain) else $
	rr=2*(randomu(seed,nchain)-0.5)
    if thawed[kp[ic]] eq 0 then rr=0.
    ;
    case scltyp[kp[ic]] of
      1: begin
	xx=alog(oldpar[kp[ic],ic])
	ss=sigpar[kp[ic]]/oldpar[kp[ic],ic]
      end
      -1: begin
	xx=exp(oldpar[kp[ic],ic])
	ss=oldpar[kp[ic],ic]*sigpar[kp[ic]]
      end
      10: begin
	xx=alog10(oldpar[kp[ic],ic])
	ss=sigpar[kp[ic]]/oldpar[kp[ic],ic]/alog(10.)
      end
      -10: begin
	xx=10.^(oldpar[kp[ic],ic])
	ss=alog(10.)*oldpar[kp[ic],ic]*sigpar[kp[ic]]
      end
      100: xx=logit(oldpar[kp[ic],ic],sigpar[kp[ic]],sigl=ss)
      -100: xx=unlogit(oldpar[kp[ic],ic],sigpar[kp[ic]],sigx=ss)
      2: begin
	xx=sqrt(oldpar[kp[ic],ic])
	ss=0.5*sigpar[kp[ic]]/oldpar[kp[ic],ic]
      end
      -2: begin
	xx=oldpar[kp[ic],ic]^2
	ss=2.*oldpar[kp[ic],ic]*sigpar[kp[ic]]
      end
      else: begin
	xx=oldpar[kp[ic],ic]
	ss=sigpar[kp[ic]]
      end
    endcase
    ;
    yy=xx+rr*ss
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	;Note: the new parameter is chosen as a Gaussian deviate
	;from the old one.  This is not necessarily the best choice,
	;especially for parameters that may be strongly bounded, or
	;which have highly skewed distributions.  It is recommended
	;that an appropriate transformation on the parameter (e.g.,
	;LOGIT()) be first performed, and that this be taken into
	;account within FUNCS when the model values are calculated.
	;A rudimentary scaling correction can be carried out here
	;on the fly using the keyword SCLPAR, and if absolutely
	;unavoidable, this program can be modified to include other
	;deviates; ask the author.
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    ;
    case scltyp[kp[ic]] of
      1: testpar[kp[ic],ic]=exp(yy)
      -1: testpar[kp[ic],ic]=alog(yy)
      10: testpar[kp[ic],ic]=10.^(yy)
      -10: testpar[kp[ic],ic]=alog10(yy)
      100: testpar[kp[ic],ic]=unlogit(yy)
      -100: testpar[kp[ic],ic]=logit(yy)
      2: testpar[kp[ic],ic]=yy^2
      -2: testpar[kp[ic],ic]=sqrt(yy)
      else: testpar[kp[ic],ic]=yy
    endcase
    tmp=testpar[*,ic] & adjustie,tmp,_extra=e & testpar[*,ic]=tmp
    if keyword_set(userng) then testpar[kp[ic],ic]=$
    	(testpar[kp[ic],ic] > parrng[kp[ic],0]) < parrng[kp[ic],1]
  endfor				;IC=0,NCHAIN-1}
endfor					;II=0,NNPAR-1}

return,testpar
end
