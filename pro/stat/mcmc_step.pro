function mcmc_step,x,y,par,sigpar,sigy=sigy,funcs=funcs,fnprob=fnprob,$
	nbatch=nbatch,rngpar=rngpar,thaw=thaw,singly=singly,sampar=sampar,$
	sclpar=sclpar,seed=seed, _extra=e
;+
;function	mcmc_step
;	returns an updated set of parameters in Markov-Chain Monte Carlo step
;
;syntax
;	newpar=mcmc_step(x,y,par,sigpar,sigy=sigy,funcs=funcs,fnprob=fnprob,$
;	nbatch=nbatch,rngpar=rngpar,thaw=thaw,/singly,sampar=sampar,$
;	sclpar=sclpar,seed=seed,testtyp=testtyp,ties=ties,$
;	FUNCS;	type=type,/fwhm,/norm,betap=betap,vrot=vrot,angle=angle,$
;		phase=phase,group=group,delp=delp,missing=missing,$
;	FNPROB;	ulim=ulim,/chi2,/binom,/cash,/castor)
;
;parameters
;	x	[INPUT; required] data points
;	y	[INPUT; required] Y(X)
;		* sizes of X and Y _must_ match
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
;	sigy	[INPUT] error on Y
;		* by default, assumed to be Gehrel's approximation,
;		  1+sqrt(abs(Y)+0.75)
;		* if single element, assumed to be
;		  -- fractional error on Y if +ve
;		  -- abs value is absolute error if -ve
;		* if size does not match Y, default is used
;	funcs	[INPUT] name of user defined function (actually a procedure,
;		but this name is used for compatibility with IDL's curvefit)
;		that takes as input X and A (the model parameters), and
;		returns Y(X;A) (the model)
;		* should be callable from the command line independently as
;		  FUNCS, X, PAR, YMODEL
;		* default is set to X3MODEL, which accepts keywords
;		  TYPE=TYPE,/FWHM,/NORM,BETAP=BETAP,VROT=VROT,ANGLE=ANGLE,$
;		  PHASE=PHASE,GROUP=GROUP,DELP=DELP,MISSING=MISSING
;	fnprob	[INPUT] name of user defined function that computes the
;		probability of model parameters
;		* must be a function that takes as input parameters
;		  data,model,errors and return probability, as follows:
;		  PROB = FNPROB( Y, MODEL, SIGY )
;		* default is LIKELI(), which accepts keywords
;		  ULIM=ULIM,/CHI2,/BINOM,/CASH,/CASTOR
;		* note: if using /CHI2, remember to also set TESTTYP='CHILRT'
;	nbatch	[INPUT; default=10] the number of new parameter sets to try 
;		before exiting the program.  only the last set is returned
;		as output and the rest are discarded.
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
;	singly	[INPUT] if set, does the Metropolis-Hastings check
;		for each parameter individually, before going on to
;		the next parameter
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
;	_extra	[INPUT ONLY] pass defined keywords to subroutines:
;		testtyp	[MCMCM]
;		ties	[ADJUSTIE]
;		*	[FUNCS]
;		*	[FNPROB]
;
;restrictions
;	requires user-defined procedure FUNCS and function FNPROB()
;	(see keyword descriptions above)
;
;subroutines
;	ADJUSTIE
;	LOGIT()
;	MCMCM()
;	UNLOGIT()
;	-FUNCS-
;	-FNPROB()-
;
;history
;	vinay kashyap (Jun2005)
;-

message,'WARNING: STILL IN ALPHA!',/informational

;	usage
ok='ok' & np=n_params() & nx=n_elements(x) & ny=n_elements(y)
szp=size(par) & nszp=n_elements(szp) & npar=szp[1] & npe=n_elements(sigpar)
if np lt 4 then ok='Insufficient parameters' else $
 if nx eq 0 then ok='X is undefined' else $
  if ny eq 0 then ok='Y(X) is undefined' else $
   if nx ne ny and nx+1L ne ny then ok='X and Y(X) are incompatible' else $
    if szp[0] eq 0 then ok='PAR are undefined' else $
     if szp[0] gt 2 then ok='cannot understand PAR > 2-D' else $
      if npe eq 0 then ok='SIGPAR are undefined' else $
       if npe ne npar then ok='PAR and SIGPAR are incompatible'
if ok ne 'ok' then begin
  print,'Usage: newpar=mcmc_step(x,y,par,sigpar,sigy=sigy,funcs=funcs,fnprob=fnprob,$'
  print,'       nbatch=nbatch,rngpar=rngpar,thaw=thaw,/singly,sampar=sampar,$'
  print,'       sclpar=sclpar,seed=seed,testtyp=testtyp,ties=ties)'
  print,'  returns an updated set of parameters in Markov-Chain Monte Carlo step'
  if np ne 0 then message,ok,/informational
  return,-1L
endif

;	check input
  ;
;  user-defined model generating procedure
if not keyword_set(funcs) then funcs='x3model'
  ;
;  user-defined probability calculating function
if not keyword_set(fnprob) then fnprob='likeli'
  ;
;  number of chains
if szp[0] eq 2 then nchain=szp[2] else nchain=1
  ;
;  number of batches
mbatch=10L & if keyword_set(nbatch) then mbatch=long(nbatch[0])>1
  ;
;  errors on Y
ysig=1.+sqrt(abs(y)+0.75) & nye=n_elements(ysig)
if nye eq ny then ysig=abs(sigy)
if nye ne 0 and nye ne ny then begin
  if nye eq 1 then begin
    if sigy[0] gt 0 then ysig=sigy[0]*ysig
    if sigy[0] lt 0 then ysig[*]=abs(sigy[0])
  endif
  if nye gt 1 then message,$
	"SIGY is incompatible with Y; using Gehrel's approximation",$
	/informational
endif
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
dsamp=fltarr(npar)+1./npar & nsamp=n_elements(sampar)
if nsamp gt 1 then begin
  ksamp = nsamp < npar
  dsamp[0L:ksamp-1L]=sampar[0L:ksamp-1L]/total(sampar[0L:ksamp-1L])
endif
sampfn=fltarr(npar+1L)
for i=1L,npar do sampfn[i]=sampfn[i-1L]+dsamp[i-1L]
sampfn=sampfn/max(sampfn)
rsamp=reform(randomu(seed,nbatch*npar*nchain),nbatch,npar,nchain)
rpar=long(interpol(lindgen(npar+1L),sampfn,rsamp))
  ;
;  parameter scaling transformations
scltyp=intarr(npar) & nscl=n_elements(sclpar)
if nscl eq 1 then scltyp[*]=sclpar[0]
if nscl eq npar then scltyp=sclpar

;	get new set of parameters
oldpar=par & newpar=par
for ib=0L,mbatch-1L do begin		;{batch
  testpar=oldpar & prb=fltarr(nchain)
  for ic=0L,nchain-1L do begin
    call_procedure,funcs,x,oldpar[*,ic],ymod, _extra=e
    prb[ic]=call_function(fnprob,y,ymod,ysig, _extra=e)
  endfor
  for ii=0L,nnpar-1L do begin		;{thawed parameters

    ip=opar[ii]		;look only at parameters that aren't frozen

    ;	---- if you must call MCMC_DEVIATE(), then do so here ----

    kp=ip+lonarr(nchain) & if nsamp gt 0 then kp=reform(rpar[ib,ip,*])
    for ic=0L,nchain-1L do begin	;{in each chain
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
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      ;  Note: the new parameter is chosen as a Gaussian deviate     ;
      ;  from the old one.  This is not necessarily the best choice, ;
      ;  especially for parameters that may be strongly bounded, or  ;
      ;  which have highly skewed or correlated distributions.  It   ;
      ;  is recommended that an appropriate transformation on the    ;
      ;  parameter (e.g., LOGIT()) be first performed, and that this ;
      ;  be taken into account within FUNCS when the model values    ;
      ;  are calculated.  A rudimentary scaling correction can be    ;
      ;  carried out here on the fly using the keyword SCLPAR.  If   ;
      ;  completely and absolutely unavoidable, this program can be  ;
      ;  modified to include other deviates; ask the author.         ;
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
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
    if keyword_set(singly) then begin
      testprb=prb
      for ic=0L,nchain-1L do begin
        call_procedure,funcs,x,testpar[*,ic],ymod, _extra=e
	testprb[ic]=call_function(fnprob,y,ymod,ysig, _extra=e)
      endfor
      test=mcmcm(prb,testprb,seed=seed,_extra=e)
      for ic=0L,nchain-1L do if test[ic] ne 0 then $
	oldpar[kp[ic],ic]=testpar[kp[ic],ic]
    endif
  endfor				;II=0,NNPAR-1}
  if not keyword_set(singly) then begin
    testprb=prb
    for ic=0L,nchain-1L do begin
      call_procedure,funcs,x,testpar[*,ic],ymod, _extra=e
      testprb[ic]=call_function(fnprob,y,ymod,ysig, _extra=e)
    endfor
    test=mcmcm(prb,testprb,seed=seed,_extra=e)
    for ic=0L,nchain-1L do if test[ic] ne 0 then oldpar[*,ic]=testpar[*,ic]
  endif
endfor					;IB=0,MBATCH-1}

return,newpar
end
