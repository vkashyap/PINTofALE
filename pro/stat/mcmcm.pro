function mcmcm,prb,testprb,testtyp=testtyp,seed=seed, _extra=e
;+
;function	mcmcm
;	carries out a Metropolis check on the test params
;	and either accepts (1) or rejects (0) the new set
;	based on the associated probabilities
;
;	the basic Metropolis criterion can be summarized thus
;	(see e.g., Kashyap & Drake 1998, ApJ, 503, 450):
;	if the new set of parameters are more likely, then
;	always accept them.  otherwise, accept it with a probability
;	p(TESTPRB,PRB), which ensures that the solution does not
;	get stuck in local minima.
;
;syntax
;	test=mcmcm(prb,testprb,testtyp=testtyp,seed=seed)
;
;parameters
;	prb	[INPUT; required] probability for current parameter set
;	testprb	[INPUT; required] probability for test parameter set
;		* NOTE: PRB and TESTPRB may be arrays (whose sizes
;		  _must_ match), in which case the output will be an
;		  integer array of the same size containing 0's and 1's
;		  as appropriate for each PRB,TESTPRB pair.
;
;keywords
;	testtyp	[I/O] by default, the ratio min(1,TESTPRB/PRB) is used
;		as the acceptance function.   set this keyword to
;		-- 'HELP' to print out available options
;		-- 'CHILRT' to use min(1,exp(-TESTPRB)/exp(-PRB))
;		* will be reset to 'DEFAULT' if 'HELP' is input
;	seed	[I/O] seed for random number generator RANDOMU()
;
;	_extra	[JUNK] here only to avoid crashing the program
;
;history
;	vinay kashyap (Jun2005)
;-

;	usage
ok='ok' & np=n_params() & n1=n_elements(prb) & n2=n_elements(testprb)
if np lt 2 then ok='Insufficient parameters' else $
 if n1 eq 0 then ok='PRB is undefined' else $
  if n2 eq 0 then ok='TESTPRB is undefined' else $
   if n1 ne n2 then ok='PRB and TESTPRB are incompatible'
if ok ne 'ok' then begin
  print,'Usage: test=mcmcm(prb,testprb,testtyp=testtyp)'
  print,'  accepts (1) or rejects (0) new parameter set'
  print,'  based on their associated likelihoods'
  if np ne 0 then message,ok,/informational
  return,-1L
endif

;	check keywords
typ='default'
if keyword_set(testtyp) then begin
  tt=strlowcase(testtyp[0])
  if strpos(tt,'h') eq 0 then typ='help' else $
   if strpos(tt,'chi') eq 0 then typ='chilrt'
endif
if typ eq 'help' then begin
  message,'available Acceptance functions are:',/informational
  print,'  DEFAULT:	min(1,TESTPRB,PRB)'
  print,'  CHILRT:	min(1,exp(-TESTPRB)/exp(-PRB))'
  message,'setting type to: DEFAULT',/informational
  typ='default' & testtyp=strupcase(typ)
endif

;	get random numbers
ru=randomu(seed,n1)

;	the acceptance function
case typ of
  'chilrt': begin
    accept_fn = prb-testprb < 0
    accept_fn = exp(accept_fn)
  end
  else: begin
    accept_fn = testprb/prb < 1
  end
endcase

;	the Metropolis check
if n1 eq 1 then begin
  test=0
  if finite(accept_fn) ne 0 then if ru[0] le accept_fn[0] then test=1
endif else begin
  test=intarr(n1)
  ok=where(finite(accept_fn) ne 0,mok)
  if mok gt 0 then begin
    ok2=where(ru[ok] le accept_fn[ok],mok2)
    if mok2 gt 0 then test[ok[ok2]]=1
  endif
endelse

return,test
end
