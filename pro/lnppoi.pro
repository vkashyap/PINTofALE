function lnppoi,lam,ct,cumul=cumul,verbose=verbose, _extra=e
;+
;function	lnppoi
;	computes and returns the ln(probability) of obtaining CT (or more) counts when the Poisson intensity if LAM
;
;syntax
;	lnprob=lnppoi(lambda,counts,/cumul,verbose=verbose)
;
;parameters
;	lam	[INPUT; required] Poisson intensity (can be an array)
;	ct	[INPUT] count for which to compute probability
;		* if not set, assumed to be 1
;		* if array, size must match LAM, or else only first element is used
;
;keywords
;	cumul	[INPUT] if set, computes the probability of observing CT or more counts,
;		otherwise just the probability of observing CT counts
;	verbose	[INPUT] controls chatter
;	_extra	[JUNK] here only to prevent crashing the program
;
;example
;	.run lnppoi
;	lam=      0.050000001 0.075000003 0.10000000 0.20000000 0.40000001 0.80000001 1.0000000
;	p(k.GE.1)=0.048770547 0.072256505 0.095162570 0.18126923 0.32967997 0.55067104 0.63212055
;	p(k.GE.3)=2.0027161e-05 6.6399574e-05 0.00015461445 0.0011484623 0.0079263449 0.047422647 0.080301404
;	p(k.EQ.3)=1.9817273e-05 6.5231972e-05 0.00015080623 0.0010916411 0.0071500790 0.038342736 0.061313239
;
;history
;	Vinay Kashyap (2017jul)
;-

;	usage
ok='ok' & np=n_params() & nlam=n_elements(lam) & nct=n_elements(ct)
vv=0L & if keyword_set(verbose) then vv=long(verbose[0])>1L
if np eq 0 then ok='Insufficient parameters' else $
 if nlam eq 0 then ok='LAM is not defined'
if ok ne 'ok' then begin
  print,'Usage: lnprob=lnppoi(lambda,counts,/cumul,verbose=verbose)'
  print,'  compute and return ln(Poisson probability) of obtaining CT (or more) counts for given intensity'
  if np ne 0 then message,ok,/informational
  return,-1L
endif
sct=lonarr(nlam)+1 & if nct ne 0 then begin
  if nct eq nlam then sct=ct[*] else begin
    if nct gt 1 and vv gt 0 then message,'CT and LAM sizes are not compatible; using only the first element of CT',/informational
    sct=lonarr(nlam)+ct[0]
  endelse
endif

;	output
lnprb=dblarr(nlam) & if nlam eq 1 then lnprb=0.D

;	compute the probabilities for each LAM
for i=0L,nlam-1L do begin
  zlam=lam[i] & zct=sct[i]
  if zlam eq 0 or zlam gt 50 then begin	;(in Gaussian regime
    message,'not implemented yet',/informational
  endif else begin		;ZLAM.GT.50)(Poisson regime
    if not keyword_set(cumul) then begin	;(ln(prob) for one value
      lpp=zct*alog(zlam)-zlam-lngamma(zct+1L)
      lnprb[i]=lpp
    endif else begin				;not CUMUL)(ln(prob) for all values .GE.value
      ;kmax=ceil(zlam+10*sqrt(zlam+1)) > zct
      kmax=zct
      kk=lindgen(kmax)
      lpp=kk*alog(zlam)-zlam-lngamma(kk+1L)
      lppmax=max(lpp) & lpp=lpp-lppmax & lnprb[i]=alog(1.D - total(exp(lpp))*exp(lppmax))
      if vv gt 5 then print,i,zlam,zct,lnprb[i]
    endelse					;CUMUL)
  endelse			;ZLAM.LE.50)
endfor

return,lnprb
end

;	example calling sequence

jnk=lnppoi()
lam=[0.05,0.075,0.1,0.2,0.4,0.8,1.]
print,'lam=      ',strjoin(strtrim(double(lam),2),' ')
print,'p(k.GE.1)=',strjoin(strtrim(exp(lnppoi(lam,1,/cumul)),2),' ')	;probability of seeing 1 or more counts
print,'p(k.GE.2)=',strjoin(strtrim(exp(lnppoi(lam,2,/cumul)),2),' ')	;probability of seeing 2 or more counts
print,'p(k.EQ.3)=',strjoin(strtrim(exp(lnppoi(lam,3)),2),' ')	;probability of seeing 3 counts

end

