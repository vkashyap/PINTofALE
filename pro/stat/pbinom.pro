function pbinom,q,size,prob,lotail=lotail,logp=logp,verbose=verbose, _extra=e
;+
;function	pbinom
;	for a given draw Q, and a known total sample size SIZE, computes the binomial probability of obtaining Q out of SIZE for specified probability PROB
;
;	duplicates the functionality of R's pbinom()
;
;syntax
;	Pr=pbinom(q,size,prob,/lotail,/logp,verbose=verbose)
;
;parameters
;	q	[INPUT; required] quantile values at which to compute the probability
;	size	[INPUT; required] sample size
;	prob	[INPUT; required] probability
;
;keywords
;	lotail	[INPUT; default=0] if set, computes Pr(X>q), otherwise Pr(X \le q)
;	logp	[INPUT; default=0] if set, returns ln(Pr)
;
;history
;	Vinay Kashyap (2026-04-01)
;-

ok='ok' & np=n_params() & nq=n_elements(q) & ns=n_elements(size) & nr=n_elements(prob)
if np lt 3 then ok='Insufficient parameters' else $
 if nq eq 0 then ok='Q is not defined' else $
  if ns eq 0 then ok='SIZE is not defined' else $
   if nr eq 0 then ok='PROB is not defined' else $
    ;if nq gt 1 then ok='Q cannot be an array' else $
     if ns gt 1 then ok='SIZE cannot be an array' else $
      if nr gt 1 then ok='PROB cannot be an array'
if ok ne 'ok' then begin
  print,'Usage: Pr=pbinom(q,size,prob,/lotail,/logp,verbose=verbose)'
  print,'  compute Binomial probability of obtaining Q out of sample SIZE at probability PROB'
  if np gt 0 then message,ok,/informational
  return,-1L
endif

Ntot=size[0]
Ndraw=q
pp=prob[0]

Pr=dblarr(nq)

ii=lindgen(Ntot+1L)
for i=0L,nq-1L do begin
  oi=lindgen(Ndraw[i]+1L) & if keyword_set(lotail) then oi=lindgen(Ntot-Ndraw[i]+1L)+Ndraw[i]
  wts=lngamma(Ntot+1L)-lngamma(ii[oi]+1L)-lngamma(Ntot-ii[oi]+1L)
  arr = wts + ii[oi]*alog(pp[0]) + (Ntot-ii[oi])*alog(pp[0])
  mxarr=max(arr,/nan) & arr=arr-mxarr
  Pr[i]=total(exp(arr))*exp(mxarr)
  if keyword_set(logp) then Pr[i]=alog(Pr[i])
endfor


return,Pr
end

;+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
;example call
;-------------------------------------------------------------------------------

print,pbinom()

Pr=pbinom([9,10,11],21,0.5)
print,'pbinom([9,10,11],21,0.5) = ',Pr
print,'Note: in R, pbinom(c(9,10,11),21,0.5) returns 0.3318119 0.5000000 0.6681881'

Pr=pbinom(1000,2000,0.5)
print,'pbinom(1000,2000,0.5) = ',Pr
print,'Note: in R, pbinom(1000,2000,0.5) returns 0.5089195'

end
