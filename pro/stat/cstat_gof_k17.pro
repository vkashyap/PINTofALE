function cstat_gof_k17,mval,obsstat,pval,C_e=C_e,S_nu=S_nu,C_nu=C_nu,zval=zval,verbose=verbose, _extra=e
;+
;function	cstat_gof_k17
;	compute and return the expected cstat value and its variance as
;	an array [expected cstat, variance] for a given set of model intensities
;
;	an implementation of the expressions for C_e, S_\nu and C_\nu from Sec 3
;	of Kaastra 2017 (A&A 605, A51)
;
;syntax
;	cstat_calc = cstat_gof_k17(mval,obsstat,pval,C_e=C_e,S_nu=S_nu,C_nu=C_nu,zval=zval,verbose=verbose)
;
;parameters
;	mval	[INPUT; required] model intensity values in units of [ct/bin]
;	obsstat	[INPUT] the observed value of cstat
;	pval	[OUTPUT] if OBSSTAT is given, the p-value of the observed cstat
;		* computed as the integral of N(cstat_expected,cstat_variance) for
;		  cstat_expected>cstat_observed
;		* numbers close to 0 mean the observed value is way in the tail of
;		  the expected distribution
;		* numbers close of 1 mean something has gone badly wrong, that either
;		  the observed values or the model values are incorrect
;
;keywords
;	C_e	[OUTPUT] expected contribution of cstat to each bin (Eq 4 of Kaastra 2017)
;	C_nu	[OUTPUT] expected variance of cstat in each bin (Eq 6 of Kaastra 2017)
;	S_nu	[OUTPUT] intermediate array (Eq 5 of Kaastra 2017)
;	zval	[OUTPUT] (OBSSTAT-CSTAT_EXPECTED)/sqrt(CSTAT_VARIANCE)
;	verbose	[INPUT] controls chatter
;	_extra	[JUNK] here only to prevent crashing the program
;
;history
;	Vinay Kashyap (2018-jul-29)
;	bug fix (VK; 2018-aug-30)
;-

;	usage
ok='ok' & np=n_params() & nm=n_elements(mval) & nobs=n_elements(obsstat)
if np lt 1 then ok='Insufficient parameters' else $
 if nm eq 0 then ok='Model values are not defined'
if ok ne 'ok' then begin
  print,'Usage: cstat_calc = cstat_gof_k17(mval,obsstat,pval,C_e=C_e,S_nu=S_nu,C_nu=C_nu,zval=zval,verbose=verbose)
  print,'  compute and return [expected cstat, variance]'
  print,'  and optionally a p-value for how good the observed cstat is'
  if np ne 0 then message,ok,/informational
  return,-1L
endif

vv=0L & if keyword_set(verbose) then vv=long(verbose[0])>1

;	a straightforward calculation
o00=where(mval le 0,mo00)
o08=where(mval gt 0 and mval le 0.5,mo08)
o09=where(mval gt 0.5 and mval le 2,mo09)
o10=where(mval gt 2 and mval le 5,mo10)
o11=where(mval gt 5 and mval le 10,mo11)
o12=where(mval gt 10,mo12)
;o13=where(mval gt 0 and mval le 0.1,mo13)
o13=where(mval gt 0,mo13)	;because this is only used to compute S_nu using the standard formula, so why not the whole thing?  We're gonna overwrite the c_nu for all mval>0.1 anyway
o14=where(mval gt 0.1 and mval le 0.2,mo14)
o15=where(mval gt 0.2 and mval le 0.3,mo15)
o16=where(mval gt 0.3 and mval le 0.5,mo16)
o17=where(mval gt 0.5 and mval le 1,mo17)
o18=where(mval gt 1 and mval le 2,mo18)
o19=where(mval gt 2 and mval le 3,mo19)
o20=where(mval gt 3 and mval le 5,mo20)
o21=where(mval gt 5 and mval le 10,mo21)
o22=where(mval gt 10,mo22)


c_e = 0.*mval & s_nu = c_e & c_nu = s_nu

if mo00 gt 0 then begin & c_e[o00] = 0. & s_nu[o00]=0. & c_nu[o00]=0. & endif	;this catches incidences of alog(0)'s
if mo08 gt 0 then c_e[o08] = -0.25*mval[o08]^3 + 1.38*mval[o08]^2 - 2.*mval[o08]*alog(mval[o08])
if mo09 gt 0 then c_e[o09] = -0.00335*mval[o09]^5 + 0.04259*mval[o09]^4 - 0.27331*mval[o09]^3 + 1.381*mval[o09]^2 - 2.*mval[o09]*alog(mval[o09])
if mo10 gt 0 then c_e[o10] = 1.019275 + 0.13145*exp((0.461-0.9*alog(mval[o10]))*alog(mval[o10]))
if mo11 gt 0 then c_e[o11] = 1.00624 + 0.604*exp(-1.68*alog(mval[o11]))
if mo12 gt 0 then c_e[o12] = 1. + 0.1649/mval[o12] + 0.226/mval[o12]^2

if mo13 gt 0 then begin
  for kk=0,4 do begin
    lnpoi = -mval[o13]+kk*alog(mval[o13])-lngamma(kk+1.)
    tmp=mval[o13] & if kk gt 0 then tmp=mval[o13]-kk+kk*alog(kk)-kk*alog(mval[o13])
    s_nu[o13] = s_nu[o13] + exp(lnpoi)*tmp^2
  endfor
  s_nu[o13] = 4.*s_nu[o13]
  c_nu[o13] = s_nu[o13] - c_e[o13]^2
endif

if mo14 gt 0 then c_nu[o14] = -262.*mval[o14]^4 + 195.*mval[o14]^3 - 51.24*mval[o14]^2 + 4.34*mval[o14] + 0.77005
if mo15 gt 0 then c_nu[o15] = 4.23*mval[o15]^2 - 2.9254*mval[o15] + 1.12522
if mo16 gt 0 then c_nu[o16] = -3.7*mval[o16]^3 + 7.328*mval[o16]^2 - 3.6926*mval[o16] + 1.20641
if mo17 gt 0 then c_nu[o17] = 1.28*mval[o17]^4 - 5.191*mval[o17]^3 + 7.666*mval[o17]^2 - 3.5446*mval[o17] + 1.15431
if mo18 gt 0 then c_nu[o18] = 0.1125*mval[o18]^4 - 0.641^mval[o18]^3 + 0.859*mval[o18]^2 + 1.0914*mval[o18] - 0.05748
if mo19 gt 0 then c_nu[o19] = 0.089*mval[o19]^3 - 0.872*mval[o19]^2 + 2.8422*mval[o19] - 0.67539
if mo20 gt 0 then c_nu[o20] = 2.12336 + 0.012202*exp( (5.717-2.6*alog(mval[o20]))*(alog(mval[o20])) )
if mo21 gt 0 then c_nu[o21] = 2.05159 + 0.331*exp( (1.343-alog(mval[o21]))*(alog(mval[o21])) )
if mo22 gt 0 then c_nu[o22] = 12./mval[o22]^3 + 0.79/mval[o22]^2 + 0.6747/mval[o22] + 2.

cstat_expected = total(c_e,/double)
cstat_variance = total(c_nu,/double)
cstat_stddev = sqrt(cstat_variance)

cstat_calc = [cstat_expected,cstat_variance]

if nobs gt 0 then begin
  zval = double(obsstat-cstat_expected)/cstat_stddev
  op=where(zval ge 0,mop,complement=om,ncomplement=mom)
  pval = 0.D * zval
  if mop gt 0 then pval[op]=erfc(zval[op]/sqrt(2.D))
  if mom gt 0 then pval[om]=erf(abs(zval[om])/sqrt(2.D))
  if nobs eq 1 then pval=reform(pval[0])	;convert to scalar
endif

if vv gt 1000 then stop,'type .CON to continue'

return,cstat_calc
end
