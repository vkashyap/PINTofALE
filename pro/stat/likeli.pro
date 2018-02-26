function likeli,data,model,derror,dsigma=dsigma,ulim=ulim,$
	binom=binom,chi2=chi2,rchi=rchi,cash=cash,castor=castor,$
	softlim=softlim, _extra=e
;+
;function	likeli
;	returns p(D|M), the likelihood of observing the data for the
;	given model.
;
;	the default action is to return exp(-chi^2/2)
;
;syntax
;	prob=likeli(data,model,derror,dsigma=dsigma,ulim=ulim,$
;	/binom,/chi2,/rchi,/cash,/castor,softlim=softlim)
;
;parameters
;	data	[INPUT; required] observed data
;	model	[INPUT; required] predicted model or model parameters
;	derror	[INPUT] errors on DATA
;		* if specified, overrides keyword DSIGMA
;
;keywords
;	dsigma	[INPUT] errors on DATA
;		* if not specified, taken to be 1+sqrt(abs(DATA)+0.75)
;		* used only if DERROR is not present
;	ulim	[INPUT] long-integer array specifying which elements of
;		DATA are upper limits (1: UL, 0: not)
;		* applies only to CHI2 and RCHI
;	cash	[INPUT] if set, returns the Cash (19??, ApJ 228, 939) statistic
;	castor	[INPUT] if set, returns the Castor correction to Cash
;	binom	[INPUT] if set, likelihood is computed as
;		    MODEL(0)^(d1)*(1.-MODEL(0))^(d0)
;		where
;		    d1=total(DATA(where(DATA ge BINOM))) and
;		    d0=total(DATA(where(DATA lt BINOM)))
;	chi2	[INPUT; default] if set, returns -alog(p(D|M))
;	rchi	[INPUT] if set, returns alog(p(D|M))/n_elements(D)
;		* priority: CASH supercedes CASTOR supercedes
;		  BINOM supercedes CHI2 supercedes RCHI
;	softlim	[INPUT] if set, then likelihood is unaffected when
;		model values are below data values, and drops as a
;		Gaussian centered on SOFTLIM*DSIGMA with error DSIGMA
;		otherwise
;		* if not set, upper limits are taken to be hard
;		  limits, i.e., likelihood drops to 0 if model
;		  values exceed data values
;	_extra	[INPUT] junk -- here only to avoid crashing the program
;
;history
;	vinay kashyap (Mar 97)
;	added keywords CASH and CASTOR (VK; Jul98)
;	changed keyword name SIGMA to DSIGMA (VK; MMaug)
;	had forgotten to actually include CASH and CASTOR in call (VK; May02)
;	added parameter DERROR as a dupe for DSIGMA (VK; Jun05)
;	BUGFIX to not crash if all DATA are upper limits (VK; Sep05)
;	added keyword SOFTLIM (VK; May06)
;	bug correction re SOFTLIM model estimate (VK; Jul07)
;	corrected output to be same as all the others when /BINOM,/CHI2
;	  is set or when /RCHI is set but not /CHI2 (VK; Apr08)
;	fixed a lower bound of 1e-30 for model intensities for CASH per
;	  discussion with Jan-Uwe (VK; Aug08)
;-

;	usage
nd=n_elements(data) & nm=n_elements(model)
if nd eq 0 or nm eq 0 then begin
  print,'Usage: prob=likeli(data,model,derror,dsigma=dsigma,ulim=ulim,$'
  print,'    /cash,/castor,/binom,/chi2,/rchi,softlim=softlim)'
  print,'  returns p(D|M) or appropriate statistic'
  return,0.
endif

;	cast data
x=[data] & m=[model] & s=1.+sqrt(abs(x)+0.75) & uu=lonarr(nd)
nerr=n_elements(derror) & ns=n_elements(dsigma) & nu=n_elements(ulim)
if ns gt 0 then begin
  if ns le nd then s(0)=dsigma else s(*)=dsigma(0:nd-1)
endif
if nerr gt 0 then begin
  if nerr le nd then s(0)=derror else s(*)=derror(0:nd-1)
endif
if nu gt 0 then begin
  if nu le nd then uu(0)=ulim else uu(*)=ulim(0:nd-1)
endif
u0=1. & if keyword_set(softlim) then u0=float(softlim(0))>1
if nd ne nm then begin
  message,'DATA and MODEL do not match',/info & return,0.
endif

;	if Cash...
if n_elements(cash) gt 0 then begin
  prob=0.
  o0=where(x lt 0 or m lt 0,mo0)
  oo=where(x ge 0,moo)
  if mo0 gt 0 then begin
    message,'warning.. some bins appear to have -ve D or M',/info
  endif
  if moo eq 0 then begin
    message,'no valid data points for Cash statistic -- returning 0',/info
  endif else prob=2*total(m(oo)-x(oo)*alog((m(oo))>(1d-30))+lngamma(x(oo)+1.))
  if finite(prob) eq 0 then stop
  return, prob
endif

;	if Castor...
if n_elements(castor) gt 0 then begin
  prob=0.
  oo=where(x gt 0 and m gt 0,moo)
  if moo lt nd then begin
    message,'warning.. some bins appear to have -ve D or M',/info
  endif
  if moo eq 0 then begin
    message,'no valid data points for Cash-Castor -- returning 0',/info
  endif else prob=2*total(m(oo)-x(oo)+x(oo)*(alog(x(oo))-alog(m(oo))))
  return, prob
endif

;	if binomial...
if n_elements(binom) gt 0 then begin
  oo=where(x lt binom(0),moo)
  if moo gt 0 then d0=total(x(oo)) else d0=0.
  oo=where(x ge binom(0),moo)
  if moo gt 0 then d1=total(x(oo)) else d1=0.
  theta=abs(model(0))<1
  prob=((theta^(d1) * (1.-theta)^(d0)) < 1)>1e-30
  if keyword_set(chi2) then prob=-alog(prob)
  return,prob
endif

;	default ... exp(-chi^2/2)
oo=where(s le 0,moo)
if moo gt 0 then begin
  message,'replacing 0 and -ve errors by 1+sqrt(abs(DATA)+0.75)',$
	/informational
  s(oo)=1.+sqrt(abs(x(oo))+0.75)
endif
;	note: where uu.NE.0.AND.M.lE.X, p(D|M)=1
ou1=where(uu ne 0 and m le x,mou1)
ou2=where(uu ne 0 and m gt x,mou2)
oo=where(uu eq 0 and s gt 0,moo)
if moo gt 0 then chi22=total(((x(oo)-m(oo))/s(oo))^2)/2. else chi22=0.
if keyword_set(softlim) then begin
  ;if mou2 gt 0 then m(ou2)=u0*s(ou2)	;THIS IS WRONG!  model is getting reset!
  if mou2 gt 0 then chi2u=total(((x(ou2)-m(ou2))/u0/s(ou2))^2)/2. else chi2u=0.
  chi22=chi22+chi2u
endif else begin
  if mou2 gt 0 then begin
    if keyword_set(chi2) then return,2d30 else return,1d-30
  endif
endelse
if keyword_set(rchi) then chi22=chi22/(moo>1)	;reduced chi^2
if not keyword_set(chi2) and not keyword_set(rchi) then begin
  chi22=chi22 < 69
  prob=exp(-chi22) < 1
endif else prob=chi22

;oo=where(uu ne 0 and m gt x,moo)
;if moo gt 0 then begin			;p(D|M)=0: model violates upper limit
;  if keyword_set(softlim) then begin
;    m(oo)=u0*s(oo) & uu(oo)=0
;  endif else begin
;    if keyword_set(chi2) then return,2d30 else return,1d-30
;  endelse
;endif
;oo=where(uu eq 0 and s gt 0,moo)
;if moo gt 0 then begin
;  chi22=total(((x(oo)-m(oo))/s(oo))^2)/2.
;  if keyword_set(rchi) then chi22=chi22/moo	;reduced chi^2
;  if not keyword_set(chi2) then begin
;    chi22=chi22 < 69
;    prob=exp(-chi22) < 1
;  endif else prob=chi22
;endif else begin
;  if keyword_set(chi2) then return,0. else return,1.
;endelse

return,prob
end
