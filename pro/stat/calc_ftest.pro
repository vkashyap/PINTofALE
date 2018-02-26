function calc_ftest,dof1,chi1,dof2,chi2,fstat=fstat,wdchi=wdchi,$
	_extra=e
;+
;function	calc_ftest
;	calculate the significance of a new model component using the
;	F test.  Returns the p-value that the observed change in the
;	chisq values for the given degrees of freedom can be obtained
;	as a fluctuation.
;
;syntax
;	pval=calc_ftest(dof1,chi1,dof2,chi2,fstat=fstat,/wdchi, /double)
;
;parameters
;	dof1	[INPUT; required] degrees of freedom of SIMPLE model
;	chi1	[INPUT; required] min chisq of SIMPLE model
;	dof2	[INPUT; required] degrees of freedom of COMPLEX model
;	chi2	[INPUT; required] min chisq of COMPLEX model
;		* if arrays, the sizes of all of them must match.
;		  the output will then be an array of the same size.
;		* NOTE: complex model _must_ be a subset of simple model
;		  and the parameter values that define one of the models
;		  should not lie on the boundary of the other.
;
;keywords
;	fstat	[OUTPUT] the value of the F statistic
;		* (CHI1/DOF1)/(CHI2/DOF2)
;		* larger the FSTAT, the more plausible the complex model
;	wdchi	[INPUT] if set, uses a different F statistic, defined as
;		(DELTA_CHI/DELTA_DOF)/(CHI2/DOF2) =
;			((CHI1-CHI2)/(DOF2-DOF1))/(CHI2/DOF2)
;		* WARNING: if DOF1 is smaller than DOF2, that is, if
;		  it appears that the input is already in this format,
;		  additional subtractions are not done.
;
;	_extra	[INPUT ONLY] pass defined keywords to IBETA()
;
;description
;	this is how CIAO/Sherpa's calc_ftest program works.
;	I just implemented it as a quick way to get at the
;	same results in IDL.
;
;	Avoid using the F test unless you know that it is
;	applicable for the problem you are working on.  See
;	Protassov et al., 2002, ApJ, 571, 545
;	  http://adsabs.harvard.edu/abs/2002ApJ...571..545P
;
;history
;	vinay kashyap (MMXII.VII)
;-

;	usage
ok='ok' & np=n_params()
nd1=n_elements(dof1) & nc1=n_elements(chi1)
nd2=n_elements(dof2) & nc2=n_elements(chi2)
if np lt 4 then ok='Insufficient parameters' else $
 if nd1 eq 0 then ok='DOF1 is not defined' else $
  if nc1 eq 0 then ok='CHI1 is not defined' else $
   if nd2 eq 0 then ok='DOF2 is not defined' else $
    if nc2 eq 0 then ok='CHI2 is not defined' else $
     if nd1 ne nc1 then ok='DOF1 and CHI1 are not compatible' else $
      if nd1 ne nd2 then ok='DOF1 and DOF2 are not compatible' else $
       if nd1 ne nc2 then ok='DOF1 and CHI2 are not compatible'
if ok ne 'ok' then begin
  print,'Usage: pval=calc_ftest(dof1,chi1,dof2,chi2,ftest=ftest,/wdchi,$'
  print,'       /double)'
  print,'  compute significance of F statistic'
  if np ne 0 then message,ok,/informational
  return,-1L
endif

;	outputs
fstat=1e6 & pval=0.
if nd1 gt 1 then begin
  fstat=0.*dof1+1e6 & pval=0.*dof1
endif

for i=0L,nd1-1L do begin		;{for all rows of input

  ;	inputs
  d1=dof1[i] & c1=1.0*chi1[i] & d2=dof2[i] & c2=1.0*chi2[i]

  ;	this has been commented out because of possible conflict with /WDCHI
  ;if d2 gt d1 then begin
  ;  message,'the complex model cannot have more '+$
  ;	'degrees of freedom than the simpler one',/informational
  ;  message,'reversing arguments for row '+strtrim(i,2),/informational
  ;  dd=d1 & cc=c1 & d1=d2 & c1=c2 & d2=dd & c2=cc
  ;endif

  if keyword_set(wdchi) then begin
    dd=(d1-d2) & cc=(c1-c2)	;computes (dCHI/dDOF)/(CHI2/DOF2)
    if dd lt 0 then begin
      d1=dd+d2 & c1=cc+c2
    endif else begin
      d1=dd & c1=cc
    endelse
    ;print,'pval=calc_ftest('+strtrim(d1,2)+','+strtrim(c1,2)+','+strtrim(d2,2)+','+strtrim(c2,2)+')'
  endif

  ;	and.. this is it
  if d1 gt 0 and c2 gt 0 then begin
    fstat[i] = (c1*d2)/(d1*c2)	;==(chi1/dof1)/(chi2/dof2)
    x=d2/(d1*fstat[i]+d2)
    pval[i]=ibeta(d2/2.,d1/2.,x, _extra=e)
  endif

endfor					;I=0,ND1-1}

return,pval
end
