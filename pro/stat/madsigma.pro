function madsigma,x,medianx=medianx,madx=madx
;+
;function	madsigma
;	compute the Median Absolute Deviation width of a distribution from a sample
;
;	See Xu et al. 2021, AJ 161, 184, Sec 2.4.1 and footnote 5
;
;syntax
;	sig=madsigma(x,medianx=medianx,madx=madx)
;
;parameters
;	X	[INPUT; required] sample of values defining a distribution
;
;keywords
;	medianX	[OUTPUT] median of X
;	madX	[OUTPUT] median of the absolute deviation from medianX
;
;history
;	Vinay Kashyap (MMXXVI.II)
;-

ok='ok' & np=n_params() & nx=n_elements(X)
if np eq 0 then ok='Insufficient parameters' else $
 if nx eq 0 then ok='X is not defined' else $
  if nx eq 1 then ok='X must be an array of size>1'
if ok ne 'ok' then begin
  print,'Usage: sigma=madsigma(X,MEDIANX=MEDIANX,MADX=MADX)
  print,'  compute a robust width of a distribution from a sample drawn from it'
  if np ne 0 then message,ok,/informational
  return,-1L
endif

medianx=median(x)
madx=median(abs(x-medianx))
sig=madx/0.6745

return,sig
end
