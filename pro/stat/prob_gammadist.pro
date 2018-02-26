function prob_gammadist,x,agamma=agamma,bgamma=bgamma,verbose=verbose, _extra=e
;+
;function	prob_gammadist
;	compute the gamma distribution, which is a continuous distribution
;	on the real line with density function
;		f(x) = b^a x^(a-1) exp(-bx) / Gamma(a)
;	(see eqn 4 of van Dyk et al., 2001, ApJ 548, 224)
;	this distribution has mean (a/b) and variance (a/b^2) for a,b>0
;
;	This is equivalent to a Poisson distribution, and is often
;	used as a prior function, with a and b representing our prior
;	expectation that a-1 counts should be observed in an area
;	b*src_area; a non-informative prior therefore has a=1 and b=0,
;	corresponding to 0 counts in 0 area.
;
;syntax
;	g=prob_gammadist(x,agamma=agamma,bgamma=bgamma,verbose=verbose)
;
;parameters
;	x	[INPUT; required] array of abscissa values at which to
;		compute the Gamma distribution
;		* for x < 0, the function will return NaNs
;
;keywords
;	agamma	[INPUT; default=1] "alpha" parameter for the function
;		* AGAMMA must be > 0, or else the function will return Infs
;	bgamma	[INPUT; default=0] "beta" parameter for the function
;		* BGAMMA should be .GE. 0, or else who knows what will happen
;	verbose	[INPUT; default=0] controls chatter
;	_extra	[JUNK] here only to prevent crashing the program
;
;history
;	vinay kashyap (Aug01)
;-

;	usage
ok='ok' & np=n_params() & nx=n_elements(x)
if np eq 0 then ok='Insufficiant parameters' else $
 if nx eq 0 then ok='X is not defined'
if ok ne 'ok' then begin
  print,'Usage: g=prob_gammadist(x,agamma=agamma,bgamma=bgamma,verbose=verbose)'
  print,'  compute the gamma distribution,'
  print,'  g(x) = bgamma^agamma x^(agamma-1) exp(-bgamma x) / Gamma(agamma)'
  if np ne 0 then message,ok,/info
  return,-1L
endif

;	keywords
aa=1. & bb=0.
if n_elements(agamma) ne 0 then aa=agamma[0]+0. else message,'Setting A=1',/info
if n_elements(bgamma) ne 0 then bb=bgamma[0]+0. else message,'Setting B=0',/info
vv=0 & if keyword_set(verbose) then vv=long(verbose[0]) > 1

;	reformat input
xx=[x[*]]
oo=where(xx lt 0,moo) & ok=where(xx gt 0,mok) & o0=where(xx eq 0,mo0)

;	initialize output
if moo ne 0 and vv gt 0 then message,strtrim(moo,2)+$
	' points will be set to NaN',/info
g=fltarr(nx) & if moo ne 0 then g[oo]=!values.F_NAN

;	compute function
if bb ne 0 then tmp1=aa*alog(bb) else tmp1=0.
tmp2=fltarr(nx) & if mok gt 0 then tmp2[ok]=(aa-1.)*alog(xx[ok])
tmp3=fltarr(nx) & if mok gt 0 then tmp3[ok]=-bb*xx[ok]
if mok gt 0 then g[ok]=tmp1+tmp2[ok]+tmp3[ok]-lngamma(aa)
if mo0 gt 0 then g[o0]=tmp1+tmp2[o0]+tmp3[o0]-lngamma(aa)
if mok gt 0 then g[ok]=exp(g[ok])
if mo0 gt 0 then g[o0]=exp(g[o0])
if aa ne 1 and mo0 gt 0 then g[o0]=0.

if vv gt 20 and !d.NAME eq 'X' then plot,x,g,xtitle='!3X',$
	ytitle='!3G(X;!4a,b!3)',title='!4a!3='+strtrim(aa,2)+' ; !4b!3='+$
	strtrim(bb,2)
if vv gt 100 then stop,'Halting.  type .CON to continue'

return,g
end
