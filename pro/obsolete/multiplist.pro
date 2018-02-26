function multiplist,x1,x2,xout,nbin=nbin,xmin=xmin,xmax=xmax,$
	nonorm=nonorm, _extra=e
;+
;function	multiplist
;	return the product of the frequency distributions of two lists
;
;syntax
;	m12=multiplist(x1,x2,xout,nbin=nbin,xmin=xmin,xmax=xmax,/nonorm)
;
;parameters
;	x1	[INPUT; required] list of numbers whose frequency
;		distribution must be multiplied with that of X2
;	x2	[INPUT; required] list of numbers whose frequency
;		distribution must be multiplied with that of X1
;		* it is assumed that X1 and X2 both have the same units
;		  because otherwise this kind of multiplication makes
;		  no sense
;	xout	[OUTPUT; required] the list of numbers at which the
;		output is tabulated
;		* by default, this is the sorted unique set of [X1,X2]
;		* if NBIN is set, generates a logarithmic or linear
;		  grid of that size
;
;keywords
;	nbin	[INPUT] number of bins in the output
;		* if not set, will be the sorted unique set of [X1,X2]
;		* if -ve, will produce logarithmic gridding
;	xmin	[INPUT] minimum value in the output grid
;		* by default, uses min(X1)>min(X2)
;	xmax	[INPUT] maximum value in the output grid
;		* by default, uses max(X1)<max(X2)
;	nonorm	[INPUT] if set, does not make a correction for the
;		width of the bins in XOUT.
;	_extra	[JUNK] here only to prevent crashing the program
;
;description
;	this program does not create histograms as intermediate products
;	but rather produces a multiplicaton at the best available data
;	resolution.  this is of great use in Monte Carlo calculations
;	or in creating probability distributions in MCMC.  the algorithm
;	is straightforward: build up a cdf for each list, and interpolate
;	onto a common grid, and multiply the d(cdf)'s, suitably normalized
;	for the bin widths and number of elements.
;
;history
;	vinay kashyap (Jun2005)
;-

message,'OBSOLTE; use ARITHTOGRAM() instead.',/informational

;	usage
ok='ok' & np=n_params() & n1=n_elements(x1) & n2=n_elements(x2)
if np lt 3 then ok='Insufficient parameters' else $
 if n1 eq 0 then ok='X1 is undefined' else $
  if n2 eq 0 then ok='X2 is undefined' else $
   if n1 eq 1 then ok='X1 must have at least 2 elements' else $
    if n2 eq 1 then ok='X2 must have at least 2 elements'
if ok ne 'ok' then begin
  print,'Usage: m12=multiplist(x1,x2,xout,nbin=nbin,xmin=xmin,xmax=xmax,/nonorm)'
  print,'  return the product of the frequency distributions of two lists'
  if np ne 0 then message,ok,/informational
  return,-1L
endif

;	keywords
xmin=min(X1) > min(X2)
xmax=max(X1) < max(X2)
mbin=0L & if keyword_set(nbin) then mbin=long(nbin[0])

;	output grid
if keyword_set(mbin) then begin
  if mbin lt 0 then begin
    lxmin=alog10(xmin) & lxmax=alog10(xmax)
    dx=(lxmax-lxmin)/(abs(mbin)-1L)
    xout=10.^(findgen(abs(mbin)+1L)*dx+lxmin)
  endif else begin
    dx=(xmax-xmin)/(mbin-1L)
    xout=findgen(mbin+1L)*dx+xmin
  endelse
endif else begin
  xout=[X1,X2] & xout=xout[uniq(xout,sort(xout))]
  ok=where(xout ge xmin and xout le xmax,mok)
  if mok eq 0 then message,'BUG!'
  xout=xout[ok]
endelse

;	make cdfs
o1=sort(X1) & o2=sort(X2)
c1=findgen(n1)/(n1-1L) & c2=findgen(n2)/(n2-1L)
cc1=interpol(c1,x1[o1],xout) & cc2=interpol(c2,x2[o2],xout)
dc1=n1*(cc1[1:*]-cc1) & dc2=n2*(cc2[1:*]-cc2)
dxout=xout[1:*]-xout & xout=0.5*(xout[1:*]+xout)
m12=dc1*dc2
if not keyword_set(nonorm) then m12=m12/dxout^2

return,m12
end
