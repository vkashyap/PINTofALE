function arithtogram,x1,x2,xout,oper,plus=plus,minus=minus,divide=divide,$
	w1=w1,w2=w2,nbin=nbin,xmin=xmin,xmax=xmax,nonorm=nonorm, _extra=e
;+
;function	arithtogram
;	return the result of an arithmetical operation on the
;	frequency distributions of two lists
;
;syntax
;	hx=arithtogram(x1,x2,xout,operator,/plus,/minus,/divide,$
;	w1=w1,w2=w2,nbin=nbin,xmin=xmin,xmax=xmax,/nonorm)
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
;	oper	[INPUT] must be one of '+', '-', '/', '*':
;		the arithmetical operation to carry out on the
;		frequency distributions of X1 and X2, i.e.,
;			result = X1 OPER X2
;		* if not given, assumed to be '*', unless one of
;		  keywords PLUS, MINUS, or DIVIDE are set 
;
;keywords
;	plus	[INPUT] if set, OPER is assumed to be '+'
;	minus	[INPUT] if set, OPER is assumed to be '-'
;	divide	[INPUT] if set, OPER is assumed to be '/'
;		* DIVIDE takes precedence over MINUS takes precedence
;		  over PLUS
;	w1	[INPUT] optional weight to be applied to histogram(X1)
;	w2	[INPUT] optional weight to be applied to histogram(X2)
;		* if not given, W1 and W2 are assumed to be 1
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
;	but rather produces an operation at the best available data
;	resolution.  this is of great use in Monte Carlo calculations
;	or in creating probability distributions in MCMC.  the algorithm
;	is straightforward: build up a cdf for each list, and interpolate
;	onto a common grid, and multiply, divide, add, or subtract the
;	d(cdf)'s, suitably normalized for the bin widths and number of
;	elements.
;
;example
;	x1=randomn(seed,10000L) & x2=randomn(seed,10000L)+2.
;	xmin=-2. & xmax=3. & nbin=101L
;	hh=arithtogram(x1,x2,xout,'*',xmin=xmin,xmax=xmax,nbin=nbin,/nonorm)
;	hp=arithtogram(x1,x2,xout,'+',xmin=xmin,xmax=xmax,nbin=nbin,/nonorm)
;	hd=arithtogram(x1,x2,xout,'/',xmin=xmin,xmax=xmax,nbin=nbin,/nonorm)
;	hm=arithtogram(x1,x2,xout,'-',xmin=xmin,xmax=xmax,nbin=nbin,/nonorm)
;	h1=histogram(x1,min=xmin,max=xmax,binsize=median(xout[1:*]-xout))
;	h2=histogram(x2,min=xmin,max=xmax,binsize=median(xout[1:*]-xout))
;	plot,xout,hh,psym=10 & oplot,xout,h1,col=2 & oplot,xout,h2,col=3
;		oplot,xout,h1*h2,psym=10,col=4
;	plot,xout,hp,psym=10 & oplot,xout,h1,col=2 & oplot,xout,h2,col=3
;		oplot,xout,h1+h2,psym=10,col=4
;	plot,xout,hd,psym=10 & oplot,xout,h1,col=2 & oplot,xout,h2,col=3
;		oplot,xout,float(h1)/float(h2),psym=10,col=4
;	plot,xout,hm,psym=10 & oplot,xout,h1,col=2 & oplot,xout,h2,col=3
;		oplot,xout,h1-h2,psym=10,col=4
;
;history
;	vinay kashyap (Jun2005)
;-

;	usage
ok='ok' & np=n_params() & n1=n_elements(x1) & n2=n_elements(x2)
if np lt 3 then ok='Insufficient parameters' else $
 if n1 eq 0 then ok='X1 is undefined' else $
  if n2 eq 0 then ok='X2 is undefined' else $
   if n1 eq 1 then ok='X1 must have at least 2 elements' else $
    if n2 eq 1 then ok='X2 must have at least 2 elements'
if ok ne 'ok' then begin
  print,'Usage: hx=arithtogram(x1,x2,xout,oper,/plus,/minus,/divide,$'
  print,'       w1=w1,w2=w2,nbin=nbin,xmin=xmin,xmax=xmax,/nonorm)'
  print,'  return the result of an arithmetic operation on the'
  print,'  frequency distributions of two lists'
  if np ne 0 then message,ok,/informational
  return,-1L
endif

;	operator
op='*'
if np lt 4 then oper='*'
if strpos(oper[0],'+') ge 0 or $
	strpos(strlowcase(oper[0]),'plu') ge 0 or $
	strpos(strlowcase(oper[0]),'add') ge 0 or $
	keyword_set(plus) then op='+'
if strpos(oper[0],'-') ge 0 or $
	strpos(strlowcase(oper[0]),'min') ge 0 or $
	strpos(strlowcase(oper[0]),'sub') ge 0 or $
	keyword_set(minus) then op='-'
if strpos(oper[0],'/') ge 0 or $
	strpos(strlowcase(oper[0]),'div') ge 0 or $
	keyword_set(divide) then op='/'

;	keywords
xmin=min(X1) > min(X2)
xmax=max(X1) < max(X2)
mbin=0L & if keyword_set(nbin) then mbin=long(nbin[0])
wt1=1. & wt2=1.
if keyword_set(w1) then wt1=w1[0]
if keyword_set(w2) then wt2=w2[0]

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
if not keyword_set(nonorm) then begin
  dc1=dc1/dxout & dc2=dc2/dxout
endif
case op of
  '+': hx=(wt1*dc1)+(wt2*dc2)
  '-': hx=(wt1*dc1)-(wt2*dc2)
  '/': hx=(wt1*dc1)/(wt2*dc2)
  else: hx=(wt1*dc1)*(wt2*dc2)
endcase

return,hx
end
