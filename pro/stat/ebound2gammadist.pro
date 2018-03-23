function ebound2gammadist,xdn,xup,plev=plev,ngrid=ngrid,niter=niter,eps=eps,verbose=verbose, _extra=e
;+
;function	ebound2gammadist
;	compute and return the gamma distribution parameters [alpha,beta]
;	that fit to an asymmetric error bar, [XDN,XUP]
;
;	Note that XDN, XUP must all be >0
;
;syntax
;	ab=ebound2gammadist(xdn,xup,plev=plev,ngrid=ngrid,niter=niter,eps=eps,verbose=verbose)
;
;parameters
;	xdn	[INPUT; required] lower error bound
;	xup	[INPUT; required] upper error bound
;
;keywords
;	plev	[INPUT; default=0.68] equal-tail confidence level
;		for which XDN and XUP are defined
;		* if <0, ignores -ve sign => ABS(PLEV)
;		* if .GE.1, assumes input is in percentage => PLEV/100.
;		* if .GE.100, assumes input is reciprocal => 1/PLEV
;	ngrid	[INPUT; default=1000] number of grid points to use for numerical
;		representation of gamma distribution
;	niter	[INPUT; default=51] number of iterations to use to
;		bracket for solution to alpha -- if iterations exceed
;		this number, kicks out and looks for the intersection
;	eps	[INPUT; default=1d-5] a small number
;	verbose	[INPUT] controls chatter
;	_extra	[JUNK] here only to prevent crashing the program
;
;example usage
;	.run ebound2gammadist
;
;history
;	Vinay Kashyap (2018-mar-23; based on Cook, J.D. (2010), Determining
;	  distribution parameters from quantiles, via Luis Campos)
;-

;	usage
ok='ok' & np=n_params() & nd=n_elements(xdn) & nu=n_elements(xup)
if np lt 2 then ok='Insufficient parameters' else $
 if nd eq 0 then ok='XDN are not defined' else $
  if nu eq 0 then ok='XUP are not defined' else $
   if nd ne nu then ok='XDN and XUP are incompatible'
if ok ne 'ok' then begin
  print,'Usage: ab=ebound2gammadist(xdn,xup,plev=plev,ngrid=ngrid,niter=niter,eps=eps,verbose=verbose)'
  print,'  compute and return [alpha,beta] for gamma distribution'
  print,'  that matches given asymmetrical error bounds'
  if np ne 0 then message,ok,/informational
  return,-1L
endif

;	parameters
zdn=[xdn[*]] & zup=[xup[*]]
oo=where(zdn lt 0,moo) & if moo gt 0 then begin & message,'XDN cannot have -ve values; returning',/informational & return,-1L & endif
oo=where(zup lt 0,moo) & if moo gt 0 then begin & message,'XUP cannot have -ve values; returning',/informational & return,-1L & endif

;	keywords
clev=fltarr(nu)+0.68 & nc=n_elements(plev)
cplev=clev & if nc eq 1 then cplev[*]=plev[0] & if nc gt 1 then cplev[0L:(nc<nu)-1L]=plev[0L:(nc<nu)-1L]
for i=0L,nu-1L do begin
  pp=abs(cplev[i])
  if pp ge 1 and pp lt 100 then pp=pp/100.
  if pp ge 100 then pp=1./pp
  clev[i]=pp
endfor
pdn=0.5-clev/2.
pup=0.5+clev/2.
;
mgrid=1000L
if keyword_set(ngrid) then mgrid=long(ngrid[0])>10L
;
miter=50L
if keyword_set(niter) then miter=long(niter[0])>10L
;
deps=1d-5
if keyword_set(eps) then deps=double(abs(eps[0]))>(1d-30)
;
vv=0L & if keyword_set(verbose) then vv=long(abs(verbose[0]))>1L

;	first find the value of alpha such that F^-1(pup;alpha,1)/F^-1(pdn;alpha,1) = xup/xdn
;	then beta = xup/F^-1(pup;alpha,1) = xdn/F^-1(pdn;alpha,1)

ab=fltarr(2,nu)	;array to hold the output

for i=0L,nu-1L do begin
  x1=xdn[i] & x2=xup[i] & p1=pdn[i] & p2=pup[i]
  if x2 lt x1 then begin & x1=xup[i] & x2=xdn[i] & endif	;catch user error
  rx=x2/x1
  a0=(x1+x2)/2. & b0=(x2-x1)/2.	;starting points
  xmin=(x1-10*b0)>0 & xmax=(x2+10*b0) & dx=(xmax-xmin)/mgrid & xx=findgen(mgrid+1L)*dx+xmin & if xmin eq 0 then xx=xx[1:*]	;grid for Gamma distribution

  ;	solve for alpha
  go_on=1 & kiter=0L & kmax=200L & za=fltarr(kmax+1L) & zr=za & zf=za
  while go_on do begin
    gg=exp((a0-1.)*alog(xx)-xx-lngamma(a0)) & cg=total(gg*dx,/cumul) & if max(cg) gt 1 then cg=cg/max(cg)
    ff=interpol(xx,cg,[p1,p2]) & rr=ff[1]/ff[0] 
    za[kiter]=a0 & zr[kiter]=rr
    if abs(rr-rx) lt deps then go_on=0 else begin
      fac=(((rr/rx)+(rr-rx))>0.9)<1.1
      a0=fac*a0
      zf[kiter]=fac
    endelse
    if kiter ge kmax then go_on=0
    kiter=kiter+1L
    if vv gt 100 then begin & print,kiter,rx,rr,a0,rx-rr,rr/rx & wait,0.01 & endif
  endwhile
  if kiter eq kmax+1L then begin
    os=sort(za)
    if vv gt 100 and !d.name eq 'X' then begin
      plot,za[os],zr[os] & oplot,za[os],0*zr+rx
    endif
    a0=interpol(za[os],zr[os],rx)
  endif
  ab[0,i]=a0

  ;	solve for beta
  gg=exp((a0-1.)*alog(xx)-xx-lngamma(a0)) & cg=total(gg*dx,/cumul) & if max(cg) gt 1 then cg=cg/max(cg)
  ff=interpol(xx,cg,[p1,p2]) & b01=x1/ff[0] & b02=x2/ff[1] & b0=(b01+b02)/2.
  if vv gt 100 then print,b01,b02,b0,a0
  ab[1,i]=b0
endfor

if vv gt 1000 then stop,'HALTing; type .CON to continue'

return,ab
end
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;	an example, using the upper and lower bounds based on Gehrels 1986 (ApJ 303, 336)
nn=[5,10,30]
if not keyword_set(ss) then ss=3.
xdn=nn-ss*sqrt(nn-0.25)+(ss^2-1.)/4.	;Gehrels 1986 Eqn 7
xup=nn+ss*sqrt(nn+0.75)+(ss^2+3.)/4.	;Gehrels 1986 Eqn 11

;	interpolate for plev from Table 3
sarr=[1.,	1.282,	1.645,	1.960,	2.,	2.326,	2.576,	3.,	3.090,	3.291]
parr=[0.8413,	0.9,	0.95,	0.975,	0.9772,	0.99,	0.995,	0.9987,	0.999,	0.9995]
plev=interpol(parr,sarr,ss)

if not keyword_set(ngrid) then ngrid=1000L
if not keyword_set(niter) then niter=200
if not keyword_set(eps) then eps=1e-5
if not keyword_set(verbose) then verbose=0

;	calling sequence
print,''
jnk=ebound2gammadist()
print,''

;	example based on small counts
ab=ebound2gammadist(xdn,xup,plev=plev,ngrid=ngrid,niter=niter,eps=eps,verbose=verbose)
for i=0L,n_elements(nn)-1L do print,$
'N='+strtrim(nn[i],2)+' ['+strtrim(xdn[i],2)+','+strtrim(xup[i],2)+'] @plev='+strtrim(plev,2)+' ==> gamma('+strtrim(ab[0,i],2)+','+strtrim(ab[1,i],2)+')'

end
