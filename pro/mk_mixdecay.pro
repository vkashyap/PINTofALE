function mk_mixdecay,x,x0,tau0,y0,pder,ybkg=ybkg,verbose=verbose, _extra=e
;+
;function	mk_mixdecay
;	returns a mixture of decays, \Sum_{i=0}^{N} y_i*exp(-(x-x_i)/tau_i) + ybkg
;
;syntax
;	py=mk_mixdecay(x,x0,tau0,y0,pder,ybkg=ybkg,verbose=verbose)
;
;parameters
;	x	[INPUT array; required] where G(X) must be computed
;	x0	[INPUT; required] the locations of the peaks of each component
;		* the component is assumed to be identically 0 for x<x0
;	tau0	[INPUT; required] the e-folding decay scales of each component
;	y0	[INPUT; required] the peak intensity of each component
;	pder	[OUTPUT; optional] partial derivatives of model wrt parameters
;		at each X; calculated only if 4 parameters are supplied in call.
;
;keywords
;	ybkg	[INPUT; default=0] background pedestal level
;	verbose	[INPUT; default=0] controls chatter
;	_extra	[JUNK] here only to prevent crashing
;
;usage summary
;	* call as a function
;
;subroutines
;	NONE
;
;history
;	vinay kashyap (Apr17)
;-

;	usage
ok='ok' & np=n_params() & nx=n_elements(x) & n1=n_elements(x0) & n2=n_elements(tau0) & n3=n_elements(y0)
if np lt 2 then ok='Insufficient parameters' else $
 if nx eq 0 then ok='X is not defined' else $
  if n1 eq 0 then ok='X0 not defined' else $
   if n2 eq 0 then ok='TAU0 not defined' else $
    if n3 eq 0 then ok='Y0 not defined' else $
     if n1 ne n2 then ok='X0 and TAU0 not compatible' else $
      if n1 ne n3 then ok='X0 and Y0 not compatible'
if ok ne 'ok' then begin
  print,'Usage: py=mk_mixdecay(x,x0,tau0,y0,pder,ybkg=ybkg,verbose=verbose)'
  print,'  returns a mixture of decays, \Sum_{i=0}^{N} y_i*exp(-(x-x_i)/tau_i) + ybkg'
  print,'  for aa=[[x_i,tau_i,y_i],i=1..N] in mk_3model'
  if np ne 0 then message,ok,/info
  return,-1L
endif

;	inputs
yz=fltarr(nx) & nb=n_elements(ybkg)
if nb gt 0 then begin
  if nb eq 1 then yz[*]=ybkg[0] else yz[0:(nb<nx)-1L]=ybkg[0:(nb<nx)-1L]
endif
vv=0 & if keyword_set(verbose) then vv=long(verbose[0])>1

;	compute function and partial derivatives
py=yz
if np ge 3 then pder=fltarr(nx,3*n1)
for i=0L,n1-1L do begin
  oi=where(x ge x0[i],moi)
  tmp=0.D*x & pyi=tmp
  if moi gt 0 then begin
    if tau0[i] ne 0 then tmp[oi]=(x[oi]-x0[i])/tau0[i] else tmp=0*x
    pyi[oi]=y0[i]*exp(-tmp[oi])
    if np ge 3 then begin
      ;	partial wrt X0 : y/tau0
      if tau0[i] ne 0 then pder[oi,3*i]=pyi[oi]/tau0[i]
      ;	partial wrt tau0 : [(x-x0)/tau^2]*y
      if tau0[i] ne 0 then pder[oi,3*i+1]=pyi[oi]*(x[oi]-x0[i])/tau0[i]^2
      ;	partial wrt y0 : y/y0
      pder[oi,3*i+2]=exp(-tmp[oi])
    endif
  endif
  py=py+pyi
endfor

if vv gt 50 then begin
  plot,x,py,xtitle='X',ytitle='P!d'+strtrim(n1,2)+'!n(X)'
  oplot,x,yz,line=2
  if vv gt 1000 then stop,'type .CON to continue'
endif

return,Py
end

;	example
x=findgen(1000)*0.1 & x0=[-10.,30.,40.,70.,101.] & tau0=[100.,20.,10.,50.,20.] & y0=[1.,2.,3.,4.,5.] & ybkg=0.*x+1

;	calling sequence
print,mk_mixdecay()

print,'' & print,'py=mk_mixdecay(x,x0,tau0,y0,pder,ybkg=ybkg,verbose=51)'
py=mk_mixdecay(x,x0,tau0,y0,pder,ybkg=ybkg,verbose=1001)

end
