function mk_bknpower,x,norm,gamma1,gamma2,pbreak,pder,Xo=Xo,$
	nobreak=nobreak,verbose=verbose, _extra=e
;+
;function	mk_bknpower
;	returns a broken Power-law function,
;	  p(X)	= NORM*(X/Xo)^gamma1 {X.le.PBREAK}
;		= NORM*(PBREAK/Xo)^(gamma1-gamma2)*(X/Xo)^gamma2 {X.ge.PBREAK}
;
;syntax
;	p=mk_bknpower(x,norm,gamma1,gamma2,pbreak,pder,Xo=Xo,/nobreak,verbose=verbose)
;
;parameters
;	X	[INPUT; required] where p(X) must be computed
;	norm	[INPUT; required] normalization
;	gamma1	[INPUT; required] power-law index 1 (valid for 0<X<PBREAK)
;	gamma2	[INPUT; necessary] power-law index 2 (valid for X>PBREAK)
;	pbreak	[INPUT; necessary] break point
;		* GAMMA2 and PBREAK are ignored if NOBREAK is set
;	pder	[OUTPUT; optional] partial derivatives of model wrt parameters
;		at each X; calculated only if 6 parameters are supplied in call.
;		* array of size [N(X),4], with columns containing the partial
;		  derivatives wrt NORM, GAMMA1, GAMMA2, and PBREAK respectively
;		* to have it return PDER as [N(X),2] for a non-broken power-law
;		  distribution, set NOBREAK (must have GAMMA2 and PBREAK in the
;		  call, but they will be ignored)
;
;keywords
;	Xo	[INPUT] X-value at which which normalization is defined
;		* default is 1.0
;	nobreak	[INPUT] if set, ignores GAMMA2 and PBREAK, essentially assuming
;		that PBREAK=$\infty$
;	verbose	[INPUT] controls chatter
;	_extra	[JUNK] here only to prevent crashing the program
;
;example
;	x=(findgen(100)+1)/10. & p=mk_bknpower(x,1.,-1,-3,1.5) & plot,x,p,/xl,/yl
;	p2=mk_bknpower(x,1.,-1,-3,1.5,/nobreak) & oplot,x,p2,line=2,thick=2
;
;history
;	Vinay Kashyap (Apr2008)
;-

ok='ok' & np=n_params() & nx=n_elements(x)
nn=n_elements(norm) & ng1=n_elements(gamma1)
yesb=1 & if keyword_set(nobreak) then yesb=0
ng2=n_elements(gamma2) & nb=n_elements(pbreak)
if np eq 0 then ok='Insufficient parameters' else $
 if np lt 3 then ok='Insufficient parameters' else $
  if keyword_set(yesb) and np lt 5 then ok='Insufficient parameters' else $
   if nx eq 0 then ok='X is undefined' else $
    if nn eq 0 then ok='NORM is undefined' else $
     if ng1 eq 0 then ok='GAMMA1 is undefined' else $
      if keyword_set(yesb) and ng2 eq 0 then ok='GAMMA2 is undefined' else $
       if keyword_set(yesb) and nb eq 0 then ok='PBREAK is undefined' else $
        if nn gt 1 then ok='NORM must be scalar' else $
	 if ng1 gt 1 then ok='GAMMA1 must be scalar' else $
	  if keyword_set(yesb) and ng2 gt 1 then ok='GAMMA2 must be scalar' else $
	   if keyword_set(yesb) and nb gt 1 then ok='PBREAK must be scalar'
if ok ne 'ok' then begin
  print, 'Usage: p=mk_bknpower(x,norm,gamma1,gamma2,pbreak,pder,Xo=Xo)
  print, '  generates a broken power-law p(X)'
  if np ne 0 then message,ok,/informational
  return,-1L
endif

;	check inputs
xx=x
y0=norm[0] & g1=gamma1[0]
g2=1.0 & xbrk=max([xx])+1.
if keyword_set(yesb) then begin
  g2=gamma2[0] & xbrk=pbreak[0]
endif

;	keywords
vv=0L & if keyword_set(verbose) then vv=long(verbose[0])>1
x0=1.0 & if keyword_set(Xo) then x0=Xo[0]*1.0

pp=fltarr(nx)
o1=where([xx] le xbrk,mo1,complement=o2,ncomplement=mo2)

;	do not multiply by NORM yet, in case it happens to be 0
if mo1 gt 0 then pp[o1]=(xx[o1]/x0)^(g1)
if mo2 gt 0 then pp[o2]=((xbrk/x0)^(g1-g2))*((xx[o2]/x0)^g2)

;	catch any NaNs and infinities and set them to 0
oy=where(finite(pp) eq 0,moy)
if moy gt 0 then pp[oy]=0.

;compute partial dervatives
;	PDER Order:  wrt norm, gamma1, gamma2, pbreak
if np ge 6 then begin
  if keyword_set(yesb) then pder=fltarr(nx,4) else pder=fltarr(nx,2)

  ;NORM
  pder[*,0]=pp

  ;GAMMA1
  if mo1 gt 0 then pder[o1,1]=y0*pp[o1]*alog(xx[o1]/x0)
  if mo2 gt 0 then pder[o2,1]=y0*pp[o2]*alog(xbrk/x0)

  if keyword_set(yesb) then begin

    ;GAMMA2
    if mo1 gt 0 then pder[o1,2]=0.
    if mo2 gt 0 then pder[o2,2]=-y0*pp[o2]*alog(xbrk/x0)+y0*pp[o2]*alog(xx[o2]/x0)

    ;PBREAK
    if mo1 gt 0 then begin
      dbrk=0.01*xbrk & if dbrk eq 0 then dbrk=0.001
      pp1=mk_bknpower(xx[o1],y0,g1,g2,xbrk+dbrk)
      pp2=mk_bknpower(xx[o1],y0,g1,g2,xbrk-dbrk)
      pder[o1,3]=(pp1-pp2)/(dbrk*2.)
    endif
    if mo2 gt 0 then pder[o2,3]=y0*pp[o2]*(g1-g2)*(x0/xbrk)

  endif
endif

return,y0*pp	;multiply by NORM and exit
end
