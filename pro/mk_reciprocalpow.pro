function mk_reciprocalpow,x,norm,gamma1,gamma2,pjoin,pder,verbose=verbose,$
	_extra=e
;+
;function	mk_reciprocalpow
;	returns a reciprocal double Power-law function,
;	  p(X)	= 1/(A1*X^GAMMA1 + A2*X^GAMMA2),
;	  with A1*J^GAMMA1 = A2*J^GAMMA2, so A2 == A1*J^(GAMMA1-GAMMA2)
;	  and A1 == NORM
;	  for X>0
;
;syntax
;	p=mk_reciprocalpow(x,norm,gamma1,gamma2,pjoin,pder,verbose=verbose)
;
;parameters
;	X	[INPUT; required] where p(X) must be computed
;	norm	[INPUT; required] normalization
;	gamma1	[INPUT; required] power-law index 1
;	gamma2	[INPUT; required] power-law index 2
;	pjoin	[INPUT; required] the joining point where the two
;		power-law components are equal
;	pder	[OUTPUT; optional] partial derivatives of model wrt parameters
;		at each X; calculated only if 6 parameters are supplied in call.
;		* array of size [N(X),4], with columns containing the partial
;		  derivatives wrt NORM, GAMMA1, GAMMA2, and PJOIN respectively
;
;keywords
;	verbose	[INPUT] controls chatter
;	_extra	[JUNK] here only to prevent crashing the program
;
;example
;	x=(findgen(100)+1)/10. & p=mk_reciprocalpow(x,1.,1.5,0.2,2) & plot,x,p,/xl,/yl
;	contrast with: oplot,x,mk_bknpower(x,0.3,-1/1.5,-5,2),line=2
;
;history
;	Vinay Kashyap (Sep2008)
;-

ok='ok' & np=n_params() & nx=n_elements(x)
nn=n_elements(norm) & ng1=n_elements(gamma1)
ng2=n_elements(gamma2) & nb=n_elements(pjoin)
if np eq 0 then ok='Insufficient parameters' else $
 if np lt 5 then ok='Not enough parameters' else $
  if keyword_set(yesb) and np lt 5 then ok='Insufficient parameters' else $
   if nx eq 0 then ok='X is undefined' else $
    if nn eq 0 then ok='NORM is undefined' else $
     if ng1 eq 0 then ok='GAMMA1 is undefined' else $
      if ng2 eq 0 then ok='GAMMA2 is undefined' else $
       if nb eq 0 then ok='PJOIN is undefined' else $
        if nn gt 1 then ok='NORM must be scalar' else $
	 if ng1 gt 1 then ok='GAMMA1 must be scalar' else $
	  if ng2 gt 1 then ok='GAMMA2 must be scalar' else $
	   if nb gt 1 then ok='PJOIN must be scalar'
if ok ne 'ok' then begin
  print, 'Usage: p=mk_reciprocalpow(x,norm,gamma1,gamma2,pjoin,pder,Xo=Xo)
  print, '  generates a reciprocal double power-law p(X)'
  if np ne 0 then message,ok,/informational
  return,-1L
endif

;	check inputs
xx=x
y0=norm[0] & g1=gamma1[0] & g2=gamma2[0] & x0=pjoin[0]
a1=y0
a2=a1 & if x0 gt 0 then a2=exp(alog(abs(a1))+(g1-g2)*alog(x0))

;	keywords
vv=0L & if keyword_set(verbose) then vv=long(verbose[0])>1

;	do not multiply by NORM yet, in case it happens to be 0
op=where(xx gt 0,mop) & pp=0.*xx & p1=pp & p2=pp
if mop gt 0 then p1[op]=exp(g1*alog(xx[op]))
if mop gt 0 then p2[op]=exp(g2*alog(xx[op]))*x0^(g1-g2)
o0=where(p1+p2 ne 0,mo0) & if mo0 gt 0 then pp[o0]=1./(p1[o0]+p2[o0])

;	catch any NaNs and infinities and set them to 0
oy=where(finite(pp) eq 0,moy)
if moy gt 0 then pp[oy]=0.

;compute partial dervatives
;	PDER Order:  wrt norm, gamma1, gamma2, pjoin
if np ge 6 then begin
  pder=fltarr(nx,4)

  ;NORM
  if a1 ne 0 then pder[*,0]=-pp/a1^2

  ;GAMMA1
  if a1 ne 0 and x0 gt 0 then pder[op,1]=-pp^2 * ( (xx[op]^g1) * alog(xx[op]) + (xx[op]^g2) * (x0^(g1-g2)) * alog(x0) ) / a1^2

  ;GAMMA2
  if a1 ne 0 and x0 gt 0 then pder[op,2]=-pp^2 * ( (xx[op]^g2) * (x0^(g1-g2)) * alog(xx[op]/x0) ) / a1^2

  ;PJOIN
  if a1 ne 0 and x0 gt 0 then pder[op,3]=-pp^2 * (g1-g2) * (xx[op]^g2) * (x0^(g1-g2-1)) / a1

endif

if a1 ne 0 then fp=pp/a1 else fp=0*pp	;multiply by NORM and exit

return,fp
end
