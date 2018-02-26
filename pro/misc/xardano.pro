function xardano,c0,c1,c2,c3,c4,iroot=iroot,indeg=indeg,eps=eps,$
	verbose=verbose, _extra=e
;+
;function	xardano
;	computes and returns the roots of a polynomial
;		c0 + c1*x + c2*x^2 + c3*x^3 + c4*x^4 = 0
;	of degree < 5
;
;syntax
;	roots=xardano(c0[,c1[,c2[,c3[,c4]]]],/iroot,indeg=indeg,eps=eps,$
;	verbose=verbose)
;
;parameters
;	c0	[INPUT; required] constant coefficient
;	c1	[INPUT] coefficient of linear term
;	c2	[INPUT] coefficient of quadratic term
;	c3	[INPUT] coefficient of cubic term
;	c4	[INPUT] coefficient of quartic term
;		* unless INDEG is set, the largest present and non-zero C
;		  determines the degree of the solution and the number of
;		  elements in the output array
;
;keywords
;	iroot	[INPUT] if set, always outputs complex roots regardless of
;		whether any of the roots are complex.
;		* if even a single root is complex, this is ignored and a
;		  complex array is output
;		* if explicitly set to 0, returns only the real roots
;		* if IROOT=0 and no real roots exist, returns [!values.F_NAN]
;	indeg	[INPUT] if set, ignores all coefficients above C{INDEG}
;		* if it turns out that C{INDEG} is zero, it gets reset to EPS
;	eps	[INPUT; default=1e-6] a small number
;	verbose	[INPUT] controls chatter
;	_extra	[JUNK] here only to prevent crashing the program
;
;description
;	INDEG=0 is here only for completeness.  Obviously, the root is either
;	[0.] (if C0=0) or [!values.F_NAN]
;	INDEG=1 returns -c0/c1
;	INDEG=2 is the familiar (b+-sqrt(b^2-4ac))/2a, with a=c2,b=c1,c=c0
;	INDEG=3, the solution comes from MacSyma
;	INDEG=4, the solution comes from https://en.wikipedia.org/wiki/Quartic_function
;etymology
;	Named after Gerolamo Cardano, because.
;
;history
;	Vinay Kashyap (MMXVI.V)
;-

;	usage
ok='ok' & np=n_params()
n0=n_elements(c0) & n1=n_elements(c1) & n2=n_elements(c2) & n3=n_elements(c3) & n4=n_elements(c4)
if np eq 0 then ok='Insufficient parameters' else $
 if n0 eq 0 then ok='C0 is undefined' else $
  if n0 gt 1 then ok='C0 must be a scalar' else $
   if n1 gt 1 then ok='C1 must be a scalar' else $
    if n2 gt 1 then ok='C2 must be a scalar' else $
     if n3 gt 1 then ok='C3 must be a scalar' else $
      if n4 gt 1 then ok='C4 must be a scalar'
if ok ne 'ok' then begin
  print,'Usage: roots=xardano(c0[,c1[,c2[,c3[,c4]]]],/iroot,indeg=indeg,eps=eps,$'
  print,'       verbose=verbose)'
  print,'  returns the roots of c0 + c1*x + c2*x^2 + c3*x^3 + c4*x^4 = 0'
  if np ne 0 then message,ok,/informational
  return,-1L
endif

;	parse inputs
vv=0L & if keyword_set(verbose) then vv=long(verbose[0])>1L
;
if not keyword_set(eps) then epsilon=1e-6 else epsilon=double(abs(eps[0]))
;
ndeg=n1+n2+n3+n4
case ndeg of
  0: begin & aa=double(c0[0]) & bb=0.D & cc=0.D & dd=0.D & ee=0.D & end	;zeros just to define the variables
  1: begin & aa=double(c1[0]) & bb=double(c0[0]) & cc=0.D & dd=0.D & ee=0.D & end
  2: begin & aa=double(c2[0]) & bb=double(c1[0]) & cc=double(c0[0]) & dd=0.D & ee=0.D & end
  3: begin & aa=double(c3[0]) & bb=double(c2[0]) & cc=double(c1[0]) & dd=double(c0[0]) & ee=0.D & end
  4: begin & aa=double(c4[0]) & bb=double(c3[0]) & cc=double(c2[0]) & dd=double(c1[0]) & ee=double(c0[0]) & end
  else: message,'cannot understand NDEG='+strtrim(ndeg,2)
endcase
;
if keyword_set(indeg) then begin
  mdeg=fix(indeg[0])>0
  if mdeg gt 4 then begin
    message,'cannot solve quintics and higher generally; limiting to quartic',/informational
    mdeg=4
  endif
endif else begin
  mdeg=ndeg
endelse
if aa eq 0 then aa=epsilon

;	degree 0
if mdeg eq 0 then begin
  ;	really, this is just a fake
  if abs(aa-epsilon) le epsilon then return,[0.] else return,[!values.F_NAN]
endif

;	degree 1
if mdeg eq 1 then return,[-bb/aa]

;	degree 2
if mdeg eq 2 then begin
  a2=2.*aa
  b2=bb^2
  ac4=4.*aa*cc
  disc=b2-ac4
  freals=bb/a2
  sqdisc2a=sqrt(abs(disc))/a2
  if disc lt 0 then begin
    if vv gt 0 then message,'complex roots',/informational
    fims=sqdisc
    xr=[complex(freals,fims),complex(freals,-fims)]
  endif else begin
    xr=freals+sqdisc*[-1,1]
    if keyword_set(iroot) then xr=complex(xr)
  endelse
endif

;	degree 3
;	in maxima, do
;	(%i1) a*x^3+b*x^2+c*x+d=0;
;	(%i2) solve(%i1,x);
if mdeg eq 3 then begin
  xr=dcomplexarr(3)
  zz1 = dcomplex(-0.5, -sqrt(3.)/2.)
  zz2 = dcomplex(-0.5,  sqrt(3.)/2.)
  f1 = dcomplex(27.*aa^2*dd^2 + (4.*bb^3-18.*aa*bb*cc)*dd + 4.*aa*cc^3 - bb^2*cc^2)
  f2 = 27.*aa^2*dd - 9.*aa*bb*cc + 2.*bb^3
  f1d= f1/4./(3.^3)/aa^4
  f2d= f2/54./aa^3
  f1z= sqrt(dcomplex(f1d))
  f3 = (f1z-f2d)^(1./3.)
  f4 = (3.*aa*cc - bb^2)
  f5 = 9.*aa^2*f3
  f6 = bb/3./aa

  xr[0] = zz1*f3 - zz2*f4/f5 - f6
  xr[1] = zz2*f3 - zz1*f4/f5 - f6
  xr[2] = f3 - f4/f5 - f6

  if n_elements(iroot) ne 0 then begin
    if iroot[0] eq 0 then begin
      zr=dblarr(3)+!values.F_NAN
      if abs(imaginary(xr[0])) lt epsilon then zr[0]=float(xr[0])
      if abs(imaginary(xr[1])) lt epsilon then zr[1]=float(xr[1])
      if abs(imaginary(xr[2])) lt epsilon then zr[2]=float(xr[2])
      xr=zr
    endif
  endif
endif

;	degree 4
if mdeg eq 4 then begin
  message,'roots for 4th degree polynomial not yet implemented.',/informational
  xr=fltarr(4)+!values.F_NAN
endif

if vv gt 1000 then stop,'halting; type .CON to continue'

return,xr
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;example
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;	cubic
if not keyword_set(c0) then c0=1.
if not keyword_set(c1) then c1=-10.
if not keyword_set(c2) then c2=10.
if not keyword_set(c3) then c3=1.
if not keyword_set(verbose) then verbose=10

print,xardano(c0,c1,c2,c3,verbose=verbose)
xx=findgen(150)*0.1-12. & yy=c0+c1*xx+c2*xx^2+c3*xx^3
plot,xx,yy & oplot,xx,0*xx,line=1

end
