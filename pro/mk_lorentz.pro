function mk_lorentz,x,mean,core,peak,pder,betap=betap,normflx=normflx,$
	missing=missing,fwhm=fwhm, _extra=e
;+
;function	mk_lorentz
;	returns a modified Lorentzian, Lb(X)
;
;syntax
;	l=mk_lorentz(x,mean,core,peak,pder,/betap,/normflx,missing=m,/fwhm)
;
;parameters
;	x	[INPUT array; required] where G(X) must be computed
;	mean	[INPUT; default: mid_point(X)] position of peak
;	core	[INPUT; default: 0.1*range(X)] core width
;	peak	[INPUT; default: 1] value of Lb(X=MEAN)
;	pder	[OUTPUT; optional] partial derivatives of model wrt
;		parameters at each X; calculated only if 5 parameters
;		are supplied in call.
;		* array of size [N(X),4], with columns containing the partial
;		  derivatives wrt MEAN, CORE, PEAK, and BETAP respectively
;
;keywords
;	betap	[INPUT; default=1] the index of the beta-profile.  default
;		is regular Lorentzian.
;		* if NORMFLX is set and BETAP.le.0.5, the default is used
;	normflx	[INPUT] if set, {\int_{-\infty}^{\infty} dX Lb(X) = PEAK}
;	missing	[INPUT] 3 element array to populate missing values of
;		MEAN, CORE, and PEAK
;	fwhm	[INPUT] if set, assumes that the input CORE is actually
;		given as full-width at half-max, and converts to the true
;		core width prior to calculation, and then converts back to
;		full-width at half-max post-calc.
;	_extra	[JUNK] here only to prevent crashing
;
;description
;	The Lorentzian is
;		L(X) = PEAK / ( 1 + ((X-MEAN)/CORE)^2 )
;	The modified Lorentzian is
;		Lb(X) = PEAK / ( 1 + ((X-MEAN)/CORE)^2 ) ^ BETAP
;
;	When integrated over the real axis (Gradshteyn & Ryzhik, 3.251,2),
;		\int dX Lb = PEAK*CORE*B(1/2, BETAP-1/2),
;	where B(x,y) is the Beta-function, with BETAP > 1/2
;
;	hence, if NORMFLX is set,
;		Lb(X) = (PEAK/CORE/B(1/2,beta-1/2))/(1+((X-MEAN)/CORE)^2)^BETAP
;
;	The hwhm is defined as that value of X=W such that
;		Lb(W)=0.5*Lb(X=MEAN)
;	i.e.,	Lb(W)=(stuff)/(1+((W-MEAN)/CORE)^2)^(BETAP)=0.5*(stuff)
;	hence,	hwhm = W-MEAN = CORE*sqrt(2^(1/BETAP)-1)
;	or,	fwhm = CORE*2*sqrt(2^(1/BETAP)-1)
;
;usage summary
;	* call as a function
;	* generates modified Lorentzian model only at specified points X
;	* needs MEAN, CORE, PEAK for complete specification
;
;examples
;	x=findgen(1000)*0.01 & peasecolr & b=2.5
;	L1=mk_lorentz(x,5,1,1) & L2=mk_lorentz(x,5,1,1,/fwhm)
;	Lb1=mk_lorentz(x,5,1,1,betap=b) & Lb2=mk_lorentz(x,5,1,1,betap=b,/fwhm)
;	Lf1=mk_lorentz(x,5,1,1,/normflx) & Lf2=mk_lorentz(x,5,1,1,betap=b,/norm)
;	plot,x,L1 & oplot,x,L2,col=2
;	oplot,x,Lb1,thick=2 & oplot,x,Lb2,thick=2,col=2
;	oplot,x,Lf1,thick=2,line=1 & oplot,x,Lf2,thick=2,col=3,line=1
;
;subroutines
;	NONE
;
;history
;	vinay kashyap (Oct98)
;	changed keyword BETA to BETAP; corrected normalization (VK; Dec98)
;	recomputed partial derivatives (VK; MarMM)
;	now works even if X are integers (VK; Jul01)
;	what if PEAK=0? (VK; Apr02)
;	converted array notation to IDL 5 (VK; Apr02)
;	added partial derivative wrt BETAP to PDER output array (VK; Jun02)
;	changed keyword NORM to NORMFLX (VK; Oct02)
;	added keyword FWHM (VK; Apr03)
;-

np=n_params()
if np lt 1 then begin
  print, 'Usage: Lb=mk_lorentz(x,mean,core,peak,pder,betap=b,missing=m,/normflx,/fwhm)'
  print, '  generates a modified Lorentzian Lb(x;beta)'
  return,[-1L]
endif

;initialize
nx=n_elements(x) & x0=x[nx/2] & mxx=max(x,min=mnx)

;	figure out the defaults
if not keyword_set(betap) then b=1. else b=betap[0]
if keyword_set(normflx) and b lt 0.5 then begin
  message,'Normalization becomes infinite!  Resetting BETAP to 0.5',/info
  b=0.5
endif
if not keyword_set(missing) then missing=[x0,0.1*(mxx-mnx),1.]
if np lt 4 then p=missing[2] else p=peak[0]
if np lt 3 then c=missing[1] else c=core[0]
if c lt 0 then c=abs(c) & if c eq 0 then c=missing[1]
if np lt 2 then m=missing[0] else m=mean[0]

;	fwhm?
if keyword_set(fwhm) then c2f=2.*sqrt(2.^(1./b)-1.) else c2f=1.
c=c/c2f		;if input is FWHM, convert to CORE

;	renorm
if keyword_set(normflx) then begin
  if b gt 0.5 then bnorm=1./beta(0.5,b-0.5)/c else bnorm=0.
endif else bnorm=1.
p=p*bnorm

;	compute function
z=(x-m+0.0)/c & z=alog10(1.+z^2) & z=-b*z
Lb=make_array(size=size(0.*x))
oz=where(z gt -30,moz) & if moz gt 0 then Lb[oz]=10.^(z[oz])

;	compute partial derivatives
if np ge 5 then begin
  pder = fltarr(nx,4)
  z=(x-m+0.0)/c & zz=(1.+z^2)
  ;	partial wrt MEAN
  Lbm=p*(zz^(-b-1))*(2.*b/c^2)*(x-m+0.0)
  pder[*,0] = Lbm[*]
  ;	partial wrt CORE
  Lbc=p*(zz^(-b-1))*2.*b*(x-m+0.0)^2/c^3
  if keyword_set(normflx) then Lbc=Lbc-p*Lb/c
  pder[*,1] = Lbc[*]
  ;	partial wrt PEAK
  Lbp=Lb*bnorm
  pder[*,2]=Lbp[*]
  ;	partial wrt BETAP
  if not keyword_set(normflx) then begin
    Lbb=Lb*p & oo=where(Lbb gt 0,moo)
    if moo gt 0 then pder[oo,3]=-Lbb[oo]*alog(Lbb[oo])
  endif else begin
    ;	because this is easier than differentiating the beta function
    delb=0.01 & b2=b+delb & b1=(b-delb)>0.5
    tmp1=mk_lorentz(x,m,c,p,betap=b1,/normflx)
    tmp2=mk_lorentz(x,m,c,p,betap=b2,/normflx)
    pder[*,3]=(tmp2-tmp1)/delb/2.
  endelse
endif
Lb=Lb*p

;	fwhm?
c=c*c2f		;if input was FWHM, convert back from CORE

;if np ge 5 then begin
;  pder = fltarr(nx,3)
;  z=(x-m)/c
;  ;z=(1.+z^2) & z=-p*(b+1)/z^(b+1) & oz=where(abs(z) gt 1e-10,moz)
;  z2=(1.+z^2) & z3=p*(-b)/z2^(b+1) & oz=where(abs(z) gt 1e-10,moz)
;  ;	partial wrt MEAN
;  Lbm=0*Lb & if moz gt 0 then Lbm(oz)=z3(oz)*(1./c^2)*(2.*(x-m)*(-1))
;  pder(*,0) = Lbm(*)
;  ;	partial wrt CORE
;  Lbc=0*Lb & if moz gt 0 then Lbc(oz)=z(oz)*(x-m)^2*(-2./c^3)
;  if keyword_set(normflx) then $
;	if moz gt 0 then Lbc(oz)=Lbc(oz)-Lb(oz)/c
;  pder(*,1) = Lbc(*)
;  ;	partial wrt PEAK
;  Lbp=Lb*bnorm/p
;  pder(*,2) = Lbp(*)
;endif

return,Lb
end
