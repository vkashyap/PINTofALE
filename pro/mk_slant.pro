function mk_slant,x,mean,core,peak,pder,angle=angle,betap=betap,$
	normflx=normflx,missing=missing, _extra=e
;+
;function	mk_slant
;	returns the product of a slanted line with a modified
;	Lorentzian, Lb(X)
;
;	the idea is to model the asymmetrical line profiles seen in
;	CXO grating data.
;
;syntax
;	Lbs=mk_slant(x,mean,core,peak,pder,angle=angle,betap=betap,$
;	/normflx,missing=missing)
;
;parameters
;	x	[INPUT array; required] where G(X) must be computed
;	mean	[INPUT; default: mid_point(X)] position of peak
;	core	[INPUT; default: 0.1*range(X)] core width
;	peak	[INPUT; default: 1] value of Lb(X=MEAN)
;	pder	[OUTPUT; optional] partial derivatives of model wrt parameters
;		at each X; calculated only if 5 parameters are supplied in call.
;		* array of size [N(X),5], with columns containing the partial
;		  derivatives wrt MEAN, CORE, PEAK, BETAP, and ANGLE respectively
;
;keywords
;	angle	[INPUT; default=0] the atan(slope) of the straight line
;		that makes the Lorentzian asymmetric
;	betap	[INPUT; default=1] the index of the beta-profile.  default
;		is regular Lorentzian.
;		* if NORMFLX is set and BETAP.le.0.5, the default is used
;	normflx	[INPUT] if set, {\int_{-\infty}^{\infty} dX Lbs(X) = PEAK}
;	missing	[INPUT] 3 element array to populate missing values of
;		MEAN, CORE, and PEAK
;	_extra	[JUNK] here only to prevent crashing
;
;description
;	The Lorentzian is
;		L(X) = PEAK / ( 1 + ((X-MEAN)/CORE)^2 )
;	The modified Lorentzian is
;		Lb(X) = PEAK / ( 1 + ((X-MEAN)/CORE)^2 ) ^ BETAP
;	The asymmetric modified Lorentzian is
;		Lbs(X) = (PEAK + PEAK*SLOPE*(X-MEAN)) / $
;			 ( 1 + ((X-MEAN)/CORE)^2 ) ^ BETAP
;	
;	When integrated over the real axis (Gradshteyn & Ryzhik, 3.251,2),
;	the second terms drops out and the integral is identical to that
;	of the modified Lorentzian itself, i.e.,
;		\int dX Lbs = CORE*PEAK*B(1/2,BETAP-1/2)
;	where B(x,y) is the Beta-function, with BETAP > 1/2
;	
;	if NORMFLX is set,
;		Lbs(X) = ((PEAK+PEAK*SLOPE*(X-MEAN))/CORE/B(1/2,beta-1/2)) / $
;			 (1+((X-MEAN)/CORE)^2)^BETAP
;
;usage summary
;	* call as a function
;	* generates modified Lorentzian model only at specified points X
;	* needs MEAN, CORE, PEAK for complete specification
;
;subroutines
;	NONE
;
;history
;	vinay kashyap (MarMM, based on MK_LORENTZ)
;	now works even if X are integers (VK; Jul01)
;	converted array notation to IDL 5 (VK; Apr02)
;	accounted for case PEAK=0, added partial derivatives of BETAP and
;	  ANGLE to PDER (VK; Jun02)
;	changed keyword NORM to NORMFLX (VK; Oct02)
;-

np=n_params()
if np lt 1 then begin
  print, 'Usage: Lbs=mk_slant(x,mean,core,peak,pder,angle=a,betap=b,missing=m,/normflx)'
  print, '  generates an asymmetric modified Lorentzian Lbs(x;beta,angle)'
  return,[-1L]
endif

;initialize
nx=n_elements(x) & x0=x[nx/2] & mxx=max(x,min=mnx)

;figure out the defaults
if not keyword_set(betap) then b=1. else b=float(betap[0])
if not keyword_set(angle) then ang=0. else ang=float(angle[0])
while ang gt 90 do ang=ang-90. & if ang eq 90. then ang=ang-1e-3
while ang lt -90 do ang=ang+90. & if ang eq -90. then ang=ang+1e-3
slope=tan(ang*!pi/180.)
if keyword_set(normflx) and b lt 0.5 then begin
  message,'Normalization becomes infinite!  Resetting BETAP to 0.5',/info
  b=0.5
endif
if not keyword_set(missing) then missing=[x0,0.1*(mxx-mnx),1.]
if np lt 4 then p=missing[2] else p=peak[0]
if np lt 3 then c=missing[1] else c=core[0]
if c lt 0 then c=abs(c) & if c eq 0 then c=missing[1]
if np lt 2 then m=missing[0] else m=mean[0]

;	renorm
if keyword_set(normflx) then begin
  if b gt 0.5 then bnorm=1./beta(0.5,b-0.5)/c else bnorm=0.
endif else bnorm=1.
p=p*bnorm

;	compute function
z=(x-m+0.0)/c & z=alog10(1.+z^2) & z=-b*z
Lbs=make_array(size=size(0.0*x))
oz=where(z gt -30,moz) & if moz gt 0 then Lbs[oz]=10.^(z[oz])
Lbs=Lbs*(1.+slope*(x-m)) > 0.

;	compute partial derivatives
if np ge 5 then begin
  pder = fltarr(nx,5)
  z=(x-m+0.0)/c & zz=(1.+z^2) & pp=p*(1.+slope*(x-m+0.0))
  ;	partial wrt MEAN
  Lbsm=pp*(zz^(-b-1))*(2.*b/c^2)*(x-m+0.0)-(zz^(-b))*p*slope
  pder[*,0] = Lbsm[*]
  ;	partial wrt CORE
  Lbsc=pp*(zz^(-b-1))*2.*b*(x-m+0.0)^2/c^3
  if keyword_set(normflx) then Lbsc=Lbsc-Lbs/c
  pder[*,1] = Lbsc[*]
  ;	partial wrt PEAK
  Lbsp=Lbs*bnorm
  pder[*,2] = Lbsp[*]
  ;	partial wrt BETAP
  if not keyword_set(normflx) then begin
    Lbsb=Lbs*p & oo=where(Lbsb gt 0,moo)
    if moo gt 0 then pder[oo,3]=-Lbsb[oo]*alog(Lbsb[oo])
  endif else begin
    ;	because this is easier than differentiating the beta function
    delb=0.01 & b2=b+delb & b1=(b-delb)>0.5
    tmp1=mk_slant(x,m,c,p,angle=ang,betap=b1,/normflx)
    tmp2=mk_slant(x,m,c,p,angle=ang,betap=b2,/normflx)
    pder[*,3]=(tmp2-tmp1)/delb/2.
  endelse
  ;	partial wrt ANGLE
  dela=0.01 & a1=ang-dela & a2=ang+dela
  tmp1=mk_slant(x,m,c,p,angle=a1,betap=b,normflx=normflx)
  tmp2=mk_slant(x,m,c,p,angle=a2,betap=b,normflx=normflx)
  pder[*,4]=(tmp2-tmp1)/dela/2.
endif
Lbs=Lbs*p

return,Lbs
end
