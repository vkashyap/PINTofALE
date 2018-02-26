function mk_rogauss,x,mean,sig,peak,pder,vrot=vrot, _extra=e
;+
;function	mk_rogauss
;	convolves a gaussian line profile with a rotational broadening function
;	and returns the broadened profile.
;
;syntax
;	rg=mk_rogauss(x,mean,sig,peak,pder,vrot=vrot,/fwhm,/normflx,$
;	missing=missing)
;
;parameters
;	x	[INPUT array; required] where rG(X) must be computed
;		* generally in a wavelength grid
;	mean	[INPUT] passed w/o comment to MK_GAUSS
;	sig	[INPUT] passed w/o comment to MK_GAUSS
;	peak	[INPUT] passed w/o comment to MK_GAUSS
;	pder	[OUTPUT; optional] partial derivatives of model wrt
;		parameters at each X
;		* calculated only if 5 parameters are supplied in call.
;		* numerically computed because of intractable convolution
;		* array of size [N(X),4], with columns containing the partial
;		  derivatives wrt MEAN, SIG, PEAK, and VROT respectively
;
;keywords
;	vrot	[INPUT; default=0] the rotational velocity, as a fraction
;		of light speed
;		* if > 1 and < 3e5, assumed to be in [km/s]
;		* if > 3e5 and < 3e10, assumed to be in [cm/s]
;		* if -ve, then convolves with a fixed rotation profile
;		  whose width is defined with reference to MEAN
;	_extra	[JUNK] pass defined keywords to MK_GAUSS
;		FWHM, NORMFLX, MISSING
;
;description
;	NORMFLX=0 && FWHM=0:
;	  G(X)=PEAK*exp((X-MEAN)^2/SIGMA^2)
;	NORMFLX=1 && FWHM=0:
;	  G(X)=(PEAK/SIGMA/sqrt(2*!PI))*exp((X-MEAN)^2/2/SIGMA^2)
;	NORMFLX=0 && FWHM=1:
;	  G(X)=PEAK*exp(2.355^2*(X-MEAN)^2/2/SIGMA^2)
;	NORMFLX=1 && FWHM=1:
;	  G(X)=(2.355*PEAK/SIGMA/sqrt(2*!PI))*exp(2.355^2*(X-MEAN)^2/2/SIGMA^2)
;
;	G(X) is then convolved with
;	  V(X';X)=(2/!PI/HW)*sqrt(1.-((X'-X)/HW)^2), where HW=VROT*X
;	(note that V(X';X) *changes* across the range of X)
;	
;usage summary
;	* call as a function
;	* needs MEAN, SIG, PEAK for complete specification
;
;subroutines
;	MK_GAUSS
;	also makes recursive call to itself if PDER is required
;
;history
;	vinay kashyap (JunMM, based on MK_SLANT)
;	converted array notation to IDL 5 (VK; Apr02)
;	added partial wrt VROT to PDER (VK; Jun02)
;	changed keyword NORM to NORMFLX (VK; Oct02)
;-

np=n_params()
if np lt 1 then begin
  print, 'Usage: rG=mk_rogauss(x,mean,sig,peak,pder,vrot=v,missing=m,/normflx,/fwhm)'
  print, '  generates a gaussian line convolved with a rotational broadening kernel'
  return,[-1L]
endif

;	keywords
if n_elements(vrot) eq 0 then rotv=0. else rotv=float(vrot[0])
if abs(rotv) gt 1 then rotv=rotv/3e5	;must have been in [km/s]
if abs(rotv) gt 1 then rotv=rotv/1e5	;must have been in [cm/s]
if abs(rotv) gt 1 then begin
  message,'what is this, a naked singularity?',/info
  message,'ignoring the alleged rotation ('+strtrim(vrot[0],2)+')',/info
  rotv=0.
endif

;	initialize
nx=n_elements(x) & x0=x[nx/2] & mxx=max(x,min=mnx)

;	make uniform grid
dx=min(abs(x[1:*]-x)) & nxx=long((mxx-mnx)/dx+0.5)+1L
xx=findgen(nxx)*dx+mnx

;	make the gaussian
g=dx*mk_gauss(xx,mean,sig,peak, _extra=e)

;	convolve with rotation
rGf=fltarr(nxx) & hwrot=rotv*xx
if rotv lt 0 then hwrot[*]=abs(rotv*mean)
if rotv ne 0 then begin
  for i=0L,nx-1L do begin
    vv=fltarr(nxx)
    if hwrot[i] ne 0 then begin
      tmp=(1.-((xx-x[i])/hwrot[i])^2) > 0
      vv=dx*(2./!pi/hwrot[i])*sqrt(tmp)
    endif else vv[i]=1.
    rGf=rGf+vv*g[i]
  endfor
endif else rGf=g

;	rebin to original grid
rGf=rGf/dx
rG=interpol(rGf,xx,x)

;	compute partial derivatives
if np ge 5 then begin
  pder = fltarr(nx,4)
  m0=mean & s0=sig & p0=peak
  ;	partial wrt MEAN
  dmean=0.05*s0 & if dmean eq 0 then dmean=dx
  pp=mk_rogauss(x,m0+dmean,s0,p0,vrot=rotv, _extra=e)
  pm=mk_rogauss(x,m0-dmean,s0,p0,vrot=rotv, _extra=e)
  rGm=0.5*(pp-pm)/dmean & pder[*,0]=rGm[*]
  ;	partial wrt SIG
  dsig=0.1*s0 & if dsig eq 0 then dsig=0.1*dx
  pp=mk_rogauss(x,m0,s0+dsig,p0,vrot=rotv, _extra=e)
  pm=mk_rogauss(x,m0,s0-dsig,p0,vrot=rotv, _extra=e)
  rGs=0.5*(pp-pm)/dsig & pder[*,1]=rGs[*]
  ;	partial wrt PEAK
  dpeak=0.01*p0 & if dpeak eq 0 then dpeak=1e-5
  pp=mk_rogauss(x,m0,s0,p0+dpeak,vrot=rotv, _extra=e)
  pm=mk_rogauss(x,m0,s0,p0-dpeak,vrot=rotv, _extra=e)
  rGp=0.5*(pp-pm)/dpeak & pder[*,2]=rGp[*]
  ;	partial wrt VROT
  dvrot=0.01*rotv & if dvrot eq 0 then dvrot=0.001
  pp=mk_rogauss(x,m0,s0,p0,vrot=rotv+dvrot, _extra=e)
  pm=mk_rogauss(x,m0,s0,p0,vrot=rotv-dvrot, _extra=e)
  rGv=0.5*(pp-pm)/dvrot & pder[*,3]=rGv[*]
endif

return,rG
end
