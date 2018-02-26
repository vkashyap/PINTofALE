function mk_powlam,x,gamma1,gamma2,normE,pder,Xo=Xo,pbreak=pbreak,$
	NH=NH,wvlar=wvlar,effar=effar,exptime=exptime,dellam=dellam,$
	_extra=e
;+
;function	mk_powlam
;	returns a power-law in WAVELENGTH G(X)
;
;syntax
;	g=mk_powlam(x,gamma1,gamma2,normE,pder,Xo=Xo,pbreak=pbreak,$
;	NH=NH,wvlar=wvlar,effar=effar,exptime=exptime,dellam=dellam)
;
;parameters
;	X	[INPUT array; required] where G(X) must be computed
;	gamma1	[INPUT; required] power-law index 1
;	gamma2	[INPUT; required] broken power-law index 2
;	normE	[INPUT; required] power-law normalization 
;		* as appropriate for an _energy [keV]_ scale
;	pder	[OUTPUT; optional] partial derivatives of model wrt parameters
;		at each X; calculated only if 5 parameters are supplied in call.
;		* array of size [N(X),4], with columns containing the partial
;		  derivatives wrt GAMMA1, GAMMA2, NORME, and PBREAK respectively
;
;keywords
;	Xo	[INPUT] Wavelength for which normalization 
;		is defined 1.0; default is 1.0 Angstroms 
;	pbreak	[INPUT] Wavelength at which power-law "breaks"
;	NH	[INPUT] Hydrogen column density for computing ISM absorption
;	wvlar	[INPUT array] wavelength grid of 
;		Chandra/HRC-S/LETG effective area
;	effar	[INPUT array] Chandra/HRC-S/LETG effective area
;	exptime	[INPUT] exposure time
;	dellam	[INPUT] wavelength binsize for computing COUNTS spectrum
;	_extra	[JUNK] 
;
;description
;
;  for X <= pbreak,
;  G(X)=exp(-tau)*norm*(X/Xo)^(gamma1-2)*effar*exptime*dellam
;
;  for X > pbreak,
;  G(X)=exp(-tau)*norm*(pbreak/Xo)^(gamma1-gamma2)*(X/Xo)^(gamma2-2)
;	*effar*exptime*dellam
;
;usage summary
;	* call as a function
;	* generates model only at specified points X
;	* needs for complete specification
;
;subroutines
;	ISMTAU
;
;history
;	Deron Pease (June 2000)
;	call ISMTAU only if NH is defined (VK; MMJul)
;	bug correction: if NH not defined, etau must be set to 1 (VK; Jun01)
;	converted array notation to IDL 5 (VK; Apr02)
;	added partial derivatives wrt PBREAK (VK; Jun02)
;	pushed partial derivatives wrt PBREAK inside if..then (VK; Nov02)
;-

np=n_params()
if np lt 1 then begin
  print, 'Usage: g=mk_powlam(x,gamma1,gamma2,norm,pder,Xo=Xo,pbreak=pbreak,NH=NH,wvlar=wvlar,effar=effar,exptime=exptime,dellam=dellam)
  print, '  generates a broken power-law G(X)' & return,[-1L]
endif

; set up stuff and locate the break
if not keyword_set(dellam) then dellam=1.
if not keyword_set(exptime) then exptime=1.
if not keyword_set(Xo) then Xo=1.0
if not keyword_set(pbreak) then pbreak=15.5

norm=normE[0]*(12.3985)^(1.-gamma1[0])		;;convert norm from nrg -> wav

nx=n_elements(X)

;ISM absorption
if keyword_set(NH) then etau=ismtau(X,NH=NH,He1=.09*NH,HeII=.01*NH) else $
  etau=0.*X+1.

;are the EFFECTIVE AREAS and WAVELENGTHS given?
nar=n_elements(effar) & nwv=n_elements(wvlar)
;;if neither are given,  if one is given,  if both but not same grid given
;;stub 2-element area, first and last X, interpolate
if nar eq 0 and nwv eq 0 then areff=0*X+1 else begin
  if nar ne 0 then areatmp=[effar] else areatmp=[1.,1.]
  if nwv ne 0 then arwvl=[wvlar] else arwvl=[min(abs(X)),max(abs(X))]
  if nar ne nwv then begin
    message,'Effective area v/s wavelength mismatch; doing nothing',/info
    areff=0*X+1.
  endif else begin
    ow=sort(arwvl) & arwvl=arwvl[ow] & areatmp=areatmp[ow]  ;sort
    areff=interpol(areatmp,arwvl,abs(X))>0       ;effective areas at WVLs
  endelse
endelse

g=fltarr(nx)
tmp=min(abs(X-pbreak),index)		;; point of break - Array Element

;compute power-law
C=exp(-etau)*areff*exptime*dellam
if index ge 0 and index le nx-1L then $
  g[0L:index]=C[0L:index]*norm*(X[0L:index]/Xo)^(gamma1[0]-2)
if index ge 0 and index lt nx-1L then $
  g[index+1L:*]=C[index+1L:*]*norm*(pbreak/Xo)^(gamma1[0]-gamma2[0]) $
		*(X[index+1L:*]/Xo)^(gamma2[0]-2)

;compute partial dervatives
;	PDER Order:  wrt gamma1, wrt gamma2, wrt norm
;	and wrt PBREAK
if np ge 5 then begin
  pder=fltarr(nx,4)
  if index ge 0 and index le nx then $
    pder[0L:index,0L]=C[0L:index]*norm*(X[0L:index]/Xo)^(gamma1[0]-2) $
		*alog(X[0L:index]/Xo) & $
  pder[0L:index,1L]=0. & $
  pder[0L:index,2L]=C[0L:index]*(X[0L:index]/Xo)^(gamma1[0]-2)
  if index ge 0 and index lt nx then $
    pder[index+1L:*,0L]=C[index+1L:*]*norm*(pbreak/Xo)^(gamma1[0]-gamma2[0]) $
		*(X[index+1L:*]/Xo)^(gamma2[0]-2)*alog(pbreak/Xo) & $
  pder[index+1L:*,1L]=C[index+1L:*]*norm*(pbreak/Xo)^(gamma1[0]-gamma2[0]) $
		*(X[index+1L:*]/Xo)^(gamma2[0]-2) $
		*( alog(X[index+1L:*]/Xo) - alog(pbreak/Xo) ) & $
  pder[index+1L:*,2L]=C[index+1L:*]*(pbreak/Xo)^(gamma1[0]-gamma2[0]) $
		*(X[index+1L:*]/Xo)^(gamma2[0]-2)
  ;	partial wrt PBREAK
  dbrk=0.01*pbreak & if dbrk eq 0 then dbrk=0.001
  pp=mk_powlam(x,gamma1,gamma2,normE,Xo=Xo,pbreak=pbreak+dbrk,$
	NH=NH,wvlar=wvlar,effar=effar,exptime=exptime,dellam=dellam)
  pm=mk_powlam(x,gamma1,gamma2,normE,Xo=Xo,pbreak=pbreak-dbrk,$
	NH=NH,wvlar=wvlar,effar=effar,exptime=exptime,dellam=dellam)
  pder[*,3]=0.5*(pp-pm)/dbrk
endif

return,g
end
