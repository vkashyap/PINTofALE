function sdss2flux,mag,uband=uband,gband=gband,rband=rband,iband=iband,zband=zband,$
	verbose=verbose, _extra=e
;+
;function	sdss2flux
;	convert SDSS asinh magnitudes to fluxes in ergs/s/cm^2/Hz
;
;syntax
;	fx=sdss2flux(mag,/gband,/rband,/iband,/zband,/uband,verbose=verbose)
;
;parameters
;	mag	[INPUT; required] SDSS (asinh) magnitude
;		* could be scalar or array
;
;keywords
;	uband	[INPUT] input is u mags
;	gband	[INPUT] input is g mags
;	rband	[INPUT] input is r mags
;	iband	[INPUT] input is i mags
;	zband	[INPUT] input is z mags
;		* default is GBAND
;		* if multiple keywords are set,
;		  GBAND > RBAND > IBAND > ZBAND > UBAND
;	verbose	[INPUT] controls chatter
;	_extra	[JUNK] here only to prevent crashing the program
;
;description
;	see SDSS photometric calibration documents
;	http://www.sdss3.org/dr10/algorithms/fluxcal.php
;	http://www.sdss3.org/dr10/algorithms/magnitudes.php#asinh
;	http://classic.sdss.org/dr6/algorithms/fluxcal.html
;
;example
;	print, 'solar luminosity in g band = ',sdss2flux(4.7,/gband)*4.*!dpi*(3.1d18*10.)^2
;
;history
;	vinay kashyap (2014sep16)
;-

;	usage
ok='ok' & np=n_params() & nm=n_elements(mag)
if np eq 0 then ok='Insufficient parameters' else $
 if nm eq 0 then ok='MAG is undefined'
if ok ne 'ok' then begin
  print,'Usage: fx=sdss2flux(mag,/gband,/rband,/iband,/zband,/uband,verbose=verbose)'
  print,'  convert SDSS (asinh) magnitudes to flux'
  if np ne 0 then message,ok,/informational
  return,-1L
endif

;	which band?
vv=0L & if keyword_set(verbose) then vv=long(verbose[0])>1L
iu=0 & if keyword_set(uband) then iu=1
iz=0 & if keyword_set(zband) then iz=1 & if iz eq 1 then begin & iu=0 & endif
ii=0 & if keyword_set(iband) then ii=1 & if ii eq 1 then begin & iu=0 & iz=0 & endif
ir=0 & if keyword_set(rband) then ir=1 & if ir eq 1 then begin & iu=0 & iz=0 & ii=0 & endif
ig=0 & if keyword_set(gband) then ig=1 & if ig eq 1 then begin & iu=0 & iz=0 & ii=0 & ir=0 & endif
if iu+ig+ir+ii+iz eq 0 then ig=1

;	initialize
inicon,fundae=ff
ubsoft = 1.4e-10 & udmag = -0.04 & ufreq0=ff.c/3524e-8 & udfreq=ufreq0*(599./3557.)
gbsoft = 0.9e-10 & gdmag =  0.00 & gfreq0=ff.c/4714e-8 & gdfreq=gfreq0*(1379./4825.)
rbsoft = 1.2e-10 & rdmag =  0.00 & rfreq0=ff.c/6182e-8 & rdfreq=rfreq0*(1382./6261.)
ibsoft = 1.8e-10 & idmag =  0.00 & ifreq0=ff.c/7592e-8 & idfreq=ifreq0*(1535./7672.)
zbsoft = 7.4e-10 & zdmag =  0.02 & zfreq0=ff.c/9003e-8 & zdfreq=zfreq0*(1370./9097.)
if iu eq 1 then begin & bsoft=ubsoft & dmag=udmag & dfreq=udfreq & endif
if ig eq 1 then begin & bsoft=gbsoft & dmag=gdmag & dfreq=gdfreq & endif
if ir eq 1 then begin & bsoft=rbsoft & dmag=rdmag & dfreq=rdfreq & endif
if ii eq 1 then begin & bsoft=ibsoft & dmag=idmag & dfreq=idfreq & endif
if iz eq 1 then begin & bsoft=zbsoft & dmag=zdmag & dfreq=zdfreq & endif

;	convert
fbyf0 = sinh(-(mag+dmag)*alog(10.)/2.5 - alog(bsoft))*2*bsoft	;f/f0, where f0 corresponds to 0 mag

fband = 3731d-23 * fbyf0 	;f0 is always 3731 Jy, 1 Jy = 1e-26 W Hz-1 m-2 = 1e-23 erg s-1 Hz-1 cm-2
fband = fband * dfreq	;[ergs/s/cm^2/Hz]*[Hz]

if vv gt 1000 then stop,'HALTing; type .CON to continue'

return,fband
end
