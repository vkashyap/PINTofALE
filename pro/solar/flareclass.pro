function flareclass,class,cgs=cgs,atSun=atSun,verbose=verbose, _extra=e
;+
;function	flareclass
;	compute the peak flux of a solar flare based on its GOES classification
;
;	It is assumed that the classification is given as a letter followed by a
;	number.  The letter is used to define the intensity level, and the number
;	is a multiple of it.  The base levels are:
;	A : 1e-8 - 1e-7 W/m^2
;	B : 1e-7 - 1e-6 W/m^2
;	C : 1e-6 - 1e-5 W/m^2
;	M : 1e-5 - 1e-4 W/m^2
;	X : >1e-4 W/m^2
;	see http://www.spaceweather.com/glossary/flareclasses.html
;
;syntax
;	fx=flareclass(class,cgs=0,/atSun,verbose=verbose)
;
;parameters
;	class	[INPUT; required] string variable that contains the GOES class
;
;keywords
;	cgs	[INPUT] by default, the output is returned in cgs, as ergs/s/cm^2 at Earth
;		* set this explicitly to 0 to return the output in SI units
;	atSun	[INPUT] if set, multiplies by 4*pi*AU^2 to return peak luminosity of flare
;		at the Sun
;	verbose	[INPUT] controls chatter
;	_extra	[JUNK] here only to prevent crashing the program
;
;history
;	Vinay Kashyap (MMXV.VI)
;-

;	usage
ok='ok' & np=n_params() & nc=n_elements(class) & szc=size(class,/type)
if np eq 0 then ok='Insufficient inputs' else $
 if nc eq 0 then ok='CLASS is undefined' else $
  if szc ne 7 then ok='CLASS must be a string variable'
if ok ne 'ok' then begin
  print,'Usage: fx=flareclass(class,cgs=0,/atSun,verbose=verbose)'
  print,'  compute peak flux of solar flare based on GOES classification'
  if np ne 0 then message,ok,/informational
  return,-1L
endif

;	keywords
vv=0L & if keyword_set(verbose) then vv=long(verbose[0])>1L
;
icgs=1	;default is to convert to cgs [ergs/s/cm^2] or [ergs/s]
	;otherwise keep in SI [W/m^2] or [J/s]
if n_elements(cgs) gt 0 then begin
  szc=size(cgs,/type)
  if szc eq 7 then begin
    ;	undocumented feature: can say cgs='n' or cgs='NO' or cgs='NOPE!' etc.
    if strpos(strlowcase(cgs[0]),'n') ge 0 then icgs=0
  endif else begin
    if long(cgs[0]) eq 0 then icgs=0
  endelse
endif
;
distcor=1.D
if keyword_set(atSun) then begin
  inicon,fundae=ff & dist=ff.AU
  if icgs eq 0 then dist=dist/1e2	;convert to [m] if output should be in SI units
  distcor=4.D*!dpi*dist^2
endif

;	define the output
szc=size(class) & fx=0.D
if szc[0] gt 0 then fx=make_array(szc[1:szc[0]],/double,value=0.D)

;	extract the first letter
cc=strupcase(strmid(class,0,1))
ci=strtrim(strmid(class,1,10),2)	;usually only need 3, maybe 4 for
					;very large X class flares, and if
					;we ever reach 10, no one will be
					;alive on Earth to call this a bug
fi=double(ci)

;	this is where a dictionary would be handy, yeah?
oa=where(cc eq 'A',moa) & if moa gt 0 then fx[oa]=1d-8*fi[oa]
ob=where(cc eq 'B',mob) & if mob gt 0 then fx[ob]=1d-7*fi[ob]
oc=where(cc eq 'C',moc) & if moc gt 0 then fx[oc]=1d-6*fi[oc]
om=where(cc eq 'M',mom) & if mom gt 0 then fx[om]=1d-5*fi[om]
ox=where(cc eq 'X',mox) & if mox gt 0 then fx[ox]=1d-4*fi[ox]

if icgs eq 1 then fx=1e3*fx	;[1 W/m^2] = [1e3 ergs/s/cm^2]
if keyword_set(atSun) then fx=fx*distcor

if vv gt 1000 then stop,'halting in FLARECLASS; type .CON to continue'

return,fx
end
