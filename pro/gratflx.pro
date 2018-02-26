function gratflx,areas,line,logT,wvl,Z,wout=wout,fout=fout, _extra=e
;+
;function	gratflx
;	computes observed flux from spectral lines using LINEFLX
;	and returns fluxes and wavelengths in a structure of structures.
;	this is essentially a wrapper routine for LINEFLX to handle
;	the case of multiple effective areas (as in different grating
;	orders)
;
;syntax
;	ostr=gratflx(areas,line,logT,wvl,Z,wout=wout,fout=fout,$
;	DEM=DEM,abund=abund,/ikev,/noph,/temp)
;
;parameters
;	areas	[INPUT; required] structure of structures containing
;		effective area information in the form:
;		  effective areas:	AREAS.(#).(0)	[cm^2]
;		  wavelengths:		AREAS.(#).(1)	[Ang]
;		  order:		AREAS.(#).(2)	[integer; default=1]
;		  comment:		AREAS.(#).(3)	[string; default='']
;
;	(the following are passed without comment to LINEFLX)
;
;	line	[INPUT; required] array of line cooling emissivities
;		in units of 1e-23 erg cm^3/s
;		* if 2D array, LINE==LINE(logT,WVL)
;		* WARNING: will return garbage if given 1D array LINE(WVL)
;		  use a for-loop to handle such a case
;		* WARNING: will get converted to 1-element vector if input
;		  as scalar
;	logT	[INPUT; required] array of log10(Temperature [K]) at
;		which emissivities are given.
;		* array size MUST match that of LINE
;	wvl	[INPUT; required] wavelength(s) [Angstrom] of lines at
;		which LINE is given
;		* array size MUST match LINE
;	Z	[INPUT] atomic number of element that generates each WVL
;		* if Z has less elements than WVL, Z is ignored
;		* if Z is to be ignored, all lines are assumed to be from
;		  same element, and abundance values are ignored
;
;keywords
;	wout	[OUTPUT] all the "translated" wavelengths
;	fout	[OUTPUT] all the fluxes
;	_extra	[INPUT] use this to pass defined keywords to subroutines
;		LINEFLX -- DEM, ABUND, IKEV, NOPH, TEMP
;
;history
;	vinay kashyap (Oct98)
;	added keywords WOUT, FOUT; ability to handle "bad" AREAS (VK;Nov98)
;-

;	usage
ok='ok'
nar=n_tags(areas)
nem=n_elements(line) & nT=n_elements(logT) & nw=n_elements(wvl)
if nar eq 0 then ok='Effective areas not in structure' else $
 if nem eq 0 then ok='emissivities not given' else $
  if nT eq 0 then ok='temperatures not given' else $
   if nw eq 0 then ok='wavelengths not given'
if nar gt 0 then begin
  tmp=aria(areas) & mar=n_tags(tmp)
  ok2='AREA structure not in right format -- see ARIA()'
  if mar ne nar then ok=ok2
  if mar eq 1 and tmp.(0).(3) eq 'help' then ok=ok2
endif
if ok ne 'ok' then begin
  print,'Usage: ostr=gratflx(areas,emissivities,logT,wvl,Z,wout=wout,fout=fout,$'
  print,'       DEM=DEM,abund=abund,/ikev,/noph,/temp)'
  print,'  returns line fluxes for each set of supplied areas'
  message,ok,/info
  return,-1L
endif

;	call LINEFLX
for i=0,nar-1 do begin			;{for each given effective area

  ;	unpack AREAS
  ea=areas.(i) & mar=n_tags(ea)
  if mar ge 1 then areff=ea.(0) else areff=ea
  if mar ge 2 then arwvl=ea.(1) else arwvl=findgen(n_elements(areff))+1.
  if mar ge 3 then order=(fix(ea.(2)))(0) else order=1
  if mar ge 4 then comm=strtrim(ea.(3)) else comm=''
  ;if mar gt 4 consider it an option to increase info carried in AREAS

  if order le 0 or strlowcase(comm) eq 'help' then begin
    ;	skip if useless
    f=fltarr(nw) & w=wvl
  endif else begin

    ;	call LINEFLX
    f=lineflx(line,logT,wvl,Z, effar=ea,wvlar=arwvl, _extra=e)

    ;	reset wavelengths
    w=wvl*order		;n\lambda=d*sin(\theta), here \theta==w, n==order

  endelse

  ;	recast in structure
  mf=n_elements(f) & mw=n_elements(w)
  if mf ne mw then f=0*w		;something was missing
  if i eq 0 then ostr=create_struct('GRATFLX'+strtrim(i+1,2),$
  	create_struct('FLUX',f,'WVL',w)) else $
	ostr=create_struct(ostr,'GRATFLX'+strtrim(i+1,2),$
	create_struct('FLUX',f,'WVL',w))

  ;	also in plainer arrays
  if i eq 0 then begin
    wout=w & fout=f
  endif else begin
    wout=[wout,w] & fout=[fout,f]
  endelse

endfor					;I=0,NAR-1}

return,ostr
end
