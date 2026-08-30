pro wcs_eqpos,skyx,skyy,header,ra,dec,CDELT=CDELT,CRPIX=CRPIX,CRVAL=CRVAL,verbose=verbose, _extra=e
;+
;procedure	wcs_eqpos
;	convert SKY coordinates to RA,DEC using the CDELT, CRPIX, CRVAL or equivalent
;	header keywords associated with RA---TAN and DEC--TAN projections
;
;syntax
;	wcs_eqpos,skyx,skyy,header,ra,dec,CDELT=CDELT,CRPIX=CRPIX,CRVAL=CRVAL,verbose=verbose
;
;parameters
;	skyx	[INPUT] x-coordinates
;	skyy	[INPUT] y-coordinates
;	header	[INPUT] FITS-like header that contains information about coordinate transformations
;	ra	[OUTPUT] RA for each SKYX [deg]
;	dec	[OUTPUT] Dec for each SKYY [deg]
;
;keywords
;	CDELT	[I/O] if supplied on input, used to do the transformation, IFF HEADER does not override it
;		* overwritten on output 
;	CRPIX	[I/O] if supplied on input, used to do the transformation, IFF HEADER does not override it
;		* overwritten on output 
;	CRVAL	[I/O] if supplied on input, used to do the transformation, IFF HEADER does not override it
;		* overwritten on output 
;	verbose	[INPUT] controls chatter
;
;	_extra	[JUNK] here only to prevent crashing the program
;
;history
;	Vinay Kashyap (2026-aug-30)
;-

ok='ok' & np=n_params() & nx=n_elements(skyx) & ny=n_elements(skyy) & nh=n_elements(header)
if nx eq 0 then ok='SKYX is not available' else $
 if ny eq 0 then ok='SKYY is not available' else $
  if nx ne ny then ok='SKYX and SKYY are not compatible'
if ok ne 'ok' then begin
  print,'Usage: wcs_eqpos,skyx,skyy,header,ra,dec,CDELT=CDELT,CRPIX=CRPIX,CRVAL=CRVAL,verbose=verbose'
  print,'  convert from skyX,Y to RA,DEC'
  if np ne 0 then message,ok,/informational
  return
endif
xx=skyx[*] & yy=skyy[*]

;	keywords
vv=0L & if keyword_set(verbose) then vv=long(verbose[0])>1L
;
;	base defaults are HRC-S
tcdltr=-3.6611111111111000d-05 & tcdltd=3.6611111111111000d-05
tcrpxr=32768.5D & tcrpxd=32768.5D
tcrvlr=0.0D & tcrvld=0.0D
if keyword_set(cdelt) then begin
  tcdltr=cdelt[0] & tcdltd=tcdltr & if n_elements(cdelt) eq 2 then tcdltd=cdelt[1]
endif
if keyword_set(crpix) then begin
  tcrpxr=crpix[0] & tcrpxd=tcrpxr & if n_elements(crpix) eq 2 then tcrpxd=crpix[1]
endif
if keyword_set(crval) then begin
  tcrvlr=crval[0] & tcrvld=tcrvlr & if n_elements(crval) eq 2 then tcrvld=crval[1]
endif

if nh gt 0 then begin
  hdr=header[*]

  fxbfind,hdr,'TCTYP',cols,vals,nf
  if nf gt 0 then begin
    o1=where(strpos(vals,'RA---TAN') ge 0,mo1)
    o2=where(strpos(vals,'DEC--TAN') ge 0,mo2)
    if mo1 ne 0 and mo2 ne 0 then begin
      fxbfind,hdr,'TCDLT',cold,vald,nfd & tcdltr=vald[o1[0]] & tcdltd=vald[o2[0]]
      fxbfind,hdr,'TCRPX',colp,valp,nfp & tcrpxr=valp[o1[0]] & tcrpxd=valp[o2[0]]
      fxbfind,hdr,'TCRVL',colv,valv,nfv & tcrvlr=valv[o1[0]] & tcrvld=valv[o2[0]]
      cdelt=[tcdltr,tcdltd]
      crpix=[tcrpxr,tcrpxd]
      crval=[tcrvlr,tcrvld]
    endif
  endif
endif

if vv gt 10 then begin
  print,'RA  = '+strtrim(tcrvlr,2)+' + tan('+strtrim(tcdltr,2)+') * (skyX - '+strtrim(tcrpxr,2)+')'
  print,'Dec = '+strtrim(tcrvld,2)+' + tan('+strtrim(tcdltd,2)+') * (skyY - '+strtrim(tcrpxd,2)+')'
endif

ra  = tcrvlr + tan(tcdltr*!dpi/180.D)*(xx-tcrpxr)
dec = tcrvld + tan(tcdltd*!dpi/180.D)*(yy-tcrpxd)

return
end
