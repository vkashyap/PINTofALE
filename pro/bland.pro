function bland,flstr,dwvl,flx=flx,wflx=wflx,arstr=arstr, _extra=e
;+
;function	bland
;	returns a "blending factor" for each line in the emissivity list,
;	defined as the ratio of the flux in the line to that within the
;	given range.  (i.e., small values mean heavily blended lines)
;
;syntax
;	b=bland(flstr,dwvl,flx=flx,wflx=wflx,arstr=arstr,$
;	dem=dem,abund=abund,effar=effar,wvlar=wvlar,$
;	/temp,/noph,/ikev,/regrid)
;
;parameters
;	flstr	[INPUT; required] structure containing emissivity and
;		wavelength information.  see RD_LINE for details.
;		* ion-balance assumed to be included
;		* if FLX and WFLX are set on input, FLSTR is ignored --
;		  still required, but ignored.
;	dwvl	[INPUT; default=median(delta(FLSTR.WVL))] the spectral
;		"resolution" -- any line nearer than this value at any
;		given wavelength is considered to be blended in to this.
;		* if array of size FLSTR.WVL (or twice that size), set to
;		  different values (with possible asymmetric ranges) at
;		  each wavelength.
;
;keywords
;	arstr	[INPUT] structure of structures containing effective
;		area information (see ARIA or GRATFLX)
;		* if given in the right format, overrides call to LINEFLX
;		  with call to GRATFLX
;	flx	[I/O] contains the calculated flux due to each line
;		* may be given on input and if WFLX is also set, FLSTR and
;		  ARSTR are ignored
;	wflx	[INPUT] if size matches that of FLX, assumed to be the
;		wavelengths of the lines.
;		* if FLX and WFLX are set on input, FLSTR and ARSTR are
;		  ignored.
;	_extra	[INPUT] pass defined keywords to
;		GRATFLX (DEM, ABUND, TEMP, NOPH, IKEV, REGRID)
;		LINEFLX (DEM, ABUND, EFFAR, WVLAR, TEMP, NOPH, IKEV, REGRID)
;
;subroutines
;	ARIA
;	GRATFLX
;	LINEFLX
;
;history
;	vinay kashyap (Sep98)
;	bug fix for when FLSTR is ignored (VK; Oct98)
;	added keyword ARSTR, added call to GRATFLX (VK; Nov98)
;	converted to IDL5 notation (VK; OctMM)
;	bug fix: crashed if FLX and WFLX set, but DWVL not set (VK;JanMMI)
;	bug fix: FLSTR can have more than 8 tags (JJD; MayMMVII)
;-

forward_function gratflx	;because GRATFLX may not exist..

n=n_tags(flstr) & n_f=n_elements(flx) & n_w=n_elements(wflx)
if n_f gt 0 and n_f eq n_w then n=8	;ignore FLSTR
if n lt 8 then begin
  print,'Usage: b=bland(line_structure,resolution,flx=lineflux,wflx=wvls,$'
  print,'       dem=dem,abund=abund,effar=effar,wvlar=wvlar,/temp,/noph,/ikev)'
  print,'  computes a blending factor for each line'
  return,-1L
endif

if n_f eq 0 or n_f ne n_w then begin		;(read FLSTR
  ;	extract info from FLSTR
  wvl=abs(flstr.WVL) & nw=n_elements(wvl)

  ;	stupid user tricks
  if nw eq 1 then return,[1.]

  ;	call LINEFLX or GRATFLX?
  callineflx=1			;by default, call LINEFLX
  if keyword_set(arstr) then begin
    nar=n_tags(arstr)
    if nar gt 0 then begin
      tmp=aria(arstr) & mar=n_tags(tmp)	;just to verify that ARSTR is good
      if mar eq nar then begin
	callineflx=0		;seems OK
	if mar eq 1 and tmp.(0).(3) eq 'help' then callineflx=1
      endif
    endif
  endif

  ;	get the fluxes
  if keyword_set(callineflx) then begin
    flx=lineflx(flstr.LINE_INT,flstr.LOGT,wvl,flstr.Z, _extra=e)
    ww=wvl
  endif else begin
    ostr=gratflx(arstr,flstr.LINE_INT,flstr.LOGT,wvl,flstr.Z,$
	fout=flx,wout=ww, _extra=e)
  endelse
endif else begin				;)(FLX and WFLX are given?
  if n_f eq 0 or n_w eq 0 then begin
    message,'FLX or WFLX not given',/info & return,[1.]
  endif
  ww=abs(wflx)
endelse						;FLX and WFLX are given)
wvl=ww & nw=n_elements(ww)

;	"resolution"
mw=n_elements(dwvl)
if mw eq 0 then begin		;(default is to use median of wvl differences
  w=ww(sort(wvl)) & dw=w[1:*]-w & dw=dw[sort(dw)] & dwvl=median(dw) & mw=1L
endif				;DWVL)
;
delp=fltarr(nw)+dwvl[0] & delm=delp	;constant at all locations
if mw eq 2 then delm[*]=dwvl[1]		;allow asymmetric ranges
if mw eq nw then begin
  delp=dwvl			;upper and lower limits identical, but may
  delm=dwvl			;vary with location
endif
if mw eq 2*nw then begin
  delp=dwvl[0L:nw-1L]		;upper and lower limits may also vary
  delm=dwvl[nw:*]		;with location
endif

;	sift through the wavelength list and get the blending factors
blend=fltarr(nw)
for i=0L,nw-1L do begin
  if i eq 50*fix(i/50) then kilroy; was here.
  dw=ww-ww[i]
  oo=where(dw lt delp[i] and dw gt -delm[i],moo)
  if moo eq 0 then message,'bug!'
  if flx[i] gt 0 then $		;what if emissivities are identically zero?
	blend[i]=flx[i]/total(flx[oo])
endfor

return,blend
end
