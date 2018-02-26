pro euve_lw,area,wvl,order=order,ardb=ardb, _extra=e
;+
;procedure	EUVE_LW
;	reads in effective area curves for the Extreme Ultra-Violet
;	Explorer's Long-Wavlength spectrometer
;
;usage
;	euve_lw,area,wvl,order=order,ardb=ardb,/help
;
;parameters
;	area	[OUTPUT] effective areas [cm^2]
;	wvl	[OUTPUT] wavelengths at which area is defined [Ang]
;
;keywords
;	order	[I/O] if input is not a vector, then returns only the
;		effective area for the specified order in AREA(WVL).
;		otherwise, on output is an array specifying the
;		spectrographic order at which each element in
;		AREA(WVL) is defined.
;		* if not specified, assumes ORDER=1
;	ardb	[INPUT] analysis reference data base directory,
;		in which to look for instrument calibration data
;		[default: /data/fubar/SCAR/ardb/]
;	_extra	[JUNK] here only to prevent crashing the program
;
;commons
;       euve_lw	{lwea, lw, lword}
;
;restrictions
;	* will crash if ARDB is incorrectly set or if lw_ea.FITS file
;	  is not present in ARDB
;	* requires IDLASTRO library
;	* requires SETSYSVAL
;
;history
;	vinay kashyap (Jan97)
;	changed CALDIR default (VK; Nov98)
;	added keyword HELP (VK; MayMM)
;	changed CALDIR to ARDB (VK; DecMM)
;	now stores in common regardless, and takes ARDB default from !ARDB;
;	  added call to SETSYSVAL (VK; Sep01)
;-

;	usage
np=n_params(0)
if np eq 1 then begin
  print,"Usage: euve_lw,area,wvl,order=order,ardb=ardb"
  print,"  read in effective area curves for EUVE's LW"
  return
endif

;	if no arguments, store in common
common euve_lw,lwea,lw,lword

;	check keywords
allord=0		;meaning extract only one order
if keyword_set(order) then begin
  szo=size(order)
  if szo(0) ne 0 then allord=1	;i.e., extract all orders
endif else order=1
;
c1='/data/fubar/SCAR/ardb/'
ivar=0 & defsysv,'!ARDB',exists=ivar	;if !ARDB exists
if ivar ne 0 then setsysval,'ARDB',c1,/getval
if not keyword_set(ardb) then ardb=c1

;	read calibration data
calfil=ardb+'lw_ea.FITS'		;file containing the calibration data
fxbopen,ulw,calfil,1,h			;open
fxbread,ulw,lw,1			;1st extension -- wavelengths [Ang]
fxbread,ulw,lwea,2			;2nd extension -- effective area [cm^2]
fxbread,ulw,lword,3			;3rd extension -- orders
fxbclose,ulw				;close

;	output
if allord eq 0 then begin
  oo=where(abs(lword) eq abs(order),moo)	;which order?
  if moo gt 0 then begin
    area=lwea(oo) & wvl=lw(oo)		;pick 'em
  endif else begin
    area=lwea & wvl=lw & order=lword	;somewhere, someone made a mistaek
  endelse
endif else order=lword

return
end
