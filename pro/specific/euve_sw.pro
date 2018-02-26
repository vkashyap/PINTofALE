pro euve_sw,area,wvl,order=order,ardb=ardb,help=help, _extra=e
;+
;procedure	EUVE_SW
;	reads in effective area curves for the Extreme Ultra-Violet
;	Explorer's Short-Wavlength spectrometer
;
;usage
;	euve_sw,area,wvl,order=order,ardb=ardb,/help
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
;	help	[INPUT] if set, prints usage and exits
;	_extra	[JUNK] here only to prevent crashing the program
;
;commons
;	euve_sw	{swea, sw, sword}
;
;restrictions
;	* will crash if ARDB is incorrectly set or if sw_ea.FITS file
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
if np eq 1 or keyword_set(help) then begin
  print,"Usage: euve_sw,area,wvl,order=order,ardb=ardb,/help"
  print,"  read in effective area curves for EUVE's SW"
  return
endif

;	if no arguments, store in common
common euve_sw,swea,sw,sword

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
calfil=ardb+'sw_ea.FITS'		;file containing the calibration data
fxbopen,usw,calfil,1,h			;open
fxbread,usw,sw,1			;1st extension -- wavelengths [Ang]
fxbread,usw,swea,2			;2nd extension -- effective area [cm^2]
fxbread,usw,sword,3			;3rd extension -- orders
fxbclose,usw				;close

;	output
if allord eq 0 then begin
  oo=where(abs(sword) eq abs(order),moo)	;which order?
  if moo gt 0 then begin
    area=swea(oo) & wvl=sw(oo)		;pick 'em
  endif else begin
    area=swea & wvl=sw & order=sword	;somewhere, someone made a miatake
  endelse
endif else order=sword

return
end
