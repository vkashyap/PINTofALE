pro euve_ds,area,wvl,ardb=ardb,help=help, _extra=e
;+
;procedure	EUVE_DS
;	reads in effective area curves for the Extreme Ultra-Violet
;	Explorer's Deep Survey instrument
;
;usage
;	euve_ds,area,wvl,ardb=ardb,/help
;
;parameters
;	area	[OUTPUT] effective areas [cm^2]
;	wvl	[OUTPUT] wavelengths at which area is defined [Ang]
;
;keywords
;	ardb	[INPUT] analysis reference data base directory,
;		in which to look for instrument calibration data
;		[default: /data/fubar/SCAR/ardb/]
;	help	[INPUT] if set, prints usage and exits
;	_extra	[JUNK] here only to prevent crashing the program
;
;commons
;	euve_ds	{dsea, ds}
;
;restrictions
;	* will crash if ARDB is incorrectly set or if ds_ea_quad1.FITS file
;	  is not present in ARDB
;	* requires IDLASTRO library
;	* requires SETSYSVAL
;
;history
;	vinay kashyap (Jul97)
;	changed CALDIR default (VK; Nov98)
;	added keyword HELP (VK; MayMM)
;	changed CALDIR to ARDB (VK; DecMM)
;	now stores in common regardless, and takes ARDB default from !ARDB
;	  added call to SETSYSVAL (VK; Sep01)
;-

;	usage
np=n_params(0)
if np eq 1 or keyword_set(help) then begin
  print,"Usage: euve_ds,area,wvl,ardb=ardb,/help"
  print,"  read in effective area curves for EUVE's DS"
  return
endif

;	if no arguments, store in common
common euve_ds,dsea,dw

;	check keywords
c1='/data/fubar/SCAR/ardb/'
ivar=0 & defsysv,'!ARDB',exists=ivar	;if !ARDB exists
if ivar ne 0 then setsysval,'ARDB',c1,/getval
if not keyword_set(ardb) then ardb=c1

;	read calibration data
calfil=ardb+'ds_ea_quad1.FITS'	;file containing the calibration data
ds=mrdfits(calfil,1,h)			;read into structure
dsea=ds.EFFAREA				;effective area [cm^2]
dw=ds.WAVELENGTH			;wavelength [Ang]

;	output
area=dsea & wvl=dw

return
end
