pro rd_pimms_file,pimmfl,effar,nrgar,acol=acol,ecol=ecol,wave=wave,$
	_extra=e
;+
;procedure	rd_pimms_file
;	read in the effective area values from PIMMS file
;
;syntax
;	rd_pimms_file,pimmfl,effar,nrgar,ecol=ecol,acol=acol,/wave
;
;parameters
;	pimmfl	[INPUT; required] scalar string specifying the full
;		pathname to PIMMS effective area file.
;		* output of GET_PIMMS_FILE
;	effar	[OUTPUT] effective areas, usually in [cm^2]
;	nrgar	[OUTPUT] energies at which EFFAR are defined
;		* assumed to be [keV]
;		* if WAVE is set, converted to [Ang]
;
;keywords
;	acol	[INPUT; default=2] column number which contains EFFAR
;	ecol	[INPUT; default=1] column number which contains NRGAR
;	wave	[INPUT] if set, assumes that ECOL contains [keV] and
;		converts to [Ang]
;		* if columns are ALREADY in [Ang], DO NOT SET THIS KEYWORD
;	_extra	[JUNK] ignore -- here only to prevent crashing the program
;
;subroutines
;	STR_2_ARR, WC, STR_SEP, SYZE
;
;history
;	vinay kashyap (Mar99)
;	changed call to STR_2_ARR from STR2ARR (VK; AprMMV)
;-

;	usage
szp=size(pimmfl) & nszp=n_elements(szp)
if szp(nszp-2) ne 7 then begin
  print,'Usage: rd_pimms_file,pimmfl,effar,nrgar,ecol=ecol,acol=acol,/wave'
  print,'  read effective area values from PIMMS file'
  return
endif

;	keywords
cola=2 & if keyword_set(acol) then cola=fix(acol) > 1
cole=1 & if keyword_set(ecol) then cole=fix(ecol) > 1
if not keyword_set(wave) then wave=0

;	error checks
np=n_elements(pimmfl) & if np gt 1 then message,$
	'Canaa handle array of filenames, Capn.  Using only first one!',/info

;	figure out dimensions of input file
nrow=wc(pimmfl(0))		;number of rows
openr,upf,pimmfl(0),/get_lun	;{number of columns
  c='' & readf,upf,c
  cc=str_2_arr(c,/r4,/squish) & ncol=n_elements(cc)
close,upf & free_lun,upf	;number of columns}

;	more error checks
if cola gt ncol or cole gt ncol then begin
  message,'There are only '+strtrim(ncol,2)+' columns in file:'+pimmfl(0),/info
  while cola gt ncol do begin
    print,'' & print,c & print,''
    read,prompt='Type column number for effective areas> ',cola
    if cola le 0 then stop,'Halting.  type .CON to continue'
  endwhile
  while cole gt ncol do begin
    print,'' & print,c & print,''
    read,prompt='Type column number for energies> ',cole
    if cole le 0 then stop,'Halting.  type .CON to continue'
  endwhile
endif

;	read the file
var=fltarr(ncol,nrow)
openr,upf,pimmfl(0),/get_lun & readf,upf,var & close,upf & free_lun,upf

;	outputs
effar=reform(var(cola-1,*)) & nrgar=reform(var(cole-1,*))

;	convert to [Ang] if asked
if keyword_set(wave) then begin
  oo=where(nrgar ne 0,moo)
  if moo gt 0 then nrgar(oo)=12.3985/nrgar(oo)
endif

return
end
