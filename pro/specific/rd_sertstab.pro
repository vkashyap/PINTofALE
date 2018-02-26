pro rd_sertstab,file,elem,wvl,flx,sigma,fwhm,fwhmsig
;+
;procedure	rd_sertstab
;	a hack to read in the ASCII table files of Brosius, J.W.,
;	Davila, J.M., Thomas, R.J., \& Monsignori-Fossi, B.C.\ 1996,
;	ApJS 106, 143 (Measuring Active and Quiet-Sun Coronal Plasma
;	Properties with Extreme-Ultraviolet Spectra from SERTS).
;
;usage
;	rd_sertstab,file,elem,wvl,flx,sigma,fwhm,fwhmsig
;
;parameters
;	file	[INPUT; required] name of file containing data
;	elem	[OUTPUT] atoms and ionic species generating the lines
;	wvl	[OUTPUT] observed wavelength [Ang]
;	flx	[OUTPUT] intensity [ergs/s/cm^2/sr]
;	sigma	[OUTPUT] 1-sigma error on FLX
;	fwhm	[OUTPUT] full-width @ half-max for lines [mAng]
;		* for active and quiet sun spectra, FWHM is 55 mA, and
;		  for limb spectra, FWHM is 49 mA
;		* the FWHM and FWHMSIG in the tables are bogus (i.e.,
;		  are only placeholders)
;	fwhmsig	[OUTPUT] 1-sigma error in FWHM
;
;keywords	NONE
;
;history
;	vinay kashyap (Mar97)
;-

;	usage
if n_elements(file) ne 1 then begin
  print,'Usage: rd_sertstab,file,elem,wvl,flx,sigma,fwhm,fwhmsig'
  print,'  read in data from ASCII tables of SERTS spectra (Brosius'
  print,'  et al. 1996)' & return
endif

;	initialize
elem=[''] & wvl=[0.] & flx=wvl & sigma=wvl & fwhm=wvl & fwhmsig=wvl
nlin=0

openr,utab,strtrim(file,2),/get_lun		;{open file

while not eof(utab) do begin			;(begin input
  line='' & readf,utab,line			;read 1 line
  flds=str_sep(line,'	') & i=0		;break into fields
  if n_elements(flds) eq 6 then begin
    elem=[elem,flds(i)] & i=i+1			;element+ion
    wvl=[wvl,float(flds(i))] & i=i+1		;wavelength
    flx=[flx,float(flds(i))] & i=i+1		;intensity
    sigma=[sigma,float(flds(i))] & i=i+1	;sigma(intensity)
    fwhm=[fwhm,float(flds(i))] & i=i+1		;fwhm
    fwhmsig=[fwhmsig,float(flds(i))]		;sigma(fwhm)
    nlin=nlin+1
  endif
endwhile					;end input)

close,utab & free_lun,utab			;close file}

;	reset
if nlin gt 0 then begin
  elem=elem(1:*) & wvl=wvl(1:*) & flx=flx(1:*) & sigma=sigma(1:*)
  fwhm=fwhm(1:*) & fwhmsig=fwhmsig(1:*)
endif

return
end
