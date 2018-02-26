function hrs2deg,hms,fsep=fsep,dec=dec, _extra=e
;+
;function	hrs2deg
;	converts sexagesimal angles to decimal degrees and returns a
;	numerical array of decimal angles
;
;parameters
;	hms	[INPUT; required] string array of sexagesimal angles
;		* may even be a string of FSEP separated sexagesimals
;
;keywords
;	fsep	[INPUT] field separator in HMS
;		* if FSEP=":", only the first 3 are used
;	dec	[INPUT] if set, does not multiply by 15 (the RA v/s Dec thing)
;		* if the H part is -ve, automatically set
;	_extra	[INPUT] allows passing defined keywords to
;		STR2ARR: SEP, SQUISH
;
;subroutines
;	STR_2_ARR
;	see also DEG2HRS
;
;history
;	vinay kashyap (Nov98)
;	now calls STR_2_ARR instead of STR2ARR (VK; Mar08)
;-

;	usage
nh=n_elements(hms) & deg=0.
if nh eq 0 then begin
  print,"Usage: deg=hrs2deg('h:m:s',fsep=fsep,/dec, sep=sep,/R8,/squish)"
  print,"  convert from sexagesimal to decimal"
  return,[0.]
endif

;	check input
hrs=hms & szh=size(hms) & nszh=n_elements(szh)
if szh(nszh-2) ne 7 then begin
  message,'input must be a string of ":" separated sexagesimal numeral',/info
  return,[0.]
endif
if keyword_set(fsep) then begin
  szf=size(fsep) & nszf=n_elements(fsep)
  if szf(nszf-2) ne 7 then fsep=' '
  if fsep ne ':' then begin
    for i=0L,nh-1L do $
      if i eq 0 then hrs=str_sep(hms,fsep) else hrs=[hrs,str_sep(hms,fsep)]
  endif
endif
nh=n_elements(hrs)

;	define output
deg=fltarr(nh)

;	convert to decimal
for i=0L,nh-1L do begin			;{for each input
  hh=hrs(i) & dd=0. & mm=0. & ss=0.
  cdec=strmid(strtrim(hh,2),0,1)	;check for initial "-"
  isign=1 & if cdec eq '-' then isign=-1

  arr=str_2_arr(hh,sep=':', _extra=e)	;split it
  narr=n_elements(arr)

  if narr ge 1 then dd=abs(float(arr(0)))
  if narr ge 2 then mm=float(arr(1))
  if narr ge 3 then ss=float(arr(2))

  if not keyword_set(dec) and cdec ne '-' then begin	;(goes from 0H to 24H
    while dd gt 24 do dd=dd-24.		;modulo 24
  endif else begin					;)(from -180 to +180
    while dd gt 180 do dd=dd-360.
  endelse				;)

  deg(i)=isign*(dd+mm/60.+ss/3600.)
  if not keyword_set(dec) and cdec ne '-' then deg(i)=deg(i)*15.

endfor					;I=0,NH-1}

return,deg
end
