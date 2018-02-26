function deg2hrs,d,h,m,s,dec=dec, _extra=e
;+
;function	deg2hrs
;	converts decimal degrees to sexagesimal system and returns a
;	string array of properly formatted output
;
;parameters
;	d	[INPUT; required] degrees
;		* may even be a string of comma separated numbers
;	h	[OUTPUT] hours
;	m	[OUTPUT] minutes
;	s	[OUTPUT] seconds
;
;keywords
;	dec	[INPUT] if set, does not divide by 15 (the RA v/s Dec thing)
;	_extra	[INPUT] allows passing defined keywords to
;		STR_2_ARR: SEP, R8, SQUISH
;
;subroutines
;	STR_2_ARR
;	see also HRS2DEG
;
;history
;	vinay kashyap (Nov98)
;	changed call to str2arr to str_2_arr (VK; Sep14)
;-

;	usage
nd=n_elements(d) & h=0 & m=0 & s=0.
if nd eq 0 then begin
  print,'Usage: hms=deg2hrs(deg,h,m,s,/dec, sep=sep,/R8,/squish)'
  print,'  convert from decimal to sexagesimal'
  return,['00:00:00.00']
endif

;	check input
deg=d & szd=size(d) & nszd=n_elements(szd)
if szd(nszd-2) eq 7 then begin
  for i=0L,nd-1L do begin
    dd=str_2_arr(d(i),/r4, _extra=e)
    if i eq 0 then deg=dd else deg=[deg,dd]
  endfor
endif
nd=n_elements(deg)

;	define outputs
hms=strarr(nd) & h=intarr(nd) & m=intarr(nd) & s=fltarr(nd)

;	convert to sexagesimal
for i=0L,nd-1L do begin			;{for each input degree
  dd=deg(i)
  if not keyword_set(dec) then begin	;(goes from 0H to 24H
    while dd lt 0 do dd=dd+360.		;can't have -ves
    while dd gt 360 do dd=dd-360.	;can't exceed one rotation
    dd=dd/15.				;24 hrs in a day
    h(i)=fix(dd) > 0
    if h(i) lt 10 then hh='0'+string(h(i),'(i1)') else hh=string(h(i),'(i2)')
  endif else begin			;)(goes from -180 to +180
    if dd lt 0 then begin
      while dd lt -180 do dd=dd+360.
    endif else begin
      while dd gt 180 do dd=dd-360.
    endelse
    h(i)=fix(abs(dd))
    if dd lt 0 then hh='-' else hh=''
    dd=abs(dd)
    if h(i) lt 10 then hh=hh+'0'+string(h(i),'(i1)') else $
      hh=hh+strtrim(string(h(i),'(i3)'),2)
  endelse				;)
  dd=(dd-h(i))*60.		;60 min in an hour
  m(i)=fix(dd) > 0
  dd=(dd-m(i))*60.		;60 sec in a minute
  s(i)=dd > 0
  if m(i) lt 10 then mm='0'+string(m(i),'(i1)') else mm=string(m(i),'(i2)')
  if s(i) lt 10 then ss='0'+string(s(i),'(f4.2)') else ss=string(s(i),'(f5.2)')
  hms(i)=hh+':'+mm+':'+ss
endfor					;I=0,ND-1}

return,hms
end
