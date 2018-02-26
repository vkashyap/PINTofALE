function date2timestamp,date,year=year, _extra=e
;+
;function	date2timestamp
;	given a decimal date (in units of [yr]), returns the calendar
;	date as a timestamp in the format YYYY-MM-DDTHH:MM:SS.S
;
;syntax
;	tstamp=date2timestamp(date,year=year)
;
;parameters
;	date	[INPUT; required] decimal date in years or fraction
;		of year (see keyword YEAR)
;
;keywords
;	year	[INPUT] if given, adds this to DATE
;		* useful for when only fraction of year is given
;	_extra	[JUNK] here only to prevent crashing the program
;
;example
;	.run date2timestamp
;	(calls date2decimal)
;
;history
;	Vinay Kashyap (2015.4863)
;-

;	usage
ok='ok' & np=n_params() & nd=n_elements(date) & sd=size(date,/type)
if np eq 0 then ok='Insufficient parameters' else $
 if nd eq 0 then ok='DATE is undefined' else $
  if sd ne 4 and sd ne 5 then ok='DATE must be floating point or double precision'
if ok ne 'ok' then begin
  print,'Usage: tstamp=date2timestamp(date,year=year)'
  print,'  given decimal year, returns YYYY-MM-DDTHH:MM:SS.S'
  if np ne 0 then message,ok,/informational
  return,-1L
endif

;	keywords
yy=date
yr0=dblarr(nd)
if keyword_set(year) then begin
  if size(year,/type) lt 6 then begin
    ny=n_elements(year) & yr0[*]=year[0]
    if ny le nd then yr0[0L:ny-1L]=year[*] else yr0[*]=year[0L:nd-1L]
  endif else message,'YEAR should be a number; ignoring',/informational
endif

;	compute Julian Days for beginning and end of the specified year
y0=long(yy)
j0=julday(1,1,y0,0,0,0) & j1=julday(12,31,y0,12,59,59) & daysperyr=ceil(j1-j0)

;	convert given date to JD and thence to calendar date
jd=j0+(yy-y0)*double(daysperyr)
jd=jd-1	;subtract 1/2 day and 1/2 again to account for conversion to and from JD
caldat,jd,mon,day,year,hour,min,sec

tstamp=	string(year,'(i4.4)')+'-'+$
	string(mon,'(i2.2)')+'-'+$
	string(day,'(i2.2)')+'T'+$
	string(hour,'(i2.2)')+':'+$
	string(min,'(i2.2)')+':'+$
	string(sec,'(f04.1)')

return,tstamp
end

;	example case

input_dates=['1776-07-04T12:02:39','1941-12-07T07:48:00','1947-08-15T00:00:00.1','2015-06-26T12:54:03',systime()]
date=date2decimal(input_dates)
tstamp=date2timestamp(date)
for i=0,n_elements(input_dates)-1 do print,input_dates[i]+string(string(byte(9b))+'->','(a10)')+string(date[i],'(f12.6)')+'	-> '+tstamp[i]

end
