function date2decimal,datestr,fod=fod,doy=doy,jday=jday, _extra=e
;+
;function	date2decimal
;		compute and return decimal representation of date
;		as fraction of the year.
;
;syntax
;	date=date2decimal(datestr,fod=fod,doy=doy,jday=jday)
;
;parameters
;	datestr	[INPUT; required] date string in format 'YYYY-MM-DDTHH:MM:SS'
;		* can be an array
;		* if given in an ununderstandable format, assumed to be today, now.
;
;keywords
;	fod	[OUTPUT] fraction of day
;	doy	[OUTPUT] day of year
;	jday	[OUTPUT] Julian day
;
;history
;	vinay kashyap (2014nov)
;	bug fix: Feb was resetting doy counter (VK; 2015jun)
;-

;	usage
ok='ok' & np=n_params() & nd=n_elements(datestr) & sd=size(datestr,/type)
if nd eq 0 then ok='Insufficient parameters' else $
 if sd ne 7 then ok='DATESTR must be a string'
if ok ne 'ok' then begin
  print,'Usage: date=date2decimal(datestr,fod=fod,doy=doy,jday=jday)'
  print,'  Given date string as YYYY-MM-DDTHH:MM:SS, computes and
  print,'  returns date in decimal years, fraction of day, day of year,
  print,'  and Julian Day'
  if np ne 0 then message,ok,/informational
  return,-1L
endif

;	initialize
mons=['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']
dayspermon=[31,28,31,30,31,30,31,31,30,31,30,31]
dayspermon_leap=dayspermon & dayspermon_leap[1]=dayspermon[1]+1

;	set up default
today=systime()
cc=strsplit(today,' ',/extract) & yr0=long(cc[4]) & day0=long(cc[2])
ok=where(strpos(mons,cc[1]) ge 0,mok) & if mok gt 0 then mon0=ok[0]+1L
hms0=cc[3] & hcc=strsplit(hms0,':',/extract) & hr0=long(hcc[0]) & min0=long(hcc[1]) & sec0=float(hcc[2])

;	outputs
date=dblarr(nd) & fod=dblarr(nd) & doy=intarr(nd) & jday=dblarr(nd)

;	parse the input
for i=0L,nd-1L do begin	;{I=0,ND-1

  cc=strsplit(datestr[i],'T',/extract)

  if n_elements(cc) eq 2 then begin	;(N(CC)=2
    cc1=strsplit(cc[0],'-',/extract) & cc2=strsplit(cc[1],':',/extract)
    if n_elements(cc1) eq 3 then begin	;(N(CC1)=3
      zyr=long(cc1[0]) & zmon=long(cc1[1]) & zday=long(cc1[2])
    endif else begin			;N(CC1)=3)(=/=3
      zyr=yr0 & zmon=mon0 & zday=day0
    endelse				;N(CC1)=/=3)
    if n_elements(cc2) eq 3 then begin	;(N(CC2)=3
      zhr=long(cc2[0]) & zmin=long(cc2[1]) & zsec=float(cc2[2])
    endif else begin			;N(CC2)=3)(=/=3
      zhr=hr0 & zmin=min0 & zsec=sec0
    endelse				;N(CC2)=/=3)
  endif else begin			;N(CC)=2)(N(CC)=/=2
    zyr=yr0 & zmon=mon0 & zday=day0 & zhr=hr0 & zmin=min0 & zsec=sec0
  endelse				;N(CC)=/=2)

  ;	is this a leap year or what?
  yesleap=0 & daymon=dayspermon	;not is the default
  z4=zyr/4. & z100=zyr/100. & z400=zyr/400.
  ;print,z4-long(z4),z100-long(z100),z400-long(z400)
  if z4-long(z4) eq 0 then begin	;(divisible by 4
    yesleap=1 & daymon=dayspermon_leap
    if z100-long(z100) eq 0 and z400-long(z400) ne 0 then begin	;(divisible by 100 but not 400
      yesleap=0 & daymon=dayspermon
    endif							;divisible by 100 but not 400)
  endif					;divisible by 4)

  ;	how many days so far?
  if zmon ge 2 then elapsed_days=total(daymon[0:zmon-2]) else elapsed_days=0
  elapsed_days=elapsed_days+zday
  doy[i]=elapsed_days

  ;	how much of the day so far?
  hours_of_day = (zhr + zmin/60.D + zsec/3600.D)
  elapsed_days=elapsed_days+hours_of_day/24.D
  fod[i]=hours_of_day/24.D

  ;	and putting it together..
  date[i]=zyr+(doy[i]+fod[i])/total(daymon)

  ;	and the Julian Days too, for completeness
  jday[i]=julday(zmon,zday,zyr,zhr,zmin,zsec)

endfor			;I=0,ND-1}

;	

return,date
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;	example

datestr=['2012-11-07T17:56:14','2014-11-07T17:56:14','2014-11-26T15:05:12','now                ']
date=date2decimal(datestr,fod=fod,doy=doy,jday=jday)
print,'input	decimal date	fraction_of_day	day_of_year	Julian Day'
for i=0,n_elements(datestr)-1 do print,datestr[i],'	',date[i],'	',fod[i],'	',doy[i],'	',jday[i]

end
