pro ti_and,tstr1,tstp1,tstr2,tstp2,tstart,tstop, _extra=e
;+
;procedure	ti_and
;	returns the intersection of two sets of time intervals as
;	defined by the start and stop times
;
;syntax
;	ti_and,tstr1,tstp1,tstr2,tstp2,tstart,tstop
;
;parameters
;	tstr1	[INPUT; required] start times for first set
;	tstp1	[INPUT; required] stop times for first set
;	tstr2	[INPUT; required] start times for second set
;	tstp2	[INPUT; required] stop times for second set
;	tstart	[OUTPUT] start times of resultanat intervals
;	tstop	[OUTPUT] stop times of resultanat intervals
;
;keywords
;	_extra	[JUNK] here only to prevent crashing the program
;
;restrictions
;	inputs are assumed to be sorted in ascending order
;
;related
;	TI_CLEAN	- cleans up intervals
;	TI_OR		- merges intervals
;	TI_WRITE	- writes intervals to file
;	TI_COVER()	- computes coverage in given bin
;	TI_FILTER()	- filter an array by GTI
;
;history
;	vinay kashyap (May98)
;	added keyword _EXTRA (VK; Nov98)
;-

;	usage
n1=n_elements(tstr1) & m1=n_elements(tstp1)
n2=n_elements(tstr2) & m2=n_elements(tstp2)
if n1 eq 0 or n2 eq 0 or m1 eq 0 or m2 eq 0 then begin
  print,'Usage: ti_and,tstr1,tstp1,tstr2,tstp2,tstart,tstop'
  print,'  returns the intersection of two sets of time intervals'
  return
endif

;	stupid user tricks
if n1 ne m1 then begin
  message,'TSTR1 and TSTP1 do not match',/info & return
endif
if n2 ne m2 then begin
  message,'TSTR2 and TSTP2 do not match',/info & return
endif
;
if tstr1(0) gt tstp1(0) then begin
  message,'min(TSTR1) > min(TSTP1)',/info & return
endif
if tstr2(0) gt tstp2(0) then begin
  message,'min(TSTR2) > min(TSTP2)',/info & return
endif
if tstr1(n1-1) gt tstp1(m1-1) then begin
  message,'max(TSTR1) > max(TSTP1)',/info & return
endif
if tstr2(n2-1) gt tstp2(m2-1) then begin
  message,'max(TSTR2) > max(TSTP2)',/info & return
endif

tzero=tstr1(0) < (tstr2(0))			;get offset

;	determine all the intervals
t11=tstr1-tzero & t12=tstp1-tzero & t21=tstr2-tzero & t22=tstp2-tzero
tt=[t11(*),t12(*),t21(*),t22(*)]
tt=tt(sort(tt)) & tt=tt(uniq(tt)) & ntt=n_elements(tt)
i1=intarr(ntt-1) & i2=intarr(ntt-1)

;	step through and determine whether to keep the interval or not
for i=0L,n1-1L do begin
  oo=where(tt ge t11(i) and tt lt t12(i),moo)
  if moo eq 0 and t12(i) ne t11(i) then message,'bug!'
  if moo gt 0 then i1(oo)=1
endfor
for i=0L,n2-1L do begin
  oo=where(tt ge t21(i) and tt lt t22(i),moo)
  if moo eq 0 and t22(i) ne t21(i) then message,'bug!'
  if moo gt 0 then i2(oo)=1
endfor

;	find the intersection
ii=intarr(ntt-1) & ii=i1*i2
oo=where(ii gt 0,moo)
if moo eq 0 then begin
  message,'no common intervals found',/info & return
endif
;
tstart=tt(oo) & tstop=tt(oo+1)
;
ti_clean,tstart,tstop

;	restore offset
tstart=tstart+tzero
tstop=tstop+tzero

return
end
