function ti_filter,time,tstart,tstop,count=count, _extra=e
;+
;function	ti_filter
;	returns an array of locations where the input list of times
;	fall within a given time interval defined by the start and
;	stop times
;
;syntax
;	it=ti_filter(time,tstart,tstop)
;
;parameters
;	time	[INPUT; required] arrival times of photons, or bin-values
;	tstart	[INPUT; required] start times of intervals
;	tstop	[INPUT; required] stop times of intervals
;
;keywords
;	count	[OUTPUT] number of elements in output array
;	_extra	[JUNK] here only to prevent crashing the program
;
;related
;	TI_CLEAN	- cleans up intervals
;	TI_OR		- merges intervals
;	TI_AND		- intersection of intervals
;	TI_WRITE	- writes intervals to file
;	TI_COVER()	- computes coverage in given bin
;
;history
;	vinay kashyap (May98)
;	added keyword _EXTRA (VK; Nov98)
;-

;	usage
count=0L
nt=n_elements(time) & nstr=n_elements(tstart) & nstp=n_elements(tstop)
if nt eq 0 or nstr eq 0 or nstp eq 0 then begin
  print,'Usage: it=ti_filter(time,tstart,tstop)'
  print,'  returns array locations where TIME lies between TSTART and TSTOP'
  return,-1L
endif

;	stupid user tricks
if nstr ne nstp then begin
  message,'TSTART does not match TSTOP',/info & return,-1L
endif
;
if tstart(0) gt tstop(0) then begin
  message,'min(TSTART) > min(TSTOP)',/info & return,-1L
endif
if tstart(nstr-1) gt tstop(nstr-1) then begin
  message,'max(TSTART) > max(TSTOP)',/info & return,-1L
endif

;	initialize
ii=intarr(nt)

;	step through each time interval
for i=0L,nstr-1L do begin
  t0=tstart(i) & t1=tstop(i)
  oo=where(time ge t0 and time le t1,moo)
  if moo gt 0 then ii(oo)=1
endfor

;	shorten the array
jj=where(ii gt 0,count)

return,jj
end
