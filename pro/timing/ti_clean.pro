pro ti_clean,tstart,tstop,leap=leap, _extra=e
;+
;procedure	ti_clean
;	given an arbitrary set of starting and stopping times of intervals,
;	first organizes them in ascending order and then merges the various
;	intervals that are contiguous.
;
;syntax
;	ti_clean,tstart,tstop,leap=leap
;
;parameters
;	tstart	[I/O; required] starting times of intervals
;	tstop	[I/O; required] ending times of intervals
;
;keywords
;	leap	[INPUT; default=1.] if the difference between a stopping
;		time and the *next* starting time of sequential intervals
;		is < LEAP, merge the two intervals into one.
;	_extra	[JUNK] here only to prevent crashing the program
;
;related
;	TI_OR		- merges intervals
;	TI_AND		- intersection of intervals
;	TI_WRITE	- writes intervals to file
;	TI_COVER()	- computes coverage in given bin
;	TI_FILTER()	- filter an array by GTI
;
;history
;	vinay kashyap (May98)
;	added keyword _EXTRA (VK; Nov98)
;	now also handles overlapping bad times (VK; Mar99)
;-

;	usage
nstr=n_elements(tstart) & nstp=n_elements(tstop)
if nstr eq 0 or nstp eq 0 then begin
  print,'ti_clean,tstart,tstop,leap=leap'
  print,'  clean up time-interval lists' & return
endif

;	keywords
if not keyword_set(leap) then jump=1. else jump=leap(0)+0.

;	stupid user tricks
if nstr ne nstp then begin
  message,'TSTART and TSTOP do not match',/info & return
endif
;
if nstr eq 1 then return	;oy vey.  don't waste our time with this.
;
tsrmn=min(tstart,max=tsrmx) & tspmn=min(tstop,max=tspmx)
if tsrmn gt tspmn then begin
  message,'min(TSTART) > min(TSTOP)??',/info & return
endif
if tsrmx gt tspmx then begin
  message,'max(TSTART) > max(TSTOP)??',/info & return
endif

;	offset
tzero=tsrmn & tstart=tstart-tzero & tstop=tstop-tzero

;	sort into ascending order
oo=sort(tstart) & tstart=tstart(oo)
oo=sort(tstop) & tstop=tstop(oo)

;	get rid of overlaps
oo=where(tstop lt tstart,moo)	;meaning previous interval has overextened
while moo gt 0 do begin		;{loop over them until no overlaps left
  ok=where(tstop ge tstart,mok)		;these are ok
  if mok gt 0 then begin
    tstart=tstart(ok) & tstop=tstop(ok)	;reset by deleting intermediates
  endif else begin
    tstart=[tsrmn] & tstop=tstart & return	;nothing left?!
  endelse
  oo=where(tstop lt tstart,moo)		;continue to check
endwhile			;MOO>0}
nstr=n_elements(tstart) & nstp=n_elements(tstp)

;	get the interval jumps
delt=tstart(1:*)-tstop
oj=where(delt le jump,moj)

;	merge the contiguous intervals
ij=intarr(nstr-1)+1 & if moj gt 0 then ij(oj)=0
ok=where(ij gt 0,mok)		;breaks in TIs are here
if mok gt 0 then begin		;(at least one gap exists
  tstr=tstart(ok+1) & tstp=tstop(ok)
  tstart=[0.,tstr] & tstop=[tstp,tspmx-tzero]
endif else begin		;)(no gaps at all
  tstart=[tsrmn-tzero] & tstop=[tspmx-tzero]
endelse				;MOK=0)

;	put back the offset
tstart=tstart+tzero
tstop=tstop+tzero

return
end
