function rd_3par,p,w,h,group,delp,pardir=pardir,root=root,missing=missing,$
	_extra=e
;+
;function 	rd_3par
;	read in parameters describing 3-parameter curve (e.g.,
;	gaussians) from parameter files.  returns the total
;	number of components read in.
;
;syntax
;	n=rd_3par(p,w,h,group,delp,pardir=pardir,root=root,missing=missing)
;
;parameters
;	p	[OUTPUT] positional parameter (e.g., mean)
;	w	[OUTPUT] width parameter (e.g., sigma)
;	h	[OUTPUT] height parameter (e.g., maximum)
;	group	[OUTPUT] grouping index
;	delp	[OUTPUT] offsets from 1st component in group
;		NOTE: these last two would be of use in generating models for
;		fitting (a given group moves around as an entity)!
;
;keywords
;	pardir	[INPUT; default: .] directory in which to look for
;		parameter files
;	root	[INPUT; default: gauss] prefix for parameter files
;	missing	[INPUT] 3-element array [p,w,h] containing defaults
;		to override hardcoded defaults to use in case any
;		of the values and/or files are missing
;	_extra	[JUNK] here only to prevent program from crashing
;
;usage summary
;	* call as a function
;	* returns number of components specified in parameter files
;	* all parameters are OUTPUT (so if any are filled prior to
;	  call, they will be overwritten!)
;			
;parameter file format
;	there are 3 files,
;	    pardir/root_pos.par (positions)
;	    pardir/root_wdt.par (widths)
;	    pardir/root_hgt.par (heights)
;	each file contains one component per line (exception -- see
;	below), and the components are linked by their line numbers.
;	in case there are more than 1 entry per line, the following
;	<sp> separated entries in _pos.par are taken to be delta_p from
;	the first value (and corresponding widths and heights in the
;	other files).  example: _pos.par may contain
;		200. 10. 20. 30.
;		215.
;		240.
;		180. -7. 21.
;	and _wdt.par may contain
;		5. 0 3.
;		
;		10.
;		5.
;	and missing values are either taken from the first column or
;	from the keyword MISSING, or from hard-coded default.
;	_hgt.par may be similar or even completely missing.
;	NOTE: the _pos.par file is assumed to be canonical, i.e., it
;	is assumed that *nothing* is missing.
;
;subroutines
;	NONE
;
;history
;	vinay kashyap (Oct 96)
;	added _EXTRA (VK; Oct98)
;-

np=n_params(0)
if np lt 3 then begin
  print, 'Usage: n=rd_3par(p,w,h,group,delp,pardir=pardir,root=root,missing=missing)'
  print, '  read in parameter files for 3-parameter curve models' & return,0
endif

;check keywords
if not keyword_set(root) then root='gauss' else root=strtrim(root,2)
if not keyword_set(missing) then begin
  missing=[0.,0.1,1.]
endif else begin
  nm=n_elements(missing)
  if nm eq 1 then missing=[missing,0.1,1.]
  if nm eq 2 then missing=[missing,1.]
endelse

;initialize
pfile=root+'_pos.par' & if keyword_set(pardir) then pfile=pardir+'/'+pfile
wfile=root+'_wdt.par' & if keyword_set(pardir) then wfile=pardir+'/'+wfile
hfile=root+'_hgt.par' & if keyword_set(pardir) then hfile=pardir+'/'+hfile
c1=findfile(pfile,count=fp)
c1=findfile(wfile,count=fw)
c1=findfile(hfile,count=fh)

;if _pos.par file does not exist, quit
if fp eq 0 then begin
  message,'no position parameter file available!',/info & return,0
endif else openr,upf,pfile,/get_lun

;read in positions of means
n=0				;number of entries
p=[-1.] & ni=[0] & group=[0] & delp=[0.] & irow=0
while not eof(upf) do begin
  c1='' & readf,upf,c1 & irow=irow+1
  pp=float(str_sep(strtrim(c1,2),' '))			;this line contains...
  np=n_elements(pp)					;...many columns
  p=[p,pp(0)] & for i=1,np-1 do p=[p,pp(0)+pp(i)]
  group=[group,intarr(np)+irow]
  delp=[delp,0.] & if np gt 1 then delp=[delp,pp(1:*)]
  n=n+np & ni=[ni,np]					;keeping track
endwhile
p=p(1:*) & ni=ni(1:*) & group=group(1:*) & delp=delp(1:*)
close,upf & free_lun,upf

;read in widths
w=0*p+missing(1) & irow=0 & i=0
if fw ne 0 then begin
  openr,uwf,wfile,/get_lun
  while not eof(uwf) do begin
    c1='' & readf,uwf,c1
    ww=float(str_sep(strtrim(c1,2),' '))		;read in columns
    nw=n_elements(ww)
    if strtrim(c1,2) eq '' then ww(0)=missing(1)	;missing line
    oo=where(ww eq 0.) & if oo(0) ne -1 then ww(oo)=ww(0)
    icol=ni(irow)
    if nw lt icol then begin				;missing columns
      w(i:i+nw-1)=ww(0:nw-1)
      w(i+nw:i+icol-1)=ww(0)
    endif else w(i:i+icol-1)=ww(0:icol-1)
    irow=irow+1 & i=i+icol
  endwhile
  close,uwf & free_lun,uwf
endif

;read in peaks (same as widths, above)
h=0*p+missing(2) & irow=0 & i=0
if fh ne 0 then begin
  openr,uhf,hfile,/get_lun
  while not eof(uhf) do begin
    c1='' & readf,uhf,c1
    hh=float(str_sep(strtrim(c1,2),' '))		;read in columns
    nh=n_elements(hh)
    if strtrim(c1,2) eq '' then hh(0)=missing(2)	;missing lines
    oo=where(oo eq 0.) & if oo(0) ne -1 then hh(oo)=hh(0)
    icol=ni(irow)
    if nh lt icol then begin				;missing columns
      h(i:i+nh-1)=hh(0:nh-1)
      h(i+nh:i+icol-1)=hh(0)
    endif else h(i:i+icol-1)=hh(0:icol-1)
    irow=irow+1 & i=i+icol
  endwhile
  close,uhf & free_lun,uhf
endif

return,n
end
