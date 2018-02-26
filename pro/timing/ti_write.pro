pro ti_write,tstart,tstop,tifile,fitsref=fitsref,exten=exten,$
	fitshdr=fitshdr,form=form,page=page, _extra=e
;+
;procedure	ti_write
;	writes out the start and stop times of a set of time intervals
;
;syntax
;	ti_write,tstart,tstop,tifile,fitsref=fitsref,exten=exten,$
;	form=form,fitshdr=fitshdr,page=page
;
;parameters
;	tstart	[INPUT; required] start times
;	tstop	[INPUT; required] stop times
;	tifile	[INPUT] name of file to write to (default: STDOUT)
;		* if name ends in .fit, .fits, .ft, or .mt then
;		  write out a FITS file with a bintable extension
;		  else an ascii file
;		* if filename contains 'none' or 'stdout', will print to STDOUT
;
;keywords
;	fitsref	[INPUT] name of existing FITS file from which to
;		steal the header
;	exten	[INPUT] extension of FITSREF from which steal header
;	fitshdr	[INPUT] explicitly specify the header.
;		* if FITSHDR is set, FITSREF and EXTEN will be ignored!
;	form	[INPUT] format in case of ASCII/STDOUT output
;		* default: '(d16.4,1x,d16.4)'
;		* ignored for FITS output
;	page	[INPUT] number of lines to display on the screen at a time
;		* default: 22
;		* ignored if not writing to STDOUT
;	_extra	[JUNK] here only to prevent crashing the program
;
;related
;	TI_CLEAN	- cleans up intervals
;	TI_OR		- merges intervals
;	TI_AND		- intersection of intervals
;	TI_COVER()	- computes coverage in given bin
;	TI_FILTER()	- filter an array by GTI
;
;history
;	vinay kashyap (May98)
;	added keyword _EXTRA (VK; Nov98)
;	bug correction with setting TIFILE='STDOUT' (VK; Sep06)
;	slightly more verbose prompt for PAGE (VK; Aug15)
;-

;	usage
nstr=n_elements(tstart) & nstp=n_elements(tstop)
if nstr eq 0 or nstp eq 0 then begin
  print,'Usage: ti_write,tstart,tstop,tifile,fitsref=fitsref,exten=exten,$'
  print,'       form=form,fitshdr=fitshdr,page=page'
  print,'  writes out the start and stop times of a set of time intervals'
  return
endif

;	stupid user tricks
if nstr ne nstp then begin
  message,'TSTART and TSTOP do not match; exiting',/info & return
endif
if n_elements(tifile) gt 1 then begin
  message,'confused.  too many output files! exiting',/info & return
endif

;	inputs
outyp='STDOUT'
if n_elements(tifile) eq 1 then begin
  i1=strpos(strlowcase(tifile),'.ft',0)
  i2=strpos(strlowcase(tifile),'.fit',0)
  i3=strpos(strlowcase(tifile),'.mt',0)
  if i1 ge 0 or i2 ge 0 or i3 ge 0 then outyp='FITS' else begin
    j1=strpos(strlowcase(tifile),'none',0)
    j2=strpos(strlowcase(tifile),'stdout',0)
    if j1 lt 0 and j2 lt 0 then outyp='ASCII'
  endelse
endif
;
if not keyword_set(form) then form='(d16.4,1x,d16.4)'

;	read/make FITS header
if outyp eq 'FITS' then begin
  if not keyword_set(fitshdr) then begin	;(no FITSHDR
    if not keyword_set(fitsref) then begin	;(no FITSREF
      mkhdr,hdr0,'',/extend
      fxbhmake,hdr1,nstr
      sxaddpar,hdr1,'TTYPE1','START','Interval start time'
      sxaddpar,hdr1,'TTYPE2','STOP','Interval end time'
      sxaddpar,hdr0,'NAXIS2',nstr,'The number of rows'
    endif else begin				;)(FITSREF
      tab=readfits(fitsref,hdr0,numrow=1)
      if keyword_set(exten) then begin		;(read in extension
        tab=readfits(fitsref,fitshdr,exten=exten,numrow=1)
      endif					;EXTEN)
      sxaddpar,hdr0,'NAXIS2',nstr,'The number of rows'
    endelse					;FITSREF)
  endif						;FITSHDR)
  sxaddpar,fitshdr,'NAXIS2',nstr,'The number of rows'
endif
      ;sxaddpar,fitshdr,'TDMIN1',min(tstart),'Smallest value present'
      ;sxaddpar,fitshdr,'TDMAX1',max(tstart),'Largest value present'
      ;sxaddpar,fitshdr,'TDMIN2',min(tstop),'Smallest value present'
      ;sxaddpar,fitshdr,'TDMAX2',max(tstop),'Largest value present'

;	output
if outyp eq 'FITS' then begin			;(FITS file
  writefits,tifile(0),0,hdr0	;write primary extension
  ;
  fxbcreate,uout,tifile(0),fitshdr	;(open binary table extension
  for i=0L,nstr-1 do fxbwrite,uout,tstart(i),1,i+1	;start times
  for i=0L,nstr-1 do fxbwrite,uout,tstop(i),2,i+1	;end times
  fxbfinish,uout			;close file)
endif else begin				;)(ASCII or STDOUT
  if outyp eq 'STDOUT' then begin
    if not keyword_set(page) then page=22
    openw,uout,'/dev/tty',/get_lun
    for i=0L,nstr-1L do begin
      printf,uout,i,tstop(i)-tstart(i),form='(i5,1x,d14.4,$)'
      printf,uout,tstart(i),tstop(i),form=form
      if i eq page*long((i-(page-1))/page)+(page-1) then begin
	;c1='' & kilroy,dot=c1 & c1=get_kbrd(1)
	c1='q/Q to return, x/z to stop, any key to continue: ' & kilroy,dot=c1 & c1=get_kbrd(1)
	if c1 eq 'q' or c1 eq 'Q' then return
	if c1 eq 'x' or c1 eq 'z' then stop
      endif
    endfor
  endif else begin
    openw,uout,tifile(0),/get_lun
    for i=0L,nstr-1L do printf,uout,tstart(i),tstop(i),form=form
  endelse
  close,uout & free_lun,uout
endelse						;OUTYP)

return
end
