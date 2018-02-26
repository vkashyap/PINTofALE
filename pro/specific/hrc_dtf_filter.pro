function hrc_dtf_filter,dtffil,tstart,tstop,gtibeg=gtibeg,gtiend=gtiend,$
	dtfthr=dtfthr,gtifil=gtifil,verbose=verbose, _extra=e
;+
;function	hrc_dtf_filter
;	computes the average dead-time correction factor
;	for a new set of good time intervals obtained after
;	filtering the existing DTF file
;
;syntax
;	dtcor=hrc_dtf_filter(dtffil,tstart,tstop,gtibeg=gtibeg,gtiend=gtiend,$
;	dtfthr=dtfthr,gtifil=gtifil,verbose=verbose,$
;	leap=leap,fitsref=fitsref,exten=exten,fitshdr=fitshdr,form=form,page=page)
;
;parameters
;	dtffil	[INPUT; required] full path name to the dtf file
;	tstart	[OUTPUT] the new start times of the GTIs
;	tstop	[OUTPUT] the new stop times of the GTIs
;
;keywords
;	gtibeg	[INPUT] beginning times for any GTI that should also
;		be included
;	gtiend	[INPUT] ending times for GTIs corresponding to GTIBEG
;		* GTIBEG and GTIEND may be arrays
;		* the sizes of GTIBEG and GTIEND must match or else
;		  both are ignored
;	dtfthr	[INPUT] the threshold dead-time correction factor
;		below which to ignore time bins
;		* default=0.99
;	gtifil	[INPUT] full path name of the output file that will
;		contain the list of TSTART and TSTOP -- this can be
;		used as a filter in dmcopy to do the actual filtering
;		* if not set or the name contains "none", will not
;		  be written out
;		* if set to number, assumes that the file should be
;		  named 'hrc_dtf_filter_gti_GTIFIL.txt'
;		* if name ends in .fit, .fits, .ft, or .mt then
;		  writes out a FITS file, otherwise an ASCII file
;	verbose	[INPUT] controls chatter
;	_extra	[INPUT ONLY] pass defined keywords to subroutines
;		TI_CLEAN: LEAP
;		TI_WRITE: FITSREF, EXTEN, FITSHDR, FORM, PAGE
;
;subroutines
;	mrdfits
;	ti_clean
;	ti_write
;	ti_and
;	ti_filter
;
;history
;	vinay kashyap (Sep2006)
;-

;	usage
ok='ok' & np=n_params() & nf=n_elements(dtffil) & szf=size(dtffil,/type)
if np eq 0 then ok='Insufficient parameters' else $
 if nf eq 0 then ok='DTFFIL is undefined' else $
  if nf gt 1 then ok='DTFFIL must be a scalar string' else $
   if szf ne 7 then ok='DTFFIL must be the name of a file'
if ok ne 'ok' then begin
  print,'Usage: dtcor=hrc_dtf_filter(dtffil,tstart,tstop,$'
  print,'       dtfthr=dtfthr,gtifil=gtifil,verbose=verbose,$'
  print,'       leap=leap,fitsref=fitsref,exten=exten,fitshdr=fitshdr,$'
  print,'       form=form,page=page)'
  print,'  computes the average dead-time correction factor for a new'
  print,'  set of good time intervals for a Chandra/HRC observation'
  if np ne 0 then message,ok,/informational
  return,-1L
endif

;	keywords
vv=0L & if keyword_set(verbose) then vv=long(verbose[0])>1
;
dtfthresh=0.99 & if keyword_set(dtfthr) then dtfthresh=float(dtfthr[0])
;
outfil='STDOUT'
ngf=n_elements(gtifil) & szgf=size(gtifil,/type)
if ngf ne 0 then begin
  if szgf eq 7 then outfil=gtifil[0] else $
    outfil='hrc_dtf_filter_gti_'+strtrim(gtifil[0],2)+'.txt'
  if vv gt 0 then message,'GTIs will be written to file: '+outfil,/informational
endif
if strpos(strlowcase(outfil),'none') ge 0 then begin
  if vv gt 0 then message,'GTIs will not be written to file',/informational
  outfil='STDOUT'
endif

;	read in dtf file
if not keyword_set(dtfext) then dtfext=1
t=mrdfits(dtffil[0],dtfext,dtfhdr)
time=t.TIME & dtf=t.DTF
tzero=double(sxpar(dtfhdr,'TSTART')) > min(time)

;	GTIBEG and GTIEND
ok='ok' & ngb=n_elements(gtibeg) & nge=n_elements(gtiend)
if ngb eq 0 then ok='GTIBEG is undefined' else $
 if nge eq 0 then ok='GTIEND is undefined' else $
  if ngb ne nge then ok='GTIBEG and GTIEND do not match'
if vv gt 0 and (ngb ne 0 or nge ne 0) then message,ok,/informational
if ok eq 'ok' then begin
  gstart=gtibeg[*] & gstop=gtiend[*]
endif else gstart=min(time,max=gstop)

;	filter
ot=where(dtf ge dtfthresh,mot)
if mot eq 0 then begin
  message,'No times selected; exiting',/informational
  return,0.
endif

;	clean up GTIs and write to file
dt=time[1:*]-time & delta_t=median(dt) ;& dt=[dt,delta_t]
dstart=time[ot]-tzero & dstop=dstart+delta_t	;this assumes that the
						;dtf binning doesn't change
ti_clean,dstart,dstop, _extra=e
if vv gt 1 then begin
  if vv gt 2 then message,strtrim(mot,2)+' of '+$
	strtrim(n_elements(time),2)+' bins selected',/informational
  message,strtrim(n_elements(dstart),2)+' GTIs',/informational
endif
dstart=dstart+tzero
dstop=dstop+tzero
ti_and,dstart,dstop,gstart,gstop,tstart,tstop
ti_write,tstart,tstop,outfil, _extra=e

;	compute the average DTCOR
ot=ti_filter(dtf,tstart,tstop,count=mot)
if mot gt 0 then dtcor=mean(dtf[ot]) else dtcor=0.

return,dtcor
end
