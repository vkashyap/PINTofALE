pro exaped,filroot,emis,tlog,wvl,edens,Z,ion,desig,econf,jon,src,$
	llist=llist,coco=coco,tlist=tlist,noerg=noerg,okeV=okeV,$
	atomdb=atomdb,logT=logT,verbose=verbose, _extra=e
;+
;procedure	exaped
;	EX APED AD IDL
;	extract line and/or continuum emissivity data from an APED file
;
;syntax
;	exaped,filroot,emis,tlog,wvl,edens,Z,ion,desig,econf,jon,src,$
;	/llist,/coco,/tlist,/noerg,/okeV,atomdb=atomdb,verbose=verbose
;
;parameters
;	filroot	[INPUT; required] APED fits file root
;		* if this includes a ".fits", the suffix is assumed
;		  to be already included and nothing extra is added.
;		* the suffix '_linelist.fits' is automatically applied,
;		  so don't specify that part.
;		  (this is inflexible at this point because we only
;		  support reading in this file.  sometime later on
;		  the other APED format will also be supported.)
;		* looks for ".gz" automatically if the filename
;		  specified without it does not exist
;	emis	[OUTPUT; required] emissivity array, of the format
;		(NTLOG,NWVL,NEDENS)
;		* APED units are [ph cm^3/s], which will be converted
;		  into [1e-23 ergs cm^3/s] unless the NOERG keyword
;		  is set
;	tlog	[OUTPUT; required] log10(Temperature [K]) grid on
;		which EMIS are placed
;	wvl	[OUTPUT; required] wavelength array
;		* units are in [Ang] unless the OKEV keyword is set,
;		  in which case output in [keV]
;	edens	[OUTPUT; required] electron densities at which
;		EMIS have been calculated
;	Z	[OUTPUT; required] atomic numbers of species producing
;		the given line
;	ion	[OUTPUT; required] ionic state of species corresponding
;		to the line
;	desig	[OUTPUT] (2,NWVL) string array of lower and upper level
;		designation
;		* currently, just the lower and upper levels
;	econf	[OUTPUT] (2,NWVL) string array of electron configuration
;		of lower and upper levels of the transition in question
;		* not implemented yet
;	jon	[OUTPUT] the ionic state that matters to the ion balance
;		* not implemented yet
;	src	[OUTPUT] numeric code specifiying that this is an APED
;		output, set to 3
;
;keywords
;	llist	[INPUT] if set, reads from FILROOT_linelist.fits
;		* currently this is the only option supported
;	coco	[INPUT] if set, reads from FILROOT_coco.fits
;		* NOT YET IMPLEMENTED
;	tlist	[INPUT] if set, reads from FILROOT_line.fits
;		* NOT YET IMPLEMENTED
;	noerg	[INPUT] if set, EMIS will be output in same units as
;		in the APED files, i.e., [ph cm^3/s]
;	okeV	[INPUT] if set, WVL will be converted from [Ang] to [keV]
;	atomdb	[INPUT] directory containing input file
;		* default is !ATOMDB
;		* hardcoded default is /data/atomdb/
;		* prepended to FILROOT only if FILROOT does not begin with
;		  a '/' (UNIX), ':' (MacOS), '\' (Windows), or '[' (VMS?)
;	logT	[INPUT] if defined as an array, then interpolates/extrapolates
;		EMIS(NTLOG,NWVL) to a new grid EMIS(NLOGT,NWVL) prior to
;		output.  TLOG will be overwritten by LOGT.
;	verbose	[INPUT] controls chatter
;	_extra	[JUNK] ignore -- here only to prevent crashing the program
;
;subroutines
;	SETSYSVAL
;	ARRAYEQ
;	REBINX
;
;description
;	The ATOMDB file format is described in
;		http://cxc.harvard.edu/atomdb/features_formats.html
;	In this case, the line list (grouped by line) file is read, and
;	in particular, only the block EMISSIVITY is read.
;
;history
;	vinay kashyap (Jul02)
;	bug correction when multiple densities are included (VK; Mar03)
;	for v2+, always go through the multiple densities code and ignore
;	  the short cut (VK; May11)
;-

;	usage
ok='ok' & np=n_params() & nf=n_elements(filroot)
szf=size(filroot) & nszf=n_elements(szf)
if np lt 7 then ok='Insufficient parameters' else $
 if nf eq 0 then ok='FILROOT is undefined' else $
  if nf gt 1 then ok='FILROOT must be a scalar' else $
   if szf[nszf-2] ne 7 then ok='FILROOT must be a string'
if ok ne 'ok' then begin
  print,'Usage: exaped,filroot,emis,tlog,wvl,edens,Z,ion,desig,econf,jon,src,$'
  print,'       /llist,/coco,/tlist,/noerg,/okeV,atomdb=atomdb,logT=logT,$'
  print,'       verbose=verbose'
  print,'  EX APED AD IDL: extract emissivities from APED file'
  if np ne 0 then message,ok,/info
  return
endif

;	keywords
filsufx='_linelist.fits'	;by default, read from _linelist.fits
if keyword_set(tlist) then filsufx='_line.fits'
if keyword_set(coco) then filsufx='_coco.fits'
if keyword_set(llist) then filsufx='_linelist.fits'
if not keyword_set(noerg) then noerg=0
if not keyword_set(okev) then okev=0
vv=0 & if keyword_set(verbose) then vv=long(verbose[0])>1
apeddir='/data/atomdb/'
ivar=0 & defsysv,'!ATOMDB',exists=ivar	;if !ATOMDB exists
if ivar ne 0 then setsysval,'ATOMDB',apeddir,/getval
if keyword_set(atomdb) then apeddir=strtrim(atomdb,2)

;	check input
infil=filroot
if filsufx eq '_line.fits' then begin
  message,filsufx+': Sorry, this option has not yet been implemented',/info
  return
endif
if filsufx eq '_coco.fits' then begin
  message,filsufx+': Sorry, this option has not yet been implemented',/info
  return
endif
if strpos(filroot,'.fits',0) ge 0 then begin
  if strpos(filroot,filsufx,0) lt 0 then begin
    message,filroot+': very likely not in acceptable format; exiting',/info
    return
  endif
endif else infil=filroot+filsufx
c=strmid(infil,0,1)
if c ne '/' and c ne ':' and c ne '\' and c ne '[' then $
	infil=filepath(infil,root_dir=apeddir)
;
fil=findfile(infil,count=ifil)
if ifil eq 0 then begin
  if vv gt 2 then message,infil+': not found; looking for '+infil+'.gz',/info
  fil=findfile(infil+'.gz',count=ifil)
  if ifil gt 0 then begin
    ;infil=infil+'.gz'	;this is not necessary for headfits() and mrdfits()?
    if vv gt 0 then message,'Reading from '+infil,/info
  endif else begin
    message,infil+': not found; exiting',/info
    return
  endelse
endif else if vv gt 3 then message,'Reading from '+infil,/info

if filsufx eq '_linelist.fits' then begin	;(read _linelist.fits

  ;	figure out which extension and read data
  iext=0 & extnam=''
  while extnam ne 'EMISSIVITY' do begin
    iext=iext+1 & hdr=headfits(infil,exten=iext)
    if strtrim(hdr[0],2) eq '-1' then begin
      message,infile+': EMISSIVITY block not found; exiting',/info
      return
    endif
    extnam=sxpar(hdr,'EXTNAME')
  endwhile
  t=mrdfits(infil,iext)
  if vv gt 10 then help,t,/struct

  ;	at which temperatures?
  ttemp=t.TEMPERATURE
  log_temp=sxpar(hdr,'LOG_TEMP')
  num_temp=sxpar(hdr,'NUM_TEMP')
  tempstep=sxpar(hdr,'TEMPSTEP')
  tempstrt=sxpar(hdr,'TEMPSTRT')
  tlog=findgen(num_temp)*tempstep+tempstrt
  if not keyword_set(log_temp) then $
	tlog=alog10(findgen(num_temp)*tempstep+tempstrt)
  ntlog=n_elements(tlog)

  ;	at which wavelengths?
  wvl=t.Lambda & nwvl=n_elements(wvl)

  ;	at which densities?
  tdens=t.DENSITY & nmat=n_elements(tdens)
  edens=reform(tdens[*,0])
  edens=edens[uniq(edens,sort(edens))] & idens=bytarr(nmat)
  for i=0L,n_elements(edens)-1L do begin
    oo=where(tdens eq edens[i],moo) & if moo gt 0 then idens[oo]=1
    if vv gt 10 then kilroy,dot='n_e='+strtrim(edens[i],2)+' '
  endfor
  oo=where(idens eq 0,moo)
  while moo gt 0 do begin
    edens=[edens,tdens[oo[0]]] & idens[oo]=1
    oo=where(idens eq 0,moo)
    if vv gt 10 then kilroy,dot='n_e='+strtrim(tdens[oo[0]],2)+' '
  endwhile

  ok=where(edens ne 0,nedens)
  if nedens eq 0 then begin
    message,'Cannot understand DENSITY column; exiting',/info
    return
  endif else edens=edens[ok]

  ;	from which elements?
  Z=t.Element & nZ=n_elements(Z)

  ;	what ionic states?
  ion=t.Ion & nIon=n_elements(Ion)

  ;	levels
  desig=strarr(2,nwvl)
  desig[0,*]=strtrim(t.LowerLev,2)
  desig[1,*]=strtrim(t.UpperLev,2)

  ;	e configurations
  econf=strarr(2,nwvl)		;hmm.  nothing to put here yet.

  ;	which ionic states, really?
  jon=ion			;see below for correction due to DR lines

  ;	data source
  src=intarr(nwvl)+3

  ;	oh, and the emissivities
  temis=t.EMISSIVITY
  emis=dblarr(ntlog,nwvl,nedens)
  n_chan=t.arraysize & f_chan=intarr(nwvl) & t_chan=(t.Temperature)[0,*]
  if max(t_chan) gt 99 then t_chan=alog10(t_chan)

  if nedens eq 1 and strpos(atomdb,'v1.') ge 0 then begin	;(this is the fast version

    for i=0L,nwvl-1L do begin
      dT=min(abs(tlog-t_chan[i]),j) & f_chan[i]=j
      if dT gt tempstep then message,'BUG!'
    endfor
    i0=0*n_chan & i1=n_chan
    for iw=0L,nwvl-1 do $
	emis[f_chan[iw]:f_chan[iw]+n_chan[iw]-1L,iw,0]=temis[i0[iw]:i1[iw]-1L,iw]

  endif else begin		;)(painstaking, step by step

    for iw=0L,nwvl-1L do begin		;{for each wavelength
      ee=reform(temis[*,iw])
      dd=reform(alog10(tdens[*,iw])) & tt=reform(alog10(ttemp[*,iw]))
      for id=0L,nedens-1L do begin	;{for each density
        oid=where(abs(dd-alog10(edens[id])) lt 1e-3,moid)
        if moid gt 0 then begin
          for it=0L,ntlog-1L do begin	;{for each temperature
    	    oit=where(abs(tt[oid]-tlog[it]) lt 1e-3*tempstep,moit)
    	    if moit gt 0 then emis[it,iw,id]=ee[oid[oit[0]]]
          endfor			;IT=0,NTLOG-1}
        endif else if vv gt 100 then kilroy,dot='x'
      endfor				;ID=0,NEDENS-1}
      if vv gt 10 then kilroy else $
       if vv gt 0 and iw eq 100L*long(iw/100L) then kilroy
      if vv-1000L eq iw then stop,'halting; type .CON to continue'
    endfor				;IW=0,NWVL-1}

  endelse			;NEDENS?1 for v1.3 and earlier)

  ;for i=0L,nwvl-1L do begin
    ;dT=min(abs(tlog-t_chan[i]),j) & f_chan[i]=j
    ;if dT gt tempstep then message,'BUG!'
  ;endfor
  ;for id=0,nedens-1 do begin
    ;i0=id*n_chan & i1=(id+1L)*n_chan
    ;for iw=0L,nwvl-1 do begin
      ;;emis[f_chan[iw]:f_chan[iw]+n_chan[iw]-1L,iw,id]=temis[i0[iw]:i1[iw]-1L,iw]
      ;j0=f_chan[iw]+id*nedens & j1=j0+n_chan[iw]/nedens-1L
      ;emis[j0:j1,iw,id]=temis[i0[iw]:i1[iw]-1L,iw]
    ;endfor
  ;endfor

  ;	DR lines have upper levels > 10000
  odr=where(t.UpperLev gt 10000,modr)
  if modr gt 0 then begin
    desig[0,odr]=desig[0,odr]+' (DR)'
    jon[odr]=jon[odr]+1
  endif

endif						;read _linelist.fits)

if filsufx eq '_coco.fits' then begin		;(read _coco.fits
  message,'not implemented yet, sorry',/info
endif						;read _coco.fits)

if filsufx eq '_line.fits' then begin		;(read _line.fits
  message,'not implemented yet, sorry',/info
endif						;read _line.fits)

;	convert from [ph cm^3/s] to [1e-23 ergs cm^3/s]
h=6.626176e-27 & c=2.9979e10	;[ergs s] & [cm/s]
nrg=h*c*1e8/abs(wvl)		;[ergs/ph]
if not keyword_set(noerg) then ph2erg=nrg*1e23 else $	;[1e-23 ergs/ph]
	ph2erg=1.
for i=0L,ntlog-1L do for k=0L,nedens-1 do emis[i,*,k]=emis[i,*,k]*ph2erg[*]

;	regrid, if necessary
nlogT=n_elements(logT)
if nlogT gt 1 then begin
  ieq=arrayeq(Tlog,logT)
  if ieq ne 1 then begin
    if vv gt 0 then message,'Regridding the emissivities to match LOGT',/info
    tmp=dblarr(nlogT,nwvl,nedens)
    for i=0,nedens-1 do begin
      if vv gt 1 then message,'EDENS='+strtrim(edens[i],2),/info
      tmp1=reform(emis[*,*,i])
      tmp2=rebinx(tmp1,Tlog,logT,verbose=vv)
      tmp[*,*,i]=tmp2[*,*]
    endfor
    emis=tmp & Tlog=logT
  endif
endif else begin
  if nlogT eq 1 and vv gt 1 then begin
    if logT[0] gt 0 then message,$
      'LOGT must have at least 2 elements; ignoring',/info else $
      message,'LOGT=0? cannot understand',/info
  endif
endelse

if keyword_set(okeV) then wvl=nrg/1.6021892e-9	;[ergs]/[ergs/keV]

return
end
