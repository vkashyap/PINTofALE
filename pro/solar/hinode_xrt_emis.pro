function hinode_xrt_emis,filter,tgrid,ldbdir=ldbdir,cdbdir=cdbdir,$
	ioneqf=ioneqf,abund=abund,cieroot=cieroot,n_e=n_e,NH=NH,$
	xrteff=xrteff,lstr=lstr,cstr=cstr,toel=toel,toph=toph,toDN=toDN,$
	EM0=EM0,pixsize=pixsize,xrtl=xrtl,xrtc=xrtc,verbose=verbose,$
	contam_time=contam_time, _extra=e
;+
;function	hinode_xrt_emis
;	compute and return the combined emissivities of Hinode/XRT
;	filters in an array of form [Ntemp,Nfilt] in units of
;	[1d-23 ergs cm^5/s] for each pixel
;
;syntax
;	xrtemis=hinode_xrt_emis(filter,tgrid,ldbdir=ldbdir,$
;	cdbdir=cdbdir,ioneqf=ioneqf,abund=abund,cieroot=cieroot,$
;	n_e=n_e,NH=NH,xrteff=xrteff,lstr=lstr,cstr=cstr,$
;	/toel,/toph,/toDN,EM0=EM0,pixsize=pixsize,verbose=verbose,$
;	contam_time=contam_time,chidir=chidir,fH2=fH2,He1=He1,HeII=HeII)
;
;parameters
;	filter	[INPUT; required] scalar or array of filter names
;		* for Hinode/XRT, acceptable filter names are
;		  Al-mesh, Al-poly, C-poly, Ti-poly, Be-thin,
;		  Be-med, Al-med, Al-thick, Be-thick,
;		  Al-poly/Al-mesh, Al-poly/Ti-poly, Al-poly/Al-thick,
;		  Al-poly/Be-thick, C-poly/Ti-poly
;		* set to 'all' to get everything
;	tgrid	[I/O] the log(T[K]) grid over which the output
;		is to be defined
;		* if not specified on input, will be set to be
;		  the same as that in the emissivity tables,
;		  the PoA default, findgen(81)*0.05+4
;
;keywords
;	ldbdir	[INPUT; '$CHIANTI'] line emissivity database directory
;	cdbdir	[INPUT; '$CONT'] continuum emissivity database directory
;	ioneqf	[INPUT; 'ioneq/mazzotta_etal.ioneq'] location of the
;		ion fraction tables
;	abund	[INPUT; getabund('Grevesse et al')] abundances
;	cieroot	[INPUT; 'cie'] root name for the files in $CONT
;	n_e	[INPUT; 1e9] electron number density in the plasma
;	NH	[INPUT; 0] H column density to apply
;	xrteff	[I/O] the XRT effective areas; input can be any of the
;		following:-
;		- a structure that contains the following fields --
;		  {TYPE, CHANNEL_NAME, WAVE, EFF_AREA, LENGTH}
;		  which can be generated using the SSW routine
;		  xrteff = MAKE_XRT_WAVE_RESP(contam_time=contam_time)
;		- the name of a save file that contains the previously
;		  generated structure, named either XRTEFF or EFF
;		- a flag that indicates that this must be calculated
;		  in situ by calling MAKE_XRT_WAVE_RESP() (requires SSW)
;		  - this can be accomplished by setting /XRTEFF, but
;		    that will prevent it from being returned up.  instead,
;		    it is better to do something like
;		    xrteff=1 & xrtemis=hinode_xrt_emis(...,xrteff=xrteff,...)
;		* on output, will always contain the structure that
;		  was read in or calculated
;	lstr	[OUTPUT] line emissivities structure, as read in
;		from RD_LINE()
;	cstr	[OUTPUT] continuum emissivities structure, as read in
;		from RD_CONT()
;	toel	[INPUT; default=0] if set, returns the output in units
;		of [el cm^5/s]
;	toph	[INPUT; default=0] if set, returns the output in units
;		of [ph cm^5/s]
;	toDN	[INPUT; default=0] if set, returns the output in units
;		of [DN cm^5/pix]
;		* NOTE: TODN overrides TOPH overrides TOEL
;	EM0	[INPUT; default=1] if set, multiples the output by
;		this value to derive [ph/s] or [ph/cm^2/s] depending
;		on whether the input has units [cm^-3] or [cm^-5]
;	pixsize	[INPUT; default=1.0024 arcsec] the size of a pixel on
;		the detector
;	xrtl	[OUTPUT] the response due solely to lines
;	xrtc	[OUTPUT] the response due solely to continuum
;		* both XRTL and XRTC are identical in size and units
;		  to the primary output, XRTEMIS
;	verbose	[INPUT; default=0] controls chatter
;	contam_time	[INPUT] passed straight to MAKE_XRT_WAVE_RESP()
;	_extra	[INPUT ONLY] pass defined keywords to subroutines
;		FOLD_IONEQ : CHIDIR
;		ISMTAU : FH2, HE1, HEII
;		MAKE_XRT_WAVE_RESP : INDEX, CONTAM_THICK, CONTAM_TIME
;
;history
;	vinay kashyap (Jun2007)
;	added keywords TOEL, PIXSIZE, XRTL, XRTC (VK; Mar2008)
;	changed call from CALC_XRT_EFFAREA() to MAKE_XRT_WAVE_RESP() (VK; Apr2009)
;	bugfix: FILTER was failing to match to XRTEFF.NAME when specified
;	  (VK; Nov2009)
;-

forward_function make_xrt_wave_resp

;	usage
ok='ok' & np=n_params() & nf=n_elements(filter) & szf=size(filter,/type)
if np eq 0 then ok='Insufficient parameters' else $
 if nf eq 0 then ok='FILTER is undefined' else $
  if szf ne 7 then ok='FILTER must be a char string'
if ok ne 'ok' then begin
  print,'Usage: xrtemis=hinode_xrt_emis(filter,tgrid,ldbdir=ldbdir,$'
  print,'       cdbdir=cdbdir,ioneqf=ioneqf,abund=abund,cieroot=cieroot,$'
  print,'       n_e=n_e,NH=NH,xrteff=xrteff,lstr=lstr,cstr=cstr,$'
  print,'	/toel,/toph,/toDN,EM0=EM0,pixsize=pixsize,verbose=verbose,$'
  print,'	contam_time=contam_time,$'
  print,'       chidir=chidir,fH2=fH2,He1=He1,HeII=HeII)'
  print,'  compute and return emissivity functions for Hinode/XRT filters'
  if np gt 0 then begin
    message,ok,/informational
    print,'Al-mesh, Al-poly, C-poly, Ti-poly, Be-thin, Be-med, Al-med, '+$
	  'Al-thick, Be-thick, Al-poly/Al-mesh, Al-poly/Ti-poly, '+$
	  'Al-poly/Al-thick, Al-poly/Be-thick, C-poly/Ti-poly'
  endif
  return,-1L
endif

;	keywords
vv=0L & if keyword_set(verbose) then vv=long(verbose[0])>1
if not keyword_set(ldbdir) then ldbdir='$CHIANTI'
if not keyword_set(ioneqf) then ioneqf='ioneq/mazzotta_etal.ioneq'
if n_elements(abund) lt 30 then abund=getabund('grevesse et al.')
if not keyword_set(cdbdir) then cdbdir='$CONT'
if not keyword_set(cieroot) then cieroot='cie'
if not keyword_set(n_e) then n_e=1e9
if not keyword_set(NH) then NH=0
if not keyword_set(EM0) then EM0=1.D
pixsizx=1.0024 & pixsizy=pixsizx
if keyword_set(pixsize) then begin
  npix=n_elements(pixsize) & szp=size(pixsize,/type)
  if szp[0] le 5 then begin
    if pixsize[0] ne 0 then pixsizx=abs(pixsize[0])
    if npix ge 2 then $
     if pixsize[1] ne 0 then pixsizy=abs(pixsize[1]) else pixsizy=pixsizx
  endif else $
    message,'ignoring input PIXSIZE; using default',/informational
endif
if not keyword_set(xrteff) then begin	;(xrteff is not set
  ; look for save file in !ARDB
  defsysv,'!ARDB',ardbdir,exists=ivar	;!ARDB is defined if ivar=1
  if ivar eq 0 then begin	;(!ARDB undefined
    help,/source_files,output=cc
    ii=strpos(cc,'hinode_xrt_emis.pro') & oo=where(ii ge 0,moo)
    if moo eq 0 then begin	;(this file cannot be found?
      ardbdir='/data/fubar/SCAR/ardb'
    endif else begin		;file not found)(file found
      ccc=cc[oo[0]]
      iii=strpos(ccc,'pro') & i0=strpos(ccc,' ',/reverse_search)
      ardbdir=filepath('',root_dir=strmid(ccc,i0+1,iii-i0-1),subdir='ardb')
    endelse			;file found)
  endif else setsysval,'!ARDB',ardbdir,/getval	;!ARDB undefined)
  message,'restoring XRTEFF from '+filepath('xrt_eff.sav',root_dir=ardbdir),$
  	/informational
  restore,filepath('xrt_eff.sav',root_dir=ardbdir),verbose=(vv<1)
  if n_tags(xrteff) eq 0 then begin	;(xrteff doesn't exist
    if n_tags(eff) eq 0 then begin	;(eff also doesn't exist
      message,'XRT effective areas not found; quitting',/informational
      return,-1L
    endif else xrteff=eff		;eff)
  endif					;xrteff)
endif else begin			;xrteff is not set)(is set
  if n_tags(xrteff) eq 0 then begin	;(xrteff not a structure
    if size(xrteff,/type) eq 7 then begin	;(xrteff is a filename
      message,'restoring XRTEFF from '+xrteff[0],/informational
      restore,xrteff[0],/verbose
      if n_tags(xrteff) eq 0 then begin	;(xrteff doesn't exist
        if n_tags(eff) eq 0 then begin	;(eff also doesn't exist
          message,'XRT effective areas not found; quitting',/informational
          return,-1L
        endif else xrteff=eff		;eff)
      endif				;xrteff)
    endif else begin			;xrteff not filename)(calculate
      message,'calling make_xrt_wave_resp() to compute XRTEFF',$
      	/informational
      ;xrteff=calc_xrt_effarea()	;REQUIRES SSW!
      wav=make_xrt_wave_resp(contam_time=contam_time)	;REQUIRES SSW!
      xrteff=wav.EFFAR
      message,'calling make_xrt_wave_resp() to compute XRTEFF',$
      	/informational
    endelse				;calculate xrteff in situ)
  endif 				;xrteff not a structure)
endelse					;xrteff)
intoel=0 & if keyword_set(toel) then intoel=1
intoph=0 & if keyword_set(toph) then intoph=1
if keyword_set(intoph) then intoel=0	;override TOEL
intoDN=0 & if keyword_set(toDN) then intoDN=1
if keyword_set(intoDN) then intoph=0	;override TOPH
if keyword_set(intoDN) then intoel=0	;override TOEL

;	extract some useful arrays from XRTEFF
;
;<-- someday: here put in a check to make sure XRTEFF is kosher -->
;
WAVE=xrteff.WAVE
;EFF_AREA=xrteff.TRANS
;CHANNEL_NAME=xrteff.CHANNEL_NAME & nchan=n_elements(CHANNEL_NAME)
EFF_AREA=xrteff.EFF_AREA
CHANNEL_NAME=xrteff.NAME & nchan=n_elements(CHANNEL_NAME)
for i=0L,nchan-1L do channel_name[i]=(strsplit((xrteff.NAME)[i],';',/extract))[0]
LENGTH=xrteff.LENGTH & lmin=min(length,max=lmax)
wmin=min(WAVE[0L:lmin-1L,*]) & wmax=max(WAVE[0L:lmax-1L,*])
ifilt=intarr(nchan)-1
if strpos(strlowcase(filter[0]),'all') ge 0 or $
   strtrim(filter[0],2) eq '' then ifilt=indgen(nchan) else begin
  ;for i=0L,nchan-1L do for j=0L,nf-1L do ifilt[i]=strpos(strupcase(filter[j]),strupcase(CHANNEL_NAME[i]))
  for i=0L,nchan-1L do for j=0L,nf-1L do begin
    oo=where(strupcase(filter[j]) eq strupcase(CHANNEL_NAME[i]),moo)
    if moo gt 0 then ifilt[i]=i
  endfor
endelse
ofilt=where(ifilt ge 0,mofilt)
if mofilt eq 0 then begin
  message,'No filters selected; quitting',/informational
  return,-1L
endif

;	some other useful stuff
if vv ge 5 then peasecolr
inicon,fundae=fundae
fov=pixsizx*pixsizy
fov=fov*((!pi/180.)*(1./3600.))^2/(4.*!dpi)	;FOV of 1" pixels covers this many steradians/pix
ergperDN=fundae.ergev*3.65*60.	;[erg/eV]*[eV/e]*[e/DN]
if not keyword_set(toDN) then ergperDN=1.	;i.e., NO conversion to DN, even if
						;"ergperDN" is used to downstream!
if keyword_set(intoel) then ergperDN=fundae.ergev*3.65	;[erg/eV]*[eV/e]
ergperDN=1d23*ergperDN		;to account for the 1e-23 in the emissivities
ergperDN=ergperDN/fov
xtitle='log T' & ytitle='[10!u-23!n ergs cm!u5!n s!u-1!n]'
if keyword_set(intoel) then ytitle='[electrons cm!u5!n s!u-1!n]'
if keyword_set(intoph) then ytitle='[ph cm!u5!n s!u-1!n]'
if keyword_set(intoDN) then ytitle='[DN cm!u5!n s!u-1!n]'

;	read in line emissivities
if n_tags(lstr) eq 0 then begin
  lconf=rd_line(atom,n_e=n_e,wrange=[wmin,wmax],dbdir=ldbdir,fstr=lstr,verbose=verbose)
  ;	if APED, then remove the Anders & Grevesse abundances
  if strpos(strlowcase(ldbdir),'ape') ge 0 then begin
    apedance,lconf,lstr.Z & lstr.LINE_INT=lconf
  endif
  ;	if not APED, apply ion balances
  if strpos(strlowcase(ldbdir),'ape') lt 0 and strpos(strlowcase(ioneqf),'none') lt 0 then begin
    lconf=fold_ioneq(lstr.LINE_INT,lstr.Z,lstr.JON,eqfile=ioneqf,verbose=verbose, _extra=e)
    lstr.LINE_INT=lconf
  endif
endif
lwvl=abs(lstr.WVL) & Z=lstr.Z & nlw=n_elements(Z)
if keyword_set(toph) then nrgl=(fundae.h*fundae.c*1e8)/lwvl	;[erg/ph]

;	read in continuum emissivities
if n_tags(cstr) eq 0 then begin
  cconf=rd_cont(cieroot,n_e=n_e,wrange=[wmin,wmax],dbdir=cdbdir,abund=abund,fcstr=cstr,verbose=verbose,_extra=e)
endif
cww=mid2bound(cstr.midWVL) & cwvl=0.5*(cww[1:*]+cww) & ncw=n_elements(cwvl)
if keyword_set(toph) then nrgc=(fundae.h*fundae.c*1e8)/cwvl	;[erg/ph]
cdw=abs(cww[1:*]-cww)
cconf=cstr.CONT_INT & for i=0L,n_elements(cstr.LOGT)-1L do cconf[i,*]=cconf[i,*]*cdw
;	recast the continuum on the same temperature grid as the line emissivities
if n_elements(cstr.LOGT) ne n_elements(lstr.LOGT) then begin
  tmp=rebinx(cconf,cstr.LOGT,lstr.LOGT,verbose=vv)
  if vv gt 9 then stop,'HALTing; type .CON to continue'
  cconf=tmp
endif

;	any intrinsic absorption?
if keyword_set(NH) then begin
  lopt=exp(-ismtau(lwvl,NH=NH,/bam,abund=abund,verbose=(vv<1),_extra=e))
  copt=exp(-ismtau(cwvl,NH=NH,/bam,abund=abund,verbose=(vv<1),_extra=e))
endif else begin
  lopt=fltarr(nlw)+1.
  copt=fltarr(ncw)+1.
endelse
if keyword_set(toph) then lopt=lopt*(1d-23/nrgl)
if keyword_set(toph) then copt=copt*(1d-23/nrgc)

;	output
if n_elements(tgrid) eq 0 then tgrid=lstr.LOGT
nT=n_elements(tgrid) & nmaxT=n_elements(lstr.LOGT) & jT=findex(lstr.LOGT,tgrid)
nf=mofilt & nW=lmin
xrtemis=dblarr(nT,nf) & xrtl=dblarr(nT,nf) & xrtc=dblarr(nT,nf)
totemis=dblarr(nT)
totlemis=dblarr(nT)
totcemis=dblarr(nT)
warf=WAVE[0:nW-1L,0] & larf=fltarr(nlw,nf) & carf=fltarr(ncw,nf)
for i=0,nf-1 do larf[*,i]=(interpol(((EFF_AREA)[0:nw-1L,ofilt[i]]),warf,lwvl)>0)<(max(EFF_AREA))
for i=0,nf-1 do carf[*,i]=(interpol(((EFF_AREA)[0:nw-1L,ofilt[i]]),warf,cwvl)>0)<(max(EFF_AREA))

;	compute the response
zab=abund[Z-1]
for it=0L,nT-1L do begin
  if vv gt 0 then kilroy
  jj=fix(jT[iT])
  jj0=(jj>0)<(nmaxT-1L) & jj1=(jj+1)<(nmaxT-1L)	;no extrapolations!
  djT=abs(jT[iT]-jj)
  ;
  ;	intensities
  linint=lopt*reform((lstr.LINE_INT)[jj0,*])*EM0*zab
  conint=copt*reform((cconf)[jj0,*])*EM0
  if djT gt 0.001 and $
     jj ge 0 and jj lt nmaxT-2 then begin	;(interpolate between slices
    linint=djT*linint + (1.-djT)*lopt*reform((lstr.LINE_INT)[jj1,*])*EM0*zab
    conint=djT*conint + (1.-djT)*copt*reform((cconf)[jj1,*])*EM0
  endif						;JT>0.001)
  ;	fluxes
  totemis[it]=(total(linint)+total(conint))/ergperDN
  totlemis[it]=(total(linint)+0*total(conint))/ergperDN
  totcemis[it]=(0*total(linint)+total(conint))/ergperDN
  for jf=0,nf-1 do xrtemis[it,jf]=$
	(total(linint*larf[*,jf])+total(conint*carf[*,jf]))/ergperDN
  for jf=0,nf-1 do xrtl[it,jf]=(total(linint*larf[*,jf]))/ergperDN
  for jf=0,nf-1 do xrtc[it,jf]=(total(conint*carf[*,jf]))/ergperDN
  if it gt 0 and vv ge 5 then begin
    plot,tgrid[0:it],totemis[0:it,0],/yl,yr=max(xrtemis[0:it,*])*[1e-5,2],$
	xtitle=xtitle,ytitle=ytitle,_extra=e
    for i=0,nf-1 do oplot,tgrid,xrtemis[*,i],col=1+(i mod 8)
    if vv ge 10 then begin
      dy=(!y.crange[1]-!y.crange[0])/(n_elements(CHANNEL_NAME)+5)
      for i=0,nf-1 do xyouts,!x.crange[0]+0.1,10.^(!y.crange[1]-(i+1)*dy),CHANNEL_NAME[ofilt[i]],color=1+(i mod 8)
    endif
  endif
endfor

if vv ge 5 then begin
  plot,tgrid,totemis,/yl,yr=max(xrtemis)*[1e-5,2],$
    xtit=xtitle,ytit=ytitle,_extra=e
  for i=0,nf-1 do oplot,tgrid,xrtemis[*,i],col=1+(i mod 8)
  dy=(!y.crange[1]-!y.crange[0])/(n_elements(CHANNEL_NAME)+5)
  for i=0,nf-1 do xyouts,!x.crange[0]+0.1,10.^(!y.crange[1]-(i+1)*dy),$
	CHANNEL_NAME[ofilt[i]],color=1+(i mod 8)
endif

if vv gt 1000 then stop,'HALTing; type .CON to continue'

return,xrtemis
end
